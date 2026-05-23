"""Streamlit UI for optical coefficient plots."""

from __future__ import annotations

import sqlite3
import sys
from copy import deepcopy
from pathlib import Path

import pandas as pd
import streamlit as st

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.coefficients.calc.optics import calculate_jump, calculate_max_reflection
from src.coefficients.core.constants import COEFFICIENT_LIST, HC_CONST
from src.coefficients.core.models import Coefficients
from src.coefficients.plotting.charts import delta_beta_plot, plot_chart
from src.definitions import COMPOUNDS_DB_PATH, F1F2_DIR


BASE_COEFFICIENTS = COEFFICIENT_LIST[:-2]
DEFAULT_ENERGY_RANGE = (90.0, 700.0)
DEFAULT_SLIDER_ENERGY_MAX = 5000.0
MAX_ENERGY_LIMIT = 100000.0
MAX_COMPOUND_SELECTORS = 15
MAX_MATERIAL_PAIRS = 15


def clamp(value: float, min_value: float, max_value: float) -> float:
    return min(max(float(value), float(min_value)), float(max_value))


def is_complex_formula(formula: str) -> bool:
    """Return true for hydrate/adduct/grouped formulas that should be sorted later."""
    markers = ("(", ")", "∙", "·", "⋅", "*", "cdot")
    return any(marker in formula for marker in markers)


def sort_formulas(formulas: list[str]) -> list[str]:
    return sorted(formulas, key=lambda item: (is_complex_formula(item), item.lower()))


@st.cache_data
def load_compound_formulas() -> list[str]:
    with sqlite3.connect(COMPOUNDS_DB_PATH) as connection:
        rows = connection.execute(
            """
            SELECT DISTINCT formula
            FROM Properties
            WHERE formula IS NOT NULL
            """
        ).fetchall()
    return sort_formulas([row[0] for row in rows])


@st.cache_data
def load_default_density(formula: str) -> float | None:
    with sqlite3.connect(COMPOUNDS_DB_PATH) as connection:
        row = connection.execute(
            """
            SELECT density
            FROM Properties
            WHERE formula = ? AND density IS NOT NULL
            LIMIT 1
            """,
            (formula,),
        ).fetchone()
    return None if row is None else float(row[0])


@st.cache_data
def load_f1f2_energy_range() -> tuple[float, float]:
    min_energy: float | None = None
    max_energy: float | None = None
    for file in Path(F1F2_DIR).glob("*.ff"):
        data = pd.read_csv(file, sep=r"\s+", skiprows=7, names=["energy", "f1", "f2"], usecols=["energy"])
        file_min = float(data["energy"].min())
        file_max = float(data["energy"].max())
        min_energy = file_min if min_energy is None else min(min_energy, file_min)
        max_energy = file_max if max_energy is None else max(max_energy, file_max)

    if min_energy is None or max_energy is None:
        return DEFAULT_ENERGY_RANGE
    return min_energy, min(max_energy, MAX_ENERGY_LIMIT)


def parse_float_list(value: str) -> list[float] | None:
    values = [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]
    if not values:
        return None
    return [float(item) for item in values]


def build_coefficients_data(formulas: list[str], energies: list[float] | None, densities: dict[str, float | None]) -> dict[str, pd.DataFrame]:
    data = {}
    for formula in formulas:
        model = Coefficients(formula, energies=energies, density=densities.get(formula))
        data[formula] = model.coefficients.copy()
    return data


def add_selector(key: str, max_count: int) -> None:
    st.session_state[f"{key}_count"] = min(st.session_state[f"{key}_count"] + 1, max_count)


def remove_selector(key: str, min_count: int) -> None:
    count_key = f"{key}_count"
    current_count = st.session_state[count_key]
    if current_count > min_count:
        st.session_state[f"{key}_{current_count - 1}"] = ""
        st.session_state[count_key] = current_count - 1


def clear_selectors(key: str, min_count: int, max_count: int) -> None:
    st.session_state[f"{key}_count"] = min_count
    for index in range(max_count):
        st.session_state[f"{key}_{index}"] = ""


def select_compounds(label: str, all_formulas: list[str], key: str, min_count: int = 1, max_count: int = MAX_COMPOUND_SELECTORS) -> list[str]:
    count_key = f"{key}_count"
    if count_key not in st.session_state:
        st.session_state[count_key] = min_count

    selected = []
    options = [""] + all_formulas
    for index in range(st.session_state[count_key]):
        value = st.selectbox(
            label if index == 0 else f"{label} {index + 1}",
            options,
            key=f"{key}_{index}",
            format_func=lambda value: "Выберите соединение" if value == "" else value,
        )
        if value:
            selected.append(value)

    add_col, remove_col, clear_col = st.columns([1, 2, 7])
    add_col.button("+", key=f"{key}_add", disabled=st.session_state[count_key] >= max_count, on_click=add_selector, args=(key, max_count))
    remove_col.button(
        "Удалить последнее",
        key=f"{key}_remove",
        disabled=st.session_state[count_key] <= min_count,
        on_click=remove_selector,
        args=(key, min_count),
    )
    clear_col.button("Очистить всё", key=f"{key}_clear", on_click=clear_selectors, args=(key, min_count, max_count))

    selected = list(dict.fromkeys(selected))
    if len(selected) < min_count:
        st.info(f"Выберите минимум {min_count} соединение(я).")
        st.stop()
    return selected


def density_overrides_editor(formulas: list[str]) -> dict[str, float | None]:
    if not formulas:
        return {}

    defaults = pd.DataFrame(
        {
            "Соединение": formulas,
            "Плотность": [load_default_density(formula) for formula in formulas],
        }
    )
    edited = st.data_editor(
        defaults,
        hide_index=True,
        use_container_width=True,
        disabled=["Соединение"],
        column_config={
            "Соединение": st.column_config.TextColumn("Соединение"),
            "Плотность": st.column_config.NumberColumn("Плотность", min_value=0.0, step=0.01),
        },
    )
    return {
        row.Соединение: None if pd.isna(row.Плотность) else float(row.Плотность)
        for row in edited.itertuples(index=False)
    }


def render_density_controls(formulas: list[str]) -> dict[str, float | None]:
    with st.expander("Переопределить плотность"):
        return density_overrides_editor(formulas)


def range_state_keys(key: str) -> tuple[str, str, str, str]:
    return f"{key}_min", f"{key}_max", f"{key}_ceiling", f"{key}_slider"


def init_range_state(
    key: str,
    energy_limits: tuple[float, float],
    default_range: tuple[float, float],
    default_ceiling: float = DEFAULT_SLIDER_ENERGY_MAX,
) -> None:
    min_allowed, max_allowed = energy_limits
    default_ceiling = clamp(default_ceiling, min_allowed, max_allowed)
    min_key, max_key, ceiling_key, slider_key = range_state_keys(key)

    st.session_state.setdefault(ceiling_key, default_ceiling)
    ceiling = clamp(st.session_state[ceiling_key], min_allowed, max_allowed)
    st.session_state[ceiling_key] = ceiling
    st.session_state.setdefault(min_key, clamp(default_range[0], min_allowed, ceiling))
    st.session_state.setdefault(max_key, clamp(default_range[1], min_allowed, ceiling))

    if st.session_state[min_key] > st.session_state[max_key]:
        st.session_state[min_key] = st.session_state[max_key]
    st.session_state[slider_key] = (float(st.session_state[min_key]), float(st.session_state[max_key]))


def sync_range_from_inputs(key: str, energy_limits: tuple[float, float]) -> None:
    min_allowed, max_allowed = energy_limits
    min_key, max_key, ceiling_key, slider_key = range_state_keys(key)
    ceiling = clamp(st.session_state[ceiling_key], min_allowed, max_allowed)
    st.session_state[ceiling_key] = ceiling
    min_value = clamp(st.session_state[min_key], min_allowed, ceiling)
    max_value = clamp(st.session_state[max_key], min_allowed, ceiling)
    if min_value > max_value:
        min_value = max_value
    st.session_state[min_key] = min_value
    st.session_state[max_key] = max_value
    st.session_state[slider_key] = (float(min_value), float(max_value))


def sync_range_from_slider(key: str) -> None:
    min_key, max_key, _, slider_key = range_state_keys(key)
    min_value, max_value = st.session_state[slider_key]
    st.session_state[min_key] = float(min_value)
    st.session_state[max_key] = float(max_value)


def render_energy_range_controls(
    label: str,
    key: str,
    energy_limits: tuple[float, float],
    default_range: tuple[float, float],
    value_name: str = "энергия",
    unit: str = "эВ",
    default_ceiling: float = DEFAULT_SLIDER_ENERGY_MAX,
) -> tuple[float, float]:
    init_range_state(key, energy_limits, default_range, default_ceiling=default_ceiling)
    min_allowed, max_allowed = energy_limits
    min_key, max_key, ceiling_key, slider_key = range_state_keys(key)

    min_col, max_col, ceiling_col = st.columns(3)
    step = 10.0 if unit == "эВ" else 0.1
    min_col.number_input(
        f"Min {value_name}, {unit}",
        min_value=float(min_allowed),
        max_value=float(st.session_state[ceiling_key]),
        step=step,
        key=min_key,
        on_change=sync_range_from_inputs,
        args=(key, energy_limits),
    )
    max_col.number_input(
        f"Max {value_name}, {unit}",
        min_value=float(min_allowed),
        max_value=float(st.session_state[ceiling_key]),
        step=step,
        key=max_key,
        on_change=sync_range_from_inputs,
        args=(key, energy_limits),
    )
    ceiling_col.number_input(
        f"Максимум ползунка, {unit}",
        min_value=float(min_allowed),
        max_value=float(max_allowed),
        step=100.0 if unit == "эВ" else 0.1,
        key=ceiling_key,
        on_change=sync_range_from_inputs,
        args=(key, energy_limits),
    )
    st.slider(
        label,
        min_value=float(min_allowed),
        max_value=float(st.session_state[ceiling_key]),
        key=slider_key,
        on_change=sync_range_from_slider,
        args=(key,),
    )
    return tuple(st.session_state[slider_key])


def render_wavelength_range_controls(key: str, energy_limits: tuple[float, float], default_range: tuple[float, float] = (1.0, 15.0)) -> tuple[float, float]:
    energy_min, energy_max = energy_limits
    wavelength_limits = (HC_CONST / energy_max, HC_CONST / energy_min)
    return render_energy_range_controls(
        "Диапазон X",
        key,
        wavelength_limits,
        default_range,
        value_name="длина волны",
        unit="Å",
        default_ceiling=wavelength_limits[1],
    )


def single_energy_keys(key: str) -> tuple[str, str, str]:
    return f"{key}_value", f"{key}_ceiling", f"{key}_slider"


def init_single_energy_state(key: str, energy_limits: tuple[float, float], default_value: float = 700.0) -> None:
    min_allowed, max_allowed = energy_limits
    value_key, ceiling_key, slider_key = single_energy_keys(key)
    default_ceiling = clamp(DEFAULT_SLIDER_ENERGY_MAX, min_allowed, max_allowed)
    st.session_state.setdefault(ceiling_key, default_ceiling)
    st.session_state[ceiling_key] = clamp(st.session_state[ceiling_key], min_allowed, max_allowed)
    st.session_state.setdefault(value_key, clamp(default_value, min_allowed, st.session_state[ceiling_key]))
    st.session_state[value_key] = clamp(st.session_state[value_key], min_allowed, st.session_state[ceiling_key])
    st.session_state[slider_key] = float(st.session_state[value_key])


def sync_single_from_input(key: str, energy_limits: tuple[float, float]) -> None:
    min_allowed, max_allowed = energy_limits
    value_key, ceiling_key, slider_key = single_energy_keys(key)
    st.session_state[ceiling_key] = clamp(st.session_state[ceiling_key], min_allowed, max_allowed)
    st.session_state[value_key] = clamp(st.session_state[value_key], min_allowed, st.session_state[ceiling_key])
    st.session_state[slider_key] = float(st.session_state[value_key])


def sync_single_from_slider(key: str) -> None:
    value_key, _, slider_key = single_energy_keys(key)
    st.session_state[value_key] = float(st.session_state[slider_key])


def render_single_energy_controls(key: str, energy_limits: tuple[float, float]) -> float:
    init_single_energy_state(key, energy_limits)
    min_allowed, max_allowed = energy_limits
    value_key, ceiling_key, slider_key = single_energy_keys(key)

    energy_col, ceiling_col = st.columns(2)
    energy_col.number_input(
        "Энергия, эВ",
        min_value=float(min_allowed),
        max_value=float(st.session_state[ceiling_key]),
        step=10.0,
        key=value_key,
        on_change=sync_single_from_input,
        args=(key, energy_limits),
    )
    ceiling_col.number_input(
        "Максимум ползунка, эВ",
        min_value=float(min_allowed),
        max_value=float(max_allowed),
        step=100.0,
        key=ceiling_key,
        on_change=sync_single_from_input,
        args=(key, energy_limits),
    )
    st.slider(
        "Выбор энергии, эВ",
        min_value=float(min_allowed),
        max_value=float(st.session_state[ceiling_key]),
        key=slider_key,
        on_change=sync_single_from_slider,
        args=(key,),
    )
    return float(st.session_state[slider_key])


def render_line_chart(all_formulas: list[str], energy_limits: tuple[float, float]) -> None:
    formulas = select_compounds("Соединение", all_formulas, key="shared_compounds")
    densities = render_density_controls(formulas)

    coefficient = st.selectbox("Коэффициент", BASE_COEFFICIENTS, index=BASE_COEFFICIENTS.index("beta"))
    x_scale = st.radio("Ось X", ["Энергия, эВ", "Длина волны, Å"], horizontal=True)

    if x_scale == "Энергия, эВ":
        x_min, x_max = render_energy_range_controls("Диапазон X", "line_energy_range", energy_limits, DEFAULT_ENERGY_RANGE)
    else:
        x_min, x_max = render_wavelength_range_controls("line_wavelength_range", energy_limits)

    dots_raw = st.text_input(
        "Энергии маркеров, эВ",
        value="",
        help="Значения можно вводить через запятую, например: 100, 250, 700.",
    )
    marker_energies = parse_float_list(dots_raw)

    data = build_coefficients_data(formulas, marker_energies, densities)
    fig = plot_chart(
        deepcopy(data),
        coefficient=coefficient,
        is_x_scale_en=x_scale == "Энергия, эВ",
        x_range=[x_min, x_max],
        en_for_dots_on_chart=marker_energies if x_scale == "Энергия, эВ" else None,
        show=False,
    )
    st.plotly_chart(fig, use_container_width=True)


def render_delta_beta(all_formulas: list[str], energy_limits: tuple[float, float]) -> None:
    formulas = select_compounds("Соединение", all_formulas, key="shared_compounds")
    densities = render_density_controls(formulas)
    energy = render_single_energy_controls("delta_beta_energy", energy_limits)

    data = build_coefficients_data(formulas, [energy], densities)
    fig = delta_beta_plot(data, energy=energy, show=False)
    st.plotly_chart(fig, use_container_width=True)


def add_pair(max_count: int) -> None:
    st.session_state["mirror_pairs_count"] = min(st.session_state["mirror_pairs_count"] + 1, max_count)


def remove_pair() -> None:
    current_count = st.session_state["mirror_pairs_count"]
    if current_count > 1:
        st.session_state[f"pair_absorber_{current_count - 1}"] = ""
        st.session_state[f"pair_spacer_{current_count - 1}"] = ""
        st.session_state["mirror_pairs_count"] = current_count - 1


def clear_pairs(max_count: int) -> None:
    st.session_state["mirror_pairs_count"] = 1
    for index in range(max_count):
        st.session_state[f"pair_absorber_{index}"] = ""
        st.session_state[f"pair_spacer_{index}"] = ""


def render_material_pairs(all_formulas: list[str]) -> list[tuple[str, str]]:
    if "mirror_pairs_count" not in st.session_state:
        st.session_state["mirror_pairs_count"] = 1

    pairs = []
    options = [""] + all_formulas
    for index in range(st.session_state["mirror_pairs_count"]):
        st.markdown(f"**Пара {index + 1}**")
        absorber_col, spacer_col = st.columns(2)
        absorber = absorber_col.selectbox(
            "Absorber",
            options,
            key=f"pair_absorber_{index}",
            format_func=lambda value: "Выберите соединение" if value == "" else value,
        )
        spacer = spacer_col.selectbox(
            "Spacer",
            options,
            key=f"pair_spacer_{index}",
            format_func=lambda value: "Выберите соединение" if value == "" else value,
        )
        if absorber and spacer:
            pairs.append((absorber, spacer))

    add_col, remove_col, clear_col = st.columns([1, 2, 7])
    add_col.button("+", key="pair_add", disabled=st.session_state["mirror_pairs_count"] >= MAX_MATERIAL_PAIRS, on_click=add_pair, args=(MAX_MATERIAL_PAIRS,))
    remove_col.button("Удалить последнее", key="pair_remove", disabled=st.session_state["mirror_pairs_count"] <= 1, on_click=remove_pair)
    clear_col.button("Очистить всё", key="pair_clear", on_click=clear_pairs, args=(MAX_MATERIAL_PAIRS,))

    if not pairs:
        st.info("Выберите минимум одну пару материалов.")
        st.stop()
    return pairs


def render_pair_metric(all_formulas: list[str], energy_limits: tuple[float, float]) -> None:
    pairs = render_material_pairs(all_formulas)
    formulas = list(dict.fromkeys([formula for pair in pairs for formula in pair]))
    densities = render_density_controls(formulas)

    calculation = st.selectbox("Расчёт", ["Оптический скачок", "Максимальное отражение"])
    if calculation == "Оптический скачок":
        metric = "jump"
        add_rp = False
        add_depth = False
        polarization = "circ"
        normal_angle = 0.0
        gamma_value: float | str | None = None
    else:
        available_metrics = ["R_S_max", "R_P_max", "R_circ_max", "gamma", "RP_S", "RP_P", "L_eff"]
        metric = st.selectbox("Показатель", available_metrics)
        polarization = st.selectbox("Поляризация", ["circ", "s", "p", "all"])
        normal_angle = st.number_input("Угол к нормали", value=0.0, step=1.0)
        gamma_mode = st.selectbox("Gamma", ["0.5", "optimal", "custom"])
        if gamma_mode == "optimal":
            gamma_value = "optimal"
        elif gamma_mode == "custom":
            gamma_value = st.number_input("Значение gamma", min_value=0.0, max_value=1.0, value=0.5, step=0.01)
        else:
            gamma_value = None
        add_rp = metric in {"RP_S", "RP_P"}
        add_depth = metric == "L_eff"

    x_min, x_max = render_energy_range_controls("Диапазон энергии, эВ", "mirror_energy_range", energy_limits, DEFAULT_ENERGY_RANGE)

    data = build_coefficients_data(formulas, None, densities)
    pair_data = {}
    for absorber, spacer in pairs:
        if calculation == "Оптический скачок":
            pair_df = calculate_jump(data[absorber], data[spacer], names=[absorber, spacer])
        else:
            calc_polarization = "all" if metric in {"R_S_max", "R_P_max", "R_circ_max"} else polarization
            pair_df = calculate_max_reflection(
                data[absorber],
                data[spacer],
                polarization=calc_polarization,
                normal_angle=normal_angle,
                gamma=gamma_value,
                names=[absorber, spacer],
                add_resolving_power=add_rp,
                add_penetration_depth=add_depth,
            )
        pair_data[f"{absorber}/{spacer}"] = pair_df

    fig = plot_chart(pair_data, coefficient=metric, is_x_scale_en=True, x_range=[x_min, x_max], show=False)
    st.plotly_chart(fig, use_container_width=True)


def main() -> None:
    st.set_page_config(page_title="Графики коэффициентов", layout="wide")
    st.title("Графики коэффициентов")

    all_formulas = load_compound_formulas()
    energy_limits = load_f1f2_energy_range()
    plot_type = st.selectbox(
        "Тип графика",
        [
            "Зависимость оптических констант от E(λ)",
            "Зависимость δ от β",
            "Параметры зеркала на основе пары материалов",
        ],
    )

    try:
        if plot_type == "Зависимость оптических констант от E(λ)":
            render_line_chart(all_formulas, energy_limits)
        elif plot_type == "Зависимость δ от β":
            render_delta_beta(all_formulas, energy_limits)
        else:
            render_pair_metric(all_formulas, energy_limits)
    except Exception as exc:
        st.error(str(exc))


if __name__ == "__main__":
    main()
