"""Plotting utilities for optical coefficient datasets."""

import itertools

import pandas as pd
import plotly
import plotly.graph_objects as go

from ..core.constants import DEF_X_RANGE_EN, DEF_X_RANGE_WL, HC_CONST

COLORS = plotly.colors.qualitative.Plotly


def get_y_range(x_range: list, data, col_name: str, margin_coef=0.05):
    """Calculate y-axis range for given x-window and metric."""
    if isinstance(data, dict):
        data = list(data.values())
    min_y = min(data[0].loc[x_range[0] : x_range[1], col_name])
    max_y = max(data[0].loc[x_range[0] : x_range[1], col_name])
    for df in data:
        series = df.loc[x_range[0] : x_range[1], col_name]
        min_y = min(series) if min(series) < min_y else min_y
        max_y = max(series) if max(series) > max_y else max_y
    margin = abs(max_y - min_y) * margin_coef
    return [min_y - margin, max_y + margin]


def plot_chart(data: dict, coefficient="beta", is_x_scale_en=True, x_range=None, en_for_dots_on_chart=None, show=True):
    """Render line chart for selected coefficient across compounds."""
    if x_range is None and is_x_scale_en:
        x_range = DEF_X_RANGE_EN
    elif x_range is None and not is_x_scale_en:
        x_range = DEF_X_RANGE_WL

    if not is_x_scale_en:
        for name in data:
            data[name].index = HC_CONST / data[name].index
            data[name].sort_index(inplace=True)

    y_range = get_y_range(x_range, data, col_name=coefficient)
    col_iter = itertools.cycle(COLORS)
    traces = [go.Scatter(x=data[name].index, y=data[name][coefficient], mode="lines", name=name, line=dict(color=next(col_iter))) for name in data.keys()]

    col_iter = itertools.cycle(COLORS)
    if en_for_dots_on_chart is not None:
        for name in data.keys():
            dots = data[name][data[name].index.isin(en_for_dots_on_chart)]
            traces.append(
                go.Scatter(
                    x=dots.index,
                    y=dots[coefficient],
                    mode="markers+text",
                    name=name,
                    showlegend=False,
                    text=dots[coefficient].round(3),
                    textposition="middle right",
                    textfont=dict(color="black", size=14),
                    line=dict(color=next(col_iter)),
                )
            )

    fig = go.Figure(data=traces)
    fig.update_xaxes(range=x_range, title="Энергия, эВ" if is_x_scale_en else "Длина волны, Å")
    fig.update_yaxes(range=y_range, fixedrange=False)
    y_label = {
        "delta": "Refractive index decrement",
        "beta": "Absorption index",
        "n": "Real part of the refractive index",
        "e1": "Real part of the permittivity",
        "e2": "Imaginary part of the permittivity",
        "jump": "Optical jump",
        "R_max": "R_max",
        "R_S_max": "R_S_max",
        "R_P_max": "R_P_max",
        "R_circ_max": "R_circ_max",
        "RP_S": "Resolving power (s)",
        "RP_P": "Resolving power (p)",
        "gamma": "d(a)/d",
        "L_eff": "Penetration depth (A)",
        "mu": "Absorption coefficient (1/A)",
    }
    fig.update_yaxes(title=y_label[coefficient])
    if show:
        fig.show()
    return fig


def delta_beta_plot(data: dict, energy: float, show=True):
    """Plot delta-beta scatter for all compounds at one energy."""
    x, y, compounds = [], [], []
    for elem in data.keys():
        row = data[elem][data[elem].index == energy]
        x.append(row.delta.values[0])
        y.append(row.beta.values[0])
        compounds.append(elem)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, text=compounds, textposition="bottom center", mode="markers+text"))
    fig.update_traces(marker_size=15, marker_line_width=2, marker_color="green")
    fig.update_layout(width=1000, height=800, font_size=18, font_family="Arial")
    fig.update_xaxes(title="Refractive index decrement", title_font_size=30, exponentformat="e", showexponent="last")
    fig.update_yaxes(title="Absorption index", title_font_size=30, exponentformat="e", showexponent="last")
    fig.add_annotation(
        yref="paper",
        y=0.98,
        xref="paper",
        x=0.02,
        text=f"<b>Энергия:</b> {energy} эВ",
        showarrow=False,
        bgcolor="#fdca40",
        borderpad=6,
        font=dict(size=24),
    )
    if show:
        fig.show()
    return fig
