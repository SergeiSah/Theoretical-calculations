"""Optical formulas based on absorber/spacer coefficients."""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve

from ..core.constants import HC_CONST


def merge_abs_sps(*dfs: pd.DataFrame) -> pd.DataFrame:
    """Merge absorber and spacer tables on energy with interpolation."""
    assert len(dfs) == 2, "There are must be two DataFrames: 1) for absorber, 2) for spacer"
    a = dfs[0].loc[:, ["delta", "beta"]]
    s = dfs[1].loc[:, ["delta", "beta"]]
    return a.merge(s, how="outer", on="energy", suffixes=["_a", "_s"]).sort_index().interpolate(method="index").dropna(axis=0)


def calculate_jump(absorber: pd.DataFrame, spacer: pd.DataFrame, names=None) -> pd.DataFrame:
    """Compute optical jump between absorber and spacer."""
    merged_df = merge_abs_sps(absorber, spacer)
    merged_df["jump"] = abs(merged_df.delta_a - merged_df.delta_s) / abs(merged_df.beta_a)
    if isinstance(names, list) and len(names) == 2:
        merged_df.columns = [f"delta_{names[0]}", f"beta_{names[0]}", f"delta_{names[1]}", f"beta_{names[1]}", "jump"]
    return merged_df


def polarization_func(pol: str, angle: float):
    """Return polarization factor for s or p polarization."""
    if pol.lower() == "s":
        return 1.0
    if pol.lower() == "p":
        ang_rad = angle / 180 * np.pi
        return np.cos(2 * ang_rad)
    raise TypeError("Wrong polarization type")


def eq_for_optimal_gamma(_gamma, _a):
    """Equation whose root gives optimal gamma."""
    return np.tan(np.pi * _gamma) - np.pi * _gamma - np.pi * _a


def calculate_resolving_power(abs_sps: pd.DataFrame, pol: str, **kwargs):
    """Compute resolving power for selected polarization."""
    f = kwargs["f"]
    ang = kwargs["normal_angle"]
    ang_rad = ang / 180 * np.pi
    pol_value = polarization_func(pol, ang)
    y = kwargs["y"]
    c = 2 * np.sqrt(2) / 3 * np.pi
    im_mu = abs_sps.gamma * abs_sps.beta_a + (1 - abs_sps.gamma) * abs_sps.beta_s
    abs_sps["RP"] = pd.NA
    for i, en in enumerate(abs_sps.index.values):
        xi = abs(f.iloc[i] * pol_value) / y.iloc[i]
        if xi >= 10:
            abs_sps.iloc[i, 5] = c / np.sin(np.pi * abs_sps.loc[en, "gamma"]) * np.power(np.cos(ang_rad), 2) / abs(pol_value) / (
                abs_sps.loc[en, "delta_a"] - abs_sps.loc[en, "delta_s"]
            )
        else:
            abs_sps.iloc[i, 5] = np.power(np.cos(ang_rad), 2) / im_mu.iloc[i] / np.power(
                (1 - np.power(pol_value / y.iloc[i], 2)) * (1 + np.power(f.iloc[i] * pol_value / y.iloc[i], 2)),
                0.75,
            )
    return abs_sps["RP"]


def calculate_penetration_depth(abs_sps: pd.DataFrame, normal_angle: float, **kwargs):
    """Compute effective penetration depth."""
    wavelengths = pd.DataFrame(HC_CONST / abs_sps.index.values, columns=["wv"]).set_index(abs_sps.index)
    im_mu = abs_sps.gamma * abs_sps.beta_a + (1 - abs_sps.gamma) * abs_sps.beta_s
    return wavelengths.wv * np.cos(normal_angle) / np.pi / im_mu / np.sqrt(kwargs["one_parenthesis"] * kwargs["two_parenthesis"])


def calculate_max_reflection(
    absorber: pd.DataFrame,
    spacer: pd.DataFrame,
    polarization="circ",
    normal_angle=0,
    gamma=None,
    names=None,
    add_resolving_power=False,
    add_penetration_depth=False,
) -> pd.DataFrame:
    """Compute maximal reflectivity and optional derived metrics."""
    df = merge_abs_sps(absorber, spacer)
    f = (df.delta_a - df.delta_s) / (df.beta_a - df.beta_s)
    g = df.beta_s / (df.beta_a - df.beta_s)

    if gamma is None:
        df["gamma"] = 0.5
    elif str(gamma).lower() in ["o", "opt", "optimal"]:
        df["gamma"] = pd.NA
        for i, _ in enumerate(df.index.values):
            df.iloc[i, 4] = fsolve(eq_for_optimal_gamma, np.array([0.499999]), (g.iloc[i],))
        df.gamma = df.gamma.astype(np.float64)
    elif isinstance(gamma, float):
        df["gamma"] = gamma
    else:
        raise TypeError

    df.loc[df.gamma < 0, "gamma"] = df.loc[df.gamma < 0, "gamma"].apply(lambda x: 1 + x)
    y = np.pi * (df.gamma + g) / np.sin(np.pi * df.gamma)
    s_pol = polarization_func("s", normal_angle)
    p_pol = polarization_func("p", normal_angle)
    u_s = np.sqrt((1 - np.power(s_pol / y, 2)) / (1 + np.power(f * s_pol / y, 2)))
    u_p = np.sqrt((1 - np.power(p_pol / y, 2)) / (1 + np.power(f * p_pol / y, 2)))
    r_s_max = (1 - u_s) / (1 + u_s)
    r_p_max = (1 - u_p) / (1 + u_p)

    if add_resolving_power:
        df["RP_S"] = calculate_resolving_power(df, "s", f=f, normal_angle=normal_angle, y=y)
        df["RP_P"] = calculate_resolving_power(df, "p", f=f, normal_angle=normal_angle, y=y)

    if add_penetration_depth:
        df["L_eff"] = calculate_penetration_depth(
            df,
            normal_angle,
            one_parenthesis=(1 - np.power(s_pol / y, 2)),
            two_parenthesis=(1 + np.power(f * s_pol / y, 2)),
        )

    if polarization.lower() == "s":
        df["R_S_max"] = r_s_max
    elif polarization.lower() == "p":
        df["R_P_max"] = r_p_max
    elif polarization.lower() in ["c", "circ", "circular"]:
        df["R_circ_max"] = (r_s_max + r_p_max) / 2
    else:
        df["R_S_max"] = r_s_max
        df["R_P_max"] = r_p_max
        df["R_circ_max"] = (r_s_max + r_p_max) / 2

    if names is not None:
        df.columns = [k.replace("_a", f"_{names[0]}").replace("_s", f"_{names[1]}") for k in df.columns]
    return df
