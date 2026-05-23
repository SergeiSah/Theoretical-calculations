"""Photoionization cross-sections and asymmetry data service."""

import itertools

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objects as go
from sqlalchemy import MetaData, Table, create_engine, select

from src.definitions import PICS_DB_PATH


class PICS:
    """Provide access to PICS and asymmetry factors."""

    __PICS_COLUMNS = ["sub_level", "j", "energy", "sigma"]
    __ASYM_COLUMNS = ["sub_level", "j", "energy", "beta", "delta", "gamma"]
    __DB_SHORT_NAMES = {
        **dict.fromkeys(["T", "Trzh", "Trzhaskovskaya"], "Trzhaskovskaya"),
        **dict.fromkeys(["Y", "Yeh"], "Yeh"),
        **dict.fromkeys(["S", "Scof", "Scofield"], "Scofield"),
    }
    __MAX_X_AXES = 10000
    __COLORS = plotly.colors.qualitative.Light24

    def __init__(self, element: str, db_author="T", year: int = None):
        """Initialize PICS data for one element and source."""
        self.element = element
        self.db_author = self.__DB_SHORT_NAMES[db_author]
        self.year = year
        pics, asym = self._get_data_from_pics_db()
        self.pics_df = pd.DataFrame(pics, columns=self.__PICS_COLUMNS).set_index("energy")
        self.asym_df = pd.DataFrame(asym, columns=self.__ASYM_COLUMNS).set_index("energy")
        self.sub_levels = pd.unique(self.pics_df.iloc[:, [0, 1]].apply(lambda x: "".join(x), axis=1))
        self.pics_dict = {
            sl: self.pics_df[self.pics_df.iloc[:, [0, 1]].apply(lambda x: "".join(x), axis=1) == sl].iloc[:, 2:]
            for sl in self.sub_levels
        }
        self.asym_dict = {
            sl: self.asym_df[self.asym_df.iloc[:, [0, 1]].apply(lambda x: "".join(x), axis=1) == sl].iloc[:, 2:]
            for sl in self.sub_levels
        }

    def _get_data_from_pics_db(self):
        """Query PICS and asymmetry tables for current element/source/year."""
        metadata = MetaData()
        engine = create_engine(f"sqlite+pysqlite:///{PICS_DB_PATH}")
        connection = engine.connect()
        sigmas = Table("pics_view", metadata, autoload_with=engine)
        asymmetries = Table("asymmetries_view", metadata, autoload_with=engine)

        if self.year is None:
            self.year = min(connection.execute(select(sigmas.c.pub_year).where(sigmas.c.author == self.db_author)).fetchall())[0]

        query_pics = (
            select(sigmas.c.sub_level, sigmas.c.j, sigmas.c.energy, sigmas.c.sigma)
            .where(sigmas.c.element == self.element)
            .where(sigmas.c.author == self.db_author)
            .where(sigmas.c.pub_year == self.year)
        )
        query_asym = (
            select(
                asymmetries.c.sub_level,
                asymmetries.c.j,
                asymmetries.c.energy,
                asymmetries.c.beta,
                asymmetries.c.gamma,
                asymmetries.c.delta,
            )
            .where(asymmetries.c.element == self.element)
            .where(asymmetries.c.author == self.db_author)
            .where(asymmetries.c.pub_year == self.year)
        )
        pics = connection.execute(query_pics).fetchall()
        asym = connection.execute(query_asym).fetchall()
        connection.close()
        engine.dispose()
        return pics, asym

    def add_value_for_energy(self, energy: float) -> None:
        """Interpolate one energy point for all subshell tables."""
        for data_map in [self.pics_dict, self.asym_dict]:
            for subshell in self.sub_levels:
                self.is_in_range(subshell, energy)
                if energy not in data_map[subshell].index:
                    try:
                        data_map[subshell] = (
                            pd.concat([data_map[subshell], pd.DataFrame(index=[energy], columns=data_map[subshell].columns)])
                            .sort_index()
                            .interpolate(method="pchip", limit_direction="backward")
                        )
                    except TypeError:
                        pass

    def is_in_range(self, subshell: str, energy: float) -> bool:
        """Validate that requested energy is inside stored range."""
        en_range = [min(self.pics_dict[subshell].index), max(self.pics_dict[subshell].index)]
        assert en_range[0] < energy < en_range[1], f'Energy "{energy}" eV is out of available range: {en_range[0]} eV - {en_range[1]}'
        return True

    def plot_subshell_pics(self, line: str) -> None:
        """Plot cross section curve for one subshell."""
        line_to_plot = self.add_scatters_for_line(self.pics_dict[line])
        scatters = self.pics_dict[line]
        fig = go.Figure()
        color = next(itertools.cycle(self.__COLORS))
        fig.add_trace(go.Scatter(x=line_to_plot.index, y=line_to_plot.sigma, line={"color": color}, name=line))
        fig.add_trace(go.Scatter(x=scatters.index, y=scatters.sigma, mode="markers", showlegend=False, line={"color": color}, name=line))
        fig.update_yaxes(title="Photoionization cross sections (kb)", tickformat=".2e", type="log")
        fig.update_xaxes(title="Energy of photons (eV)", tickformat="g")
        fig.show()

    def plot_all_pics(self) -> None:
        """Plot cross section curves for all subshells."""
        fig = go.Figure()
        min_energy = 1e10
        max_energy = 0
        col_palette_iter = itertools.cycle(self.__COLORS)
        for subshell in self.pics_dict:
            color = next(col_palette_iter)
            data = self.add_scatters_for_line(self.pics_dict[subshell])
            fig.add_trace(go.Scatter(x=data.index, y=data.sigma, name=subshell, line={"color": color}))
            fig.add_trace(
                go.Scatter(
                    x=self.pics_dict[subshell].index,
                    y=self.pics_dict[subshell].sigma,
                    mode="markers",
                    name=subshell,
                    line={"color": color},
                    showlegend=False,
                )
            )
            min_energy = min(data.index) if min_energy > min(data.index) else min_energy
            max_energy = max(data.index) if max(data.index) < self.__MAX_X_AXES else self.__MAX_X_AXES
        fig.update_yaxes(type="log", title="Photoionization cross sections (kb)", tickformat=".2e")
        fig.update_xaxes(title="Energy of photons (eV)", tickformat="g", range=[min_energy, max_energy])
        fig.show()

    @staticmethod
    def add_scatters_for_line(data, step=100):
        """Generate dense interpolated curve for plotting."""
        energies = np.arange(min(data.index), 10000, step)
        to_plot = data
        for energy in energies:
            to_plot = pd.concat([to_plot, pd.DataFrame(index=[energy], columns=to_plot.columns)])
        return to_plot.sort_index().interpolate(method="pchip", limit_direction="backward").drop_duplicates()

    def get_pics(self, line: str, energy: float) -> pd.DataFrame:
        """Return cross section value at requested energy."""
        self.add_value_for_energy(energy)
        return self.pics_dict[line][self.pics_dict[line].index == energy]

    def get_asymmetries(self, line: str, energy: float) -> pd.DataFrame:
        """Return asymmetry coefficients at requested energy."""
        self.add_value_for_energy(energy)
        return self.asym_dict[line][self.asym_dict[line].index == energy]

    def get_asymmetry_factor(self, line: str, energy: float, angle: float) -> pd.DataFrame:
        """Compute angular asymmetry factor for one subshell."""
        radians = angle / 180 * np.pi
        beta, gamma, delta = self.get_asymmetries(line, energy).values[0]
        if pd.isna(beta):
            return pd.DataFrame(index=[energy], data={"asymmetry factor": pd.NA})
        gamma = 0.0 if pd.isna(gamma) else gamma
        delta = 0.0 if pd.isna(delta) else delta
        as_factor = 1 + 0.5 * beta * (1.5 * np.power(np.sin(radians), 2) - 1) + (gamma / 2 * np.power(np.sin(radians), 2) + delta) * np.cos(radians)
        return pd.DataFrame(index=[energy], data={"asymmetry factor": as_factor})

    def get_rsf(self, line: str, energy: float):
        """Placeholder for relative sensitivity factor calculation."""
        pass
