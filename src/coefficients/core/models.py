"""Domain models for optical coefficients."""

import numpy as np
import pandas as pd

from ...materials.compounds import CompoundsManager
from ..data.f1f2_repository import F1F2Repository
from .constants import COEFFICIENT_LIST, HC_CONST


class Coefficients(CompoundsManager):
    """Compute optical coefficients for a compound."""

    def __init__(self, compound_formula: str, energies=None, density=None, repository: F1F2Repository | None = None):
        """Initialize coefficient model and precompute coefficient table."""
        super().__init__(compound_formula)
        self.elements = list(self.stoichiometry.keys())
        self.energies = energies
        self._repository = repository or F1F2Repository()
        if density is not None:
            self.density = density
        self.coefficients = self._calculate_coefficients()
        if self.energies:
            self.add_energies()

    def _calculate_coefficients(self):
        """Calculate coefficient dataframe from f1/f2 factors."""
        df_f1f2 = self._repository.load_weighted_factors(self.stoichiometry)
        delta = 2.7008645e-6 * self.density / self.mol_weight * np.power(df_f1f2.wavelength, 2) * df_f1f2.f1
        beta = 2.7008645e-6 * self.density / self.mol_weight * np.power(df_f1f2.wavelength, 2) * df_f1f2.f2
        n = 1 - delta
        e1 = np.power(n, 2) + np.power(beta, 2)
        e2 = 2 * n * beta
        mu = 4 * np.pi * beta / df_f1f2.wavelength
        df_all_prop = pd.DataFrame()
        for col, val in zip(COEFFICIENT_LIST[::-1], [df_f1f2.f2, df_f1f2.f1, mu, e2, e1, n, beta, delta]):
            df_all_prop.insert(loc=0, column=col, value=val)
        return df_all_prop

    def add_energies(self):
        """Add and interpolate values for requested energies."""
        df = pd.DataFrame(np.nan, index=self.energies, columns=COEFFICIENT_LIST)
        df.index.set_names("energy", inplace=True)
        self.coefficients = pd.concat([self.coefficients, df])
        self.coefficients = self.coefficients.sort_values(by="energy").interpolate(method="index")

    def get_coefficients_at_energies(self):
        """Return coefficient rows for selected energies."""
        return self.coefficients[self.coefficients.index.isin(self.energies)]

    def change_scale(self):
        """Convert index between energy and wavelength scales."""
        self.coefficients.index = HC_CONST / self.coefficients.index
        names = ["energy", "wavelength"]
        self.coefficients.index.name = [name for name in names if name != self.coefficients.index.name][0]

    def get_absorption_coefficients(self):
        """Calculate absorption coefficient from beta and wavelength."""
        ind = self.coefficients.index
        wl = ind if self.coefficients.index.name == "wavelength" else HC_CONST / ind
        return 4 * np.pi * self.coefficients["beta"] / wl

