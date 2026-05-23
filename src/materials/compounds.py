"""Compound properties service."""

import pandas as pd
from chemparse import parse_formula
from mp_api.client import MPRester
from sqlalchemy import select

from src.definitions import COMPOUNDS_DB_PATH, get_materials_project_credentials
from .db import connect_sqlite_table, disconnect
from .elements import ElementsManager

pd.set_option("display.max_columns", None)


class CompoundsManager:
    """Read and derive compound properties from local DB or Materials Project."""

    def __init__(self, compound_formula: str, fill_with_mp=True, row_num_in_properties=0, avg_val=False):
        """Initialize manager for one chemical formula."""
        self.compound_formula = compound_formula
        self.__row_in_prop = row_num_in_properties
        self.__avg_val = avg_val
        self.__fill_with_mp = fill_with_mp
        self._density_override = None

    @property
    def stoichiometry(self):
        """Return parsed stoichiometry dictionary."""
        return parse_formula(self.compound_formula)

    @property
    def properties(self):
        """Return candidate property rows for the formula."""
        res = self._get_properties_table(self.stoichiometry)
        if res.empty and self.__fill_with_mp:
            return self._download_from_mp()
        return res

    @property
    def density(self):
        """Return selected or overridden density value."""
        if self._density_override is not None:
            return self._density_override
        if self.properties.empty:
            return None
        if self.__avg_val:
            return self.properties["density"].mean()
        return self.properties.loc[self.__row_in_prop, "density"]

    @density.setter
    def density(self, dens: float):
        """Override density value for calculations."""
        self._density_override = dens

    @property
    def mol_weight(self):
        """Return molar mass for selected row."""
        if self.properties.empty:
            return None
        return self.properties.loc[self.__row_in_prop, "mol_weight"]

    @property
    def chem_name(self):
        """Return chemical name for selected row."""
        if self.properties.empty:
            return None
        return self.properties.loc[self.__row_in_prop, "chem_name"]

    @property
    def row_num_in_df(self):
        """Return currently selected row index."""
        return self.__row_in_prop

    @row_num_in_df.setter
    def row_num_in_df(self, num: int):
        """Set selected row index if it exists."""
        if num > self.properties.shape[0] - 1:
            print("Entered number is higher than maximum possible row number in the dataframe")
        else:
            self.__row_in_prop = num

    @property
    def avg_density(self):
        """Return flag that controls density averaging."""
        return self.__avg_val

    @avg_density.setter
    def avg_density(self, avg_val: bool):
        """Enable or disable averaging of densities."""
        self.__avg_val = avg_val

    def _get_properties_table(self, stoichiometry: dict) -> pd.DataFrame:
        """Read all property rows matching stoichiometry from local DB view."""
        engine, connection, properties_table = connect_sqlite_table(COMPOUNDS_DB_PATH, "Properties")
        query = select("*").select_from(properties_table)
        properties_all = pd.DataFrame(
            connection.execute(query).fetchall(),
            columns=["formula", "chem_name", "density", "mol_weight", "source_name", "id"],
        )
        disconnect(engine, connection)
        return properties_all[properties_all["formula"].apply(parse_formula) == stoichiometry].reset_index(drop=True)

    def _download_from_mp(self):
        """Download fallback properties from Materials Project API."""
        api_key = get_materials_project_credentials().api_key
        if not api_key:
            raise RuntimeError("Materials Project API key is not set. Add MAPI_KEY to the project root .env file.")

        with MPRester(api_key) as mpr:
            results = mpr.summary.search(formula=[self.compound_formula])[0]

        return pd.DataFrame(
            {
                "formula": [self.compound_formula],
                "chem_name": [""],
                "density": [results.density],
                "mol_weight": [self._calculate_mol_weight()],
                "source_name": ["Materials Project"],
            }
        )

    def _calculate_mol_weight(self):
        """Calculate molar mass from element stoichiometry."""
        sum_mol_weight = 0
        for elem in self.stoichiometry.keys():
            sum_mol_weight += ElementsManager(elem).at_weight * self.stoichiometry[elem]
        return sum_mol_weight
