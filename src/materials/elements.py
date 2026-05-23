"""Element properties service."""

import pandas as pd
from sqlalchemy import select

from src.definitions import ELEMENTS_DB_PATH
from .db import connect_sqlite_table, disconnect


class ElementsManager:
    """Load elemental properties from local database."""

    def __init__(self, element: str):
        """Initialize manager for a chemical element."""
        self.element = element
        self.properties = self._get_properties()
        self.z = self.properties["Z"].values[0]
        self.density = self.properties["Density"].values[0]
        self.at_weight = self.properties["At_weight"].values[0]
        self.boil_point = self.properties["Boil_point"].values[0]
        self.melt_point = self.properties["Melt_point"].values[0]

    def _get_properties(self):
        """Fetch one element row from `Elements.db`."""
        engine, connection, element_table = connect_sqlite_table(ELEMENTS_DB_PATH, "Element")
        query = select("*").select_from(element_table).where(element_table.c.Element == self.element)
        props = pd.DataFrame(
            [connection.execute(query).fetchall()[0][1:]],
            columns=["Z", "Element", "At_weight", "Density", "Melt_point", "Boil_point"],
        )
        disconnect(engine, connection)
        return props
