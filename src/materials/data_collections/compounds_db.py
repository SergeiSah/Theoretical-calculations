"""API for writing compound records into local `Compounds.db`."""

import pandas as pd
from chemparse import parse_formula
from sqlalchemy import MetaData, Table, create_engine, insert, select

from src.definitions import COMPOUNDS_DB_PATH


class CompoundsDB:
    """Insert normalized compound records into database tables."""

    def __init__(self):
        """Open database and reflect required tables."""
        metadata = MetaData()
        self.engine = create_engine(f"sqlite+pysqlite:///{COMPOUNDS_DB_PATH}")
        self.connection = self.engine.connect()
        self.compounds = Table("Compound", metadata, autoload_with=self.engine)
        self.sources = Table("Source", metadata, autoload_with=self.engine)
        self.densities = Table("Density", metadata, autoload_with=self.engine)

    def add_record(self, properties: dict):
        """Insert one compound density record if it does not exist."""
        if self.is_id_in_db(properties["id"]):
            print("The record is already in the database")
            return

        if not self.is_source_in_db(properties["source"]):
            self.connection.execute(insert(self.sources).values(name=properties["source"]))

        query = select(self.sources.c.id).select_from(self.sources).where(self.sources.c.name == properties["source"])
        source_id = self.connection.execute(query).fetchone()[0]

        if not self.is_compound_in_db(properties["formula"]):
            self.connection.execute(insert(self.compounds).values(name=properties["name"], formula=properties["formula"]))

        query = select(self.compounds.c.id, self.compounds.c.formula).select_from(self.compounds)
        comps = pd.DataFrame(self.connection.execute(query).fetchall(), columns=["id", "formula"])
        compound_id = comps[comps["formula"].apply(parse_formula) == parse_formula(properties["formula"])]["id"].values[0]

        self.connection.execute(
            insert(self.densities).values(
                compound_id=int(compound_id),
                density=properties["density"],
                mol_weight=properties["mol_weight"],
                source_id=source_id,
                id_in_source=properties["id"],
            )
        )

    def is_compound_in_db(self, compound) -> bool:
        """Check whether chemically equivalent formula exists in DB."""
        query = select(self.compounds.c.formula).select_from(self.compounds)
        chem_comp = pd.DataFrame(self.connection.execute(query).fetchall(), columns=["formula"])
        return not chem_comp[chem_comp["formula"].apply(parse_formula) == parse_formula(compound)].empty

    def is_source_in_db(self, source) -> bool:
        """Check whether source exists in DB."""
        query = select(self.sources.c.name).select_from(self.sources).where(self.sources.c.name == source)
        return len(self.connection.execute(query).fetchall()) > 0

    def is_id_in_db(self, id_in_source) -> bool:
        """Check whether source-specific record id exists in DB."""
        query = select(self.densities.c.id_in_source).select_from(self.densities).where(self.densities.c.id_in_source == id_in_source)
        return len(self.connection.execute(query).fetchall()) > 0
