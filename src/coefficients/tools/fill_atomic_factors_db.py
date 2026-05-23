"""Utility for importing `.ff` data into `Atomic_factors.db`."""

import os
from pathlib import Path

import pandas as pd
from sqlalchemy import MetaData, Table, create_engine

from src.definitions import ATOMIC_FACTORS_DB_PATH, F1F2_DIR


def main():
    """Append all `.ff` records into `atomic_factors` SQLite table."""
    base_dir = Path(F1F2_DIR)
    db = create_engine(f"sqlite+pysqlite:///{ATOMIC_FACTORS_DB_PATH}")
    metadata = MetaData()
    conn = db.connect()
    Table("Atomic_factors", metadata, autoload_with=db)

    for file in os.listdir(base_dir):
        if file.endswith(".ff"):
            element = file.split(".")[0]
            df = pd.read_csv(base_dir / file, sep=r"\s+", skiprows=7, names=["energy", "f1", "f2"])
            df.insert(0, "element", element)
            df.to_sql("atomic_factors", conn, if_exists="append", index=False)

    conn.close()
    db.dispose()


if __name__ == "__main__":
    main()
