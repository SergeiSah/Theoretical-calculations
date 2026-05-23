"""Utility for duplicate-energy checks in `.ff` files."""

import os
from pathlib import Path

from src.definitions import F1F2_DIR

import pandas as pd


def main():
    """Print duplicate energy count for each `.ff` file."""
    base_dir = Path(F1F2_DIR)
    for file in os.listdir(base_dir):
        if file.endswith(".ff"):
            df = pd.read_csv(base_dir / file, sep=r"\s+", skiprows=7, names=["energy", "f1", "f2"])
            sum_dup = df.duplicated(subset=["energy"]).sum()
            print(f"{file}: {sum_dup}")


if __name__ == "__main__":
    main()
