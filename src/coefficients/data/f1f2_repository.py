"""Repository for atomic scattering factors f1/f2."""

from pathlib import Path

import pandas as pd

from src.definitions import F1F2_DIR

from ..core.constants import HC_CONST


class F1F2Repository:
    """Load and aggregate f1/f2 tables for compound stoichiometry."""

    def __init__(self, base_dir: Path | None = None):
        """Initialize repository with directory that contains `.ff` files."""
        self.base_dir = base_dir or Path(F1F2_DIR)

    def load_weighted_factors(self, stoichiometry: dict[str, float]) -> pd.DataFrame:
        """Load factors for each element and aggregate with stoichiometric weights."""
        factors = []
        for elem, st in stoichiometry.items():
            factor = pd.read_csv(
                self.base_dir / f"{elem}.ff",
                header=7,
                sep=r"\s+",
                names=["energy", "f1", "f2"],
            ).set_index("energy")
            factor["f1"] = factor["f1"] * st
            factor["f2"] = factor["f2"] * st
            factors.append(factor)

        if len(factors) == 1:
            factors[0].insert(0, "wavelength", HC_CONST / factors[0].index)
            return factors[0][factors[0].index > 50]

        merged_f = factors[0]
        for factor in factors[1:]:
            merged_f = merged_f.merge(factor, how="outer", left_index=True, right_index=True)

        factors = [factor.reindex(merged_f.index).interpolate(method="index") for factor in factors]
        sum_factors = factors[0]
        sum_factors["f1"] = sum(factor["f1"] for factor in factors)
        sum_factors["f2"] = sum(factor["f2"] for factor in factors)
        sum_factors.insert(0, "wavelength", HC_CONST / sum_factors.index)
        return sum_factors[sum_factors.index > 50]
