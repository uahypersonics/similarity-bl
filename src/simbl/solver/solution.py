"""Canonical solved similarity-profile data shared across SIMBL workflows."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray


# --------------------------------------------------
# similarity solution
# --------------------------------------------------
@dataclass(frozen=True, slots=True)
class SimilaritySolution:
    """Validated similarity profiles produced by a SIMBL solver."""

    eta: NDArray[np.float64]
    f: NDArray[np.float64]
    fp: NDArray[np.float64]
    fpp: NDArray[np.float64]
    tau: NDArray[np.float64]
    taup: NDArray[np.float64]
    metadata: dict[str, Any]
    g: NDArray[np.float64] | None = None
    gp: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        """Validate profile dimensions, coordinates, and values."""

        # validate all profile arrays against the eta coordinate
        profiles = {
            "eta": self.eta,
            "f": self.f,
            "fp": self.fp,
            "fpp": self.fpp,
            "tau": self.tau,
            "taup": self.taup,
        }
        for name, values in profiles.items():
            if values.ndim != 1:
                raise ValueError(f"profiles.{name} must be one-dimensional")
            if values.shape != self.eta.shape:
                raise ValueError(f"profiles.{name} must have the same shape as eta")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"profiles.{name} must contain only finite values")

        # validate optional crossflow profiles as one complete pair
        if (self.g is None) != (self.gp is None):
            raise ValueError("profiles.g and profiles.gp must be provided together")
        if self.g is not None and self.gp is not None:
            crossflow_profiles = {"g": self.g, "gp": self.gp}
            for name, values in crossflow_profiles.items():
                if values.ndim != 1 or values.shape != self.eta.shape:
                    raise ValueError(f"profiles.{name} must have the same shape as eta")
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"profiles.{name} must contain only finite values")

        # require enough ordered points for cubic interpolation
        if self.eta.size < 4:
            raise ValueError("similarity profile must contain at least four points")
        if not np.all(np.diff(self.eta) > 0.0):
            raise ValueError("profiles.eta must be strictly increasing")
