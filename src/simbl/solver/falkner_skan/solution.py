"""Falkner-Skan solution dataclass and builder"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.shooting import ShootingResult


# --------------------------------------------------
# solution dataclass
# --------------------------------------------------
@dataclass(frozen=True)
class FalknerSkanSolution:
    """Container for Falkner-Skan similarity solution data

    Attributes
    ----------
    eta : ndarray
        Similarity coordinate
    f : ndarray
        Stream function
    fp : ndarray
        Velocity ratio u/u_e
    fpp : ndarray
        Velocity gradient d^2f/deta^2
    tau : ndarray
        Temperature ratio T/T_e
    taup : ndarray
        Temperature gradient dg/deta
    """

    # grid
    eta: NDArray[np.float64]

    # streamwise velocity profiles
    f: NDArray[np.float64]
    fp: NDArray[np.float64]
    fpp: NDArray[np.float64]

    # temperature profiles
    tau: NDArray[np.float64]
    taup: NDArray[np.float64]


# --------------------------------------------------
# solution builder
# --------------------------------------------------
def build_solution(result: ShootingResult) -> FalknerSkanSolution:
    """Extract profiles from shooting result into solution dataclass

    Parameters
    ----------
    result : ShootingResult
        Converged shooting method result with 5-row solution array

    Returns
    -------
    FalknerSkanSolution
        Named solution profiles
        Contains:
        - eta: similarity coordinate
        - f, fp, fpp: stream function and velocity profiles
        - tau, taup: temperature profile and gradient
    """

    # map solution array rows to corresponding dataclass fields
    # # state vector: [f, fp, fpp, tau, taup]
    return FalknerSkanSolution(
        eta=result.eta,
        f=result.solution[0, :],
        fp=result.solution[1, :],
        fpp=result.solution[2, :],
        tau=result.solution[3, :],
        taup=result.solution[4, :],
    )
