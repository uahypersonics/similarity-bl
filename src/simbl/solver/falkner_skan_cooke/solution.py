"""Falkner-Skan-Cooke solution dataclass and builder"""

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
class FalknerSkanCookeSolution:
    """Container for Falkner-Skan-Cooke solution data

    Extends the 2D Falkner-Skan solution with crossflow profiles g(eta) and g'(eta),

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
        Temperature ratio T/T_e (Liu's tau, avoiding 'g' notation for clarity)
    taup : ndarray
        Temperature gradient dtau/deta
    g : ndarray
        Crossflow velocity ratio w/w_e (Liu's g, denoted g_cf to avoid confusion with temperature)
    gp : ndarray
        Crossflow velocity gradient dg_cf/deta
    """

    # grid
    eta: NDArray[np.float64]

    # streamwise velocity profiles
    f: NDArray[np.float64]
    fp: NDArray[np.float64]
    fpp: NDArray[np.float64]

    # crossflow velocity profiles (g_cf = w/w_e, avoiding notation confusion)
    g: NDArray[np.float64]
    gp: NDArray[np.float64]

    # temperature profiles (tau = T/T_e, using Liu's notation)
    tau: NDArray[np.float64]
    taup: NDArray[np.float64]


# --------------------------------------------------
# solution builder
# --------------------------------------------------
def build_solution(result: ShootingResult) -> FalknerSkanCookeSolution:
    """Extract profiles from shooting result into solution dataclass

    Parameters
    ----------
    result : ShootingResult
        Converged shooting method result with 7-row solution array

    Returns
    -------
    FalknerSkanCookeSolution
        Named solution profiles
        Contains:
        - eta: similarity coordinate
        - f, fp, fpp: stream function and velocity profiles
        - tau, taup: temperature profile and gradient
        - g, gp: crossflow velocity profile and gradient
    """

    # map solution array rows to corresponding dataclass fields
    # state vector: [f, fp, fpp, tau, taup, g, gp]
    return FalknerSkanCookeSolution(
        eta=result.eta,
        f=result.solution[0, :],
        fp=result.solution[1, :],
        fpp=result.solution[2, :],
        tau=result.solution[3, :],
        taup=result.solution[4, :],
        g=result.solution[5, :],
        gp=result.solution[6, :],
    )
