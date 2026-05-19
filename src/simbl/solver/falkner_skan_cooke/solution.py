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

    Extends the 2D solution with crossflow profiles w(eta) and w'(eta),
    following Liu (2021), Phys. Fluids 33, 126109.

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
    w : ndarray
        Crossflow velocity ratio w/u_e
    wp : ndarray
        Crossflow velocity gradient dw/deta
    g : ndarray
        Temperature ratio T/T_e
    gp : ndarray
        Temperature gradient dg/deta
    """

    # grid
    eta: NDArray[np.float64]

    # streamwise velocity profiles
    f: NDArray[np.float64]
    fp: NDArray[np.float64]
    fpp: NDArray[np.float64]

    # crossflow velocity profiles
    w: NDArray[np.float64]
    wp: NDArray[np.float64]

    # temperature profiles
    g: NDArray[np.float64]
    gp: NDArray[np.float64]


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
        - w, wp: crossflow velocity profile and gradient
        - g, gp: temperature profile and gradient
    """
    # map solution array rows to corresponding dataclass fields
    return FalknerSkanCookeSolution(
        eta=result.eta,
        f=result.solution[0, :],
        fp=result.solution[1, :],
        fpp=result.solution[2, :],
        w=result.solution[3, :],
        wp=result.solution[4, :],
        g=result.solution[5, :],
        gp=result.solution[6, :],
    )
