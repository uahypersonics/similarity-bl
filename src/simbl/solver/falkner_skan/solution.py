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

    # -- boundary layer thickness methods ---

    def delta99(self) -> float:
        """Boundary layer thickness in eta space.

        Returns the similarity coordinate eta at which u/u_e = 0.99.
        To convert to physical wall-normal distance, multiply by the
        Levy-Lees length scale L = sqrt(2 * nu_e * x / u_e).

        Returns:
            eta value where fp(eta) = 0.99.
        """
        # interpolate to find the eta crossing u/u_e = 0.99
        return float(np.interp(0.99, self.fp, self.eta))

    def delta_star(self) -> float:
        """Compressible displacement thickness in eta space.

        Computes integral_0^eta_max (tau - fp) d_eta, proportional to the
        physical displacement thickness delta*. To convert to physical units,
        multiply by L = sqrt(2 * nu_e * x / u_e).

        For the incompressible limit (tau = 1 everywhere), this reduces to
        the standard form integral (1 - fp) d_eta.

        Returns:
            Dimensionless displacement thickness in eta space.
        """
        # integrate (tau - fp): compressible correction via density rho/rho_e = 1/tau
        # absorbs into dy/d_eta = tau * L, leaving (tau - fp) as the integrand
        return float(np.trapezoid(self.tau - self.fp, self.eta))

    def theta(self) -> float:
        """Compressible momentum thickness in eta space.

        Computes integral_0^eta_max fp * (1 - fp) d_eta, proportional to the
        physical momentum thickness theta. To convert to physical units,
        multiply by L = sqrt(2 * nu_e * x / u_e).

        Note: the density weighting (rho/rho_e = 1/tau) cancels exactly with
        the dy/d_eta = tau * L factor, so the eta-space integrand is fp * (1 - fp)
        regardless of the level of compressibility.

        Returns:
            Dimensionless momentum thickness in eta space.
        """
        # integrate fp * (1 - fp) to get momentum thickness
        return float(np.trapezoid(self.fp * (1.0 - self.fp), self.eta))


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
