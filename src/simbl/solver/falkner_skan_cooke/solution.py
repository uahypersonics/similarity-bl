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

    # -- boundary layer thickness methods ---

    def delta99(self) -> float:
        """Boundary layer thickness in eta space.

        Returns the similarity coordinate eta at which u/u_e = 0.99.
        To convert to physical wall-normal distance, multiply by the
        Illingworth-Stewartson length scale appropriate for the local station.

        Returns:
            eta value where fp(eta) = 0.99.
        """
        # interpolate to find the eta crossing u/u_e = 0.99
        return float(np.interp(0.99, self.fp, self.eta))

    def delta_star(self) -> float:
        """Compressible displacement thickness in eta space.

        Computes integral_0^eta_max (tau - fp) d_eta, proportional to the
        physical displacement thickness delta*. To convert to physical units,
        multiply by the Illingworth-Stewartson length scale.

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
        multiply by the Illingworth-Stewartson length scale.

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
