"""Numerical settings for the similarity solver

These settings are purely numerical and do not affect the physical problem being solved
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import warnings
from dataclasses import dataclass

from simbl.solver.equations import EQUATIONS_ALIASES, VALID_EQUATIONS


# --------------------------------------------------
# SolverOptions dataclass (immutable)
#
# store solver configuration (numerical settings) independent of the physical problem
# --------------------------------------------------
@dataclass(frozen=True)
class SolverOptions:
    """Numerical settings for the similarity solver

    Parameters
    ----------
    eta_max : float, optional
        Maximum similarity coordinate. Default 15.0.
    n_points : int, optional
        Number of grid points. Default 500.
    tolerance : float, optional
        Convergence tolerance for shooting method. Default 1e-8.
    max_iterations : int, optional
        Maximum shooting iterations. Default 100.
    ode_method : str, optional
        SciPy ODE solver method. Default "LSODA".
    epsilon : float, optional
        Perturbation step size for finite-difference Jacobian. Default 0.005.
    relaxation_factor : float, optional
        Factor to reduce epsilon by when the Jacobian is singular. Default 0.5.
    equations : str, optional
        Governing equations: 'falkner_skan' (2D) or 'falkner_skan_cooke' (3D). Default 'falkner_skan'.

    Examples
    --------
    Default options (suitable for most cases):

    >>> options = SolverOptions()

    Higher resolution and tighter tolerance:

    >>> options = SolverOptions(
    ...     n_points=500,
    ...     eta_max=20.0,
    ...     tolerance=1e-10,
    ... )
    """

    # --------------------------------------------------
    # define parameters with type annotations and default values
    # --------------------------------------------------

    # maximum similarity coordinate
    eta_max: float = 15.0
    # number of points to discretize the similarity coordinate eta
    n_points: int = 500
    # convergence tolerance for shooting method
    tolerance: float = 1e-8
    # maximum number of iterations for shooting method
    max_iterations: int = 100
    # ODE solver method (passed to SciPy's solve_ivp)
    ode_method: str = "LSODA"
    # perturbation step size for finite-difference Jacobian
    epsilon: float = 0.005
    # factor to reduce epsilon by when the Jacobian is singular
    relaxation_factor: float = 0.5
    # governing equations to solve
    equations: str = "falkner_skan"

    # --------------------------------------------------
    # post-init validation (runs automatically after dataclass __init__)
    #
    # sanity checks to ensure numerical settings are reasonable
    #
    # hard error checks first: prevent code from running
    # soft warnings second: alert user to potentially problematic settings but allow code to run
    # --------------------------------------------------
    def __post_init__(self) -> None:
        """Validate solver options"""

        # --------------------------------------------------
        # hard error checks: raise exceptions for invalid settings
        # --------------------------------------------------
        if self.eta_max <= 0:
            raise ValueError(f"eta_max must be positive: {self.eta_max}")
        if self.n_points <= 10:
            raise ValueError(f"n_points must be > 10: {self.n_points}")
        if self.tolerance <= 0:
            raise ValueError(f"tolerance must be positive: {self.tolerance}")
        if self.max_iterations < 1:
            raise ValueError(f"max_iterations must be >= 1: {self.max_iterations}")
        if self.epsilon <= 0:
            raise ValueError(f"epsilon must be positive: {self.epsilon}")
        if not 0 < self.relaxation_factor < 1:
            raise ValueError(f"relaxation_factor must be between 0 and 1: {self.relaxation_factor}")
        # normalize short aliases (fs, fsc) to their canonical form
        object.__setattr__(self, "equations", EQUATIONS_ALIASES.get(self.equations, self.equations))
        if self.equations not in VALID_EQUATIONS:
            raise ValueError(
                f"equations must be one of {set(VALID_EQUATIONS)}: got {self.equations!r}"
            )

        # --------------------------------------------------
        # soft warnings: alert user to potentially problematic settings but allow code to run
        # --------------------------------------------------
        if self.eta_max < 5.0:
            warnings.warn(
                f"eta_max = {self.eta_max} may be too small to reach freestream conditions",
                stacklevel=2,
            )
        if self.max_iterations < 5:
            warnings.warn(
                f"max_iterations = {self.max_iterations} is very low, solver may not converge",
                stacklevel=2,
            )
