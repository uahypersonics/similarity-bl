"""Numerical settings for the similarity solver

These settings are purely numerical and do not affect the physical problem being solved
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from simbl.solver.equations import EQUATIONS_ALIASES, VALID_EQUATIONS

if TYPE_CHECKING:
    from simbl.solver.inputs import SimilarityInputs


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
    max_solve_time : float, optional
        Maximum wall-clock time [seconds] for the entire shooting method call. Default 20.0.

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
    # maximum wall-clock time [seconds] for the entire shooting method call.
    # the Newton loop checks elapsed time after each iteration and aborts early
    # if this limit is exceeded.  default 20 s is well within the 30 s process
    # timeout used by the lookup table generator.
    max_solve_time: float = 20.0
    # primary solver to use:
    #   "shooting" -- Newton shooting method (default, fast)
    #   "bvp"      -- scipy.integrate.solve_bvp collocation (robust, slower)
    solver_method: str = "shooting"
    # when solver_method="shooting" and bvp_fallback=True, retry with the BVP
    # collocation solver if shooting fails to converge.  the BVP receives the
    # last converged profile from the caller (initial_profile) as its starting
    # mesh — a physically correct profile from an adjacent parameter point.
    bvp_fallback: bool = False
    # post-convergence physical plausibility threshold for f''(eta_max).
    # a genuine asymptotic solution has f'' -> 0 at the far field; a large
    # value indicates Newton converged to a spurious oscillatory branch.
    # set to None to disable the check.
    fpp_edge_threshold: float = 0.01
    # post-convergence physical plausibility threshold for tau'(eta_max).
    # a genuine asymptotic solution has tau' -> 0 at the far field; a large
    # value indicates the temperature profile has not relaxed to the edge
    # condition (spurious branch or eta_max too small).
    # calibrated for eta_max >= 5: at eta_max=5 physical solutions reach
    # |tau'(eta_max)| ~ 0.13-0.29 for cold walls / high Mach, so the
    # original threshold of 0.05 was incorrectly rejecting physical solutions.
    taup_edge_threshold: float = 0.5
    # divisor for the BDF fallback max_step inside shooting_method.
    # max_step = eta_max / bdf_max_step_divisor.
    # larger values → smaller steps → more stable for stiff near-wall gradients.
    # default 50 preserves the original behaviour (max_step = eta_max/50 = 0.12 at eta_max=6).
    bdf_max_step_divisor: float = 50.0
    # when True, record the shooting variable values and residual norm at every
    # Newton iteration.  stored in ShootingResult.history as a list of dicts
    # with keys "iteration", "shooting_vars", "residual_norm".
    # off by default to avoid overhead in large parameter sweeps.
    store_history: bool = False

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
        # validate solver_method
        valid_methods = {"shooting", "bvp"}
        if self.solver_method not in valid_methods:
            raise ValueError(
                f"solver_method must be one of {valid_methods}: got {self.solver_method!r}"
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


# --------------------------------------------------
# default_options: physics-informed SolverOptions
#
# derive numerical settings from the physical problem description
# so that an inexperienced user does not have to tune them manually
# --------------------------------------------------
def default_options(inputs: SimilarityInputs) -> SolverOptions:
    """Return SolverOptions with physics-informed defaults for the given inputs.

    Derives numerical settings (eta_max, equations) directly from the physical
    problem description so that an inexperienced user does not need to tune them
    manually.  All other options take the standard dataclass defaults.

    Parameters
    ----------
    inputs : SimilarityInputs
        Physical problem description (Mach, beta, wall BC, equations, etc.).

    Returns
    -------
    SolverOptions
        Solver options calibrated for the given inputs.

    Examples
    --------
    >>> from simbl import SimilarityInputs, default_options, solve_similarity
    >>> inputs = SimilarityInputs(mach_edge=0.5, temp_edge=300.0, beta=0.8)
    >>> options = default_options(inputs)
    >>> solution, result = solve_similarity(inputs, options)
    """

    # --------------------------------------------------
    # heuristic to set eta_max from beta
    # --------------------------------------------------

    # For favorable pressure gradient (beta > 0) the boundary layer is thinner
    # and the solution reaches freestream conditions at smaller eta.
    eta_max = max(6.0, 15.0 - 5.0 * inputs.beta)

    # --------------------------------------------------
    # set equations to match the physical problem
    # --------------------------------------------------

    # non-zero sweep angle means Falkner-Skan-Cooke (3-D); otherwise plain FS
    equations = "falkner_skan_cooke" if inputs.sweep_angle != 0.0 else "falkner_skan"

    # --------------------------------------------------
    # build and return options
    # --------------------------------------------------
    return SolverOptions(equations=equations, eta_max=eta_max)
