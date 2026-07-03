"""BVP collocation solver: alternative to the shooting method

Uses scipy.integrate.solve_bvp (4th-order collocation) to solve the
compressible similarity boundary layer equations.

Interface
---------
bvp_method() has the same signature as shooting_method() and returns a
ShootingResult dataclass with the same fields, so it can be used as a drop-in replacement

Requires solver_problem.bc_function to be set (done by the model builders).
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_bvp

from simbl.solver.shooting import ShootingResult, SolverProblem

if TYPE_CHECKING:
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.options import SolverOptions


# --------------------------------------------------
# bvp_method: collocation BVP solver (alternative to shooting)
# --------------------------------------------------
def bvp_method(
    solver_problem: SolverProblem,
    options: SolverOptions,
    initial_profile: ShootingResult | None = None,
    problem: SimilarityInputs | None = None,
) -> ShootingResult:
    """Solve the similarity BVP using scipy.integrate.solve_bvp (collocation)

    Parameters
    ----------
    solver_problem : SolverProblem
        Problem definition. Must have bc_function and shooting_var_indices set.
    options : SolverOptions
        Numerical settings (eta_max, n_points, tolerance used).
    initial_profile : ShootingResult, optional
        Converged solution from an adjacent parameter point (e.g. previous
        Mach) to use as the BVP initial mesh guess.  Takes priority over the
        physics-based profile.  If None, falls through to Crocco-Busemann.
    problem : SimilarityInputs, optional
        Physics specification.  When provided and initial_profile is None,
        the Crocco-Busemann profile (profiles.py) is used as the initial
        guess.  When both are None, falls back to a linear ramp.

    Returns
    -------
    ShootingResult
        Result in the same format as shooting_method. shooting_vars contains
        the wall derivative values extracted from the BVP solution.
    """

    # validate: make sure bc_function is set
    if solver_problem.bc_function is None:
        raise ValueError(
            "bvp_method requires bc_function to be set on SolverProblem. "
            "Update the model builder to supply it."
        )

    # build initial mesh
    # use exponential clustering near the wall (eta=0) where temperature and
    # velocity gradients are largest.  equidistant spacing wastes resolution
    # in the flat freestream region and starves the stiff near-wall layer.
    # solve_bvp refines adaptively from here, but it needs a good start.
    n_init = max(20, options.n_points // 20)
    eta_init = _clustered_mesh(options.eta_max, n_init)

    # build initial profile guess
    if initial_profile is not None:
        # interpolate adjacent converged solution onto the new coarse mesh
        y_init = _interpolate_profile(initial_profile.eta, initial_profile.solution, eta_init)
    elif problem is not None:
        # Crocco-Busemann physics-based profile (profiles.py)
        from simbl.solver.profiles import build_initial_profile
        y_init = build_initial_profile(problem, solver_problem, eta_init)
    else:
        # last resort: linear ramp from wall IC toward expected edge values
        y_wall = solver_problem.build_initial_condition(solver_problem.initial_guess)
        n_state = len(y_wall)
        y_init = _build_ramp_profile(y_wall, n_state, eta_init)

    # -- wrap ODE for solve_bvp's vectorized calling convention ---
    # solve_ivp calls ode_function(eta, y) with y shape (n_vars,) — one point.
    # solve_bvp calls fun(x, y) with y shape (n_vars, m) — all mesh points at once.
    # the model ODE functions were written for solve_ivp, so wrap them here.
    def ode_vectorized(x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
        n_vars, n_pts = y.shape
        dy = np.zeros_like(y)
        for i in range(n_pts):
            dy[:, i] = solver_problem.ode_function(x[i], y[:, i])
        return dy

    # -- call solve_bvp ---
    # solve_bvp struggles with tol < 1e-6 due to finite-difference Jacobian noise;
    # clamp from below but respect tighter user tolerances above 1e-6
    bvp_tol = max(options.tolerance, 1e-6)

    sol = solve_bvp(
        ode_vectorized,
        solver_problem.bc_function,
        eta_init,
        y_init,
        tol=bvp_tol,
        max_nodes=10000,
        verbose=0,
    )

    # -- evaluate on standard output grid ---
    # the BVP solver produces a continuous solution object; sample it at the
    # same grid that shooting_method would use so downstream code is identical
    eta_out = np.linspace(0.0, options.eta_max, options.n_points)
    y_out = sol.sol(eta_out)

    # -- extract shooting variable values from the wall state ---
    # shooting_var_indices tells us which state vector entries are the unknowns
    # that the shooting method would have searched for (e.g. f''(0), tau'(0))
    if solver_problem.shooting_var_indices is not None:
        shooting_vars = y_out[solver_problem.shooting_var_indices, 0]
    else:
        # fallback: guess that shooting vars are at indices 2 and 4 (FS default)
        shooting_vars = y_out[[2, 4], 0]

    # -- compute edge residual ---
    residual = solver_problem.residual_function(y_out[:, -1])

    # -- physical plausibility checks (same as shooting_method) ---
    # a genuine asymptotic solution has f'' -> 0 and tau' -> 0 at the far field;
    # a large value means the collocation converged to a non-physical branch
    converged = bool(sol.success)

    if converged and abs(y_out[2, -1]) >= options.fpp_edge_threshold:
        fpp_val = float(y_out[2, -1])
        warnings.warn(
            f"BVP: possible spurious solution — f''(eta_max) = {fpp_val:.4f} "
            f"exceeds threshold {options.fpp_edge_threshold}.",
            stacklevel=3,
        )
        converged = False

    if converged and abs(y_out[4, -1]) >= options.taup_edge_threshold:
        taup_val = float(y_out[4, -1])
        warnings.warn(
            f"BVP: possible spurious solution — tau'(eta_max) = {taup_val:.4f} "
            f"exceeds threshold {options.taup_edge_threshold}.",
            stacklevel=3,
        )
        converged = False

    return ShootingResult(
        converged=converged,
        timed_out=False,
        iterations=0,  # not applicable for collocation; field kept for compatibility
        eta=eta_out,
        solution=y_out,
        shooting_vars=shooting_vars,
        residual=residual,
    )


# --------------------------------------------------
# _interpolate_profile: map old solution onto a new eta grid
# --------------------------------------------------
def _interpolate_profile(
    eta_old: NDArray[np.float64],
    y_old: NDArray[np.float64],
    eta_new: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Interpolate a solution profile from one eta grid to another

    Parameters
    ----------
    eta_old : NDArray
        Original eta grid (1-D, sorted ascending).
    y_old : NDArray
        Solution array, shape (n_vars, len(eta_old)).
    eta_new : NDArray
        Target eta grid.

    Returns
    -------
    NDArray
        Interpolated profile, shape (n_vars, len(eta_new)).
    """
    n_vars = y_old.shape[0]
    y_new = np.zeros((n_vars, len(eta_new)))

    # interpolate each variable independently
    for i in range(n_vars):
        y_new[i, :] = np.interp(eta_new, eta_old, y_old[i, :])

    return y_new


# --------------------------------------------------
# _build_ramp_profile: simple linear initial profile guess
# --------------------------------------------------
def _build_ramp_profile(
    y_wall: NDArray[np.float64],
    n_state: int,
    eta: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Build a linear ramp profile from wall values toward expected edge values

    For each state variable, linearly interpolates from the known or guessed
    wall value to the physically expected edge value:
        f'  -> 1   (streamwise velocity reaches freestream)
        tau -> 1   (temperature reaches edge value)
        g   -> 1   (crossflow velocity reaches edge value, FSC only)
        others -> 0

    Parameters
    ----------
    y_wall : NDArray
        Initial condition vector at the wall (from build_initial_condition).
    n_state : int
        Number of state variables (5 for FS, 7 for FSC).
    eta : NDArray
        Eta grid to evaluate on.

    Returns
    -------
    NDArray
        Profile array, shape (n_state, len(eta)).
    """
    n_pts = len(eta)
    y = np.zeros((n_state, n_pts))

    # expected edge values indexed by state variable position
    # FS  (5): [f, f', f'', tau, tau']
    # FSC (7): [f, f', f'', tau, tau', g, g']
    edge_values = np.zeros(n_state)
    if n_state >= 2:
        edge_values[1] = 1.0   # f' -> 1
    if n_state >= 4:
        edge_values[3] = 1.0   # tau -> 1
    if n_state >= 6:
        edge_values[5] = 1.0   # g -> 1 (crossflow; 1 for swept, 0 for aligned, use 1 as placeholder)

    # linear ramp for each variable from wall value to expected edge value
    for i in range(n_state):
        y[i, :] = np.linspace(y_wall[i], edge_values[i], n_pts)

    return y


# --------------------------------------------------
# _clustered_mesh: exponentially stretched eta grid, dense near wall
# --------------------------------------------------
def _clustered_mesh(eta_max: float, n: int, alpha: float = 4.0) -> NDArray[np.float64]:
    """Build an exponentially clustered eta grid concentrated near eta=0

    The temperature and velocity gradients are largest near the wall (eta=0)
    and decay toward the freestream.  Clustering points there gives the BVP
    solver much better initial resolution in the stiff near-wall layer.

    The mapping is:
        eta(t) = eta_max * (exp(alpha * t) - 1) / (exp(alpha) - 1),  t in [0, 1]

    At alpha=4: ~50% of points fall in the first 15% of the domain.
    At alpha=0 (limit): recovers uniform spacing.

    Parameters
    ----------
    eta_max : float
        Outer edge of the domain.
    n : int
        Number of mesh points.
    alpha : float
        Clustering strength.  Larger = more points near the wall.

    Returns
    -------
    NDArray
        Monotonically increasing eta grid of length n, starting at 0.
    """
    # uniform parameter t in [0, 1]
    t = np.linspace(0.0, 1.0, n)

    # exponential map: clusters points near t=0 (i.e. eta=0)
    eta = eta_max * (np.exp(alpha * t) - 1.0) / (np.exp(alpha) - 1.0)

    return eta
