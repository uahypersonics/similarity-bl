"""Newton-Raphson shooting method

Solves boundary value problems using the shooting method with
finite-difference Jacobian. Model-agnostic -- operates only
through the SolverProblem abstraction.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import time
import warnings
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from simbl.solver.options import SolverOptions


# --------------------------------------------------
# SolverProblem dataclass
# --------------------------------------------------
@dataclass
class SolverProblem:
    """Internal interface: model-specific code -> generic shooter

    Built by each model's build_solver_problem() function.
    The generic shooting method operates exclusively through
    this abstraction, with no direct imports from model packages.

    Attributes
    ----------
    ode_function : Callable
        ODE right-hand side with signature f(eta, y) -> dy.
    build_initial_condition : Callable[[NDArray], NDArray]
        Maps shooting variables to initial condition vector y(0).
    residual_function : Callable[[NDArray], NDArray]
        Maps edge state y(eta_max) to residual vector.
    n_shooting : int
        Number of shooting variables.
    initial_guess : NDArray[np.float64]
        Starting values for the shooting variables.
    """

    ode_function: Callable[[float, NDArray[np.float64]], NDArray[np.float64]]
    build_initial_condition: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    residual_function: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    n_shooting: int
    initial_guess: NDArray[np.float64]
    # boundary condition function for solve_bvp: bc(ya, yb) -> residual vector
    # encodes wall BCs in ya and edge BCs in yb; set by the model builder
    bc_function: Callable[[NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]] | None = None
    # indices into the state vector y that correspond to the shooting variables
    # (the unknowns the shooting method optimises at the wall)
    # used by bvp_method to extract equivalent wall values from the BVP solution
    shooting_var_indices: list[int] | None = None


# --------------------------------------------------
# ShootingResult dataclass
# --------------------------------------------------
@dataclass
class ShootingResult:
    """Result from shooting method"""

    converged: bool
    """Indicates whether the solution converged within the specified tolerance and iteration limit"""
    timed_out: bool
    """True when the shooting method aborted because max_solve_time was exceeded"""
    iterations: int
    """Number of iterations required for convergence (or max iterations if not converged)"""
    eta: NDArray[np.float64]
    """Similarity coordinate array corresponding to the solution"""
    solution: NDArray[np.float64]
    """Solution array at each eta point, shape (n_variables, n_eta)"""
    shooting_vars: NDArray[np.float64]
    """Converged shooting variable values at the wall (e.g., f''(0), g(0) or g'(0))"""
    residual: NDArray[np.float64]
    """Final residual vector at the edge (should be close to zero if converged)"""
    history: list[dict] | None = None
    """Per-iteration convergence history when options.store_history=True.
    Each entry: {"iteration": int, "shooting_vars": NDArray, "residual_norm": float}.
    None when store_history=False (default)."""

    # special method (dunder): controls how the object is printed / repr'd
    def __repr__(self) -> str:

        # format shooting variables and residual as compact inline lists
        svars = ", ".join(f"{v:.6g}" for v in self.shooting_vars)
        resid = ", ".join(f"{v:.3e}" for v in self.residual)
        # build a readable multi-line block
        lines = [
            "ShootingResult",
            f"  converged   : {self.converged}",
            f"  timed_out   : {self.timed_out}",
            f"  iterations  : {self.iterations}",
            f"  eta range   : [{self.eta[0]:.4g}, {self.eta[-1]:.4g}]  ({len(self.eta)} points)",
            f"  solution    : shape {self.solution.shape}",
            f"  shoot vars  : [{svars}]",
            f"  residual    : [{resid}]",
        ]
        return "\n".join(lines)


# --------------------------------------------------
# Newton-Raphson shooting method
# --------------------------------------------------
def shooting_method(
    solver_problem: SolverProblem,
    options: SolverOptions,
) -> ShootingResult:
    """Solve the BVP using shooting method with Newton-Raphson iteration

    Parameters
    ----------
    solver_problem : SolverProblem
        Model-specific problem definition containing:
        - ODE function (bound to problem parameters and visc_model)
        - initial condition builder (maps shooting variables to y0)
        - residual function (maps edge values to residual vector)
        - number of shooting variables
        - initial guess for shooting variables (from lookup table or user override)
    options : SolverOptions
        Numerical settings containing:
        - eta_max: maximum value of similarity coordinate
        - n_points: number of points in eta discretization
        - tolerance: convergence tolerance for Newton-Raphson iteration
        - max_iterations: maximum number of Newton-Raphson iterations
        - ode_method: integration method for solve_ivp (e.g., 'RK45', 'BDF')
        - epsilon: perturbation step size for finite-difference Jacobian
        - relaxation_factor: factor to reduce epsilon by when the Jacobian is singular

    Returns
    -------
    ShootingResult
        Shooting method result, including convergence status, solution array, and final shooting variable values
    """

    # discretize eta (similarity coordinate) based on specified max and number of points
    eta = np.linspace(0, options.eta_max, options.n_points)

    # perturbation step size for centered finite-difference Jacobian (default 0.005, see options.py)
    epsilon = np.full(solver_problem.n_shooting, options.epsilon)

    # initialize shooting variables
    s = solver_problem.initial_guess.copy()

    # wall-clock deadline for the whole shooting call.
    # once the budget is spent, the next derivative
    # evaluation raises and aborts that integration cleanly.
    t_shoot_start = time.monotonic()
    deadline = t_shoot_start + options.max_solve_time

    class _SolveTimeout(Exception):
        """Raised from the ODE RHS when the wall-clock budget is exhausted."""

    def _timed_ode(eta_val: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        if time.monotonic() > deadline:
            raise _SolveTimeout
        return solver_problem.ode_function(eta_val, y)

    def _integrate(s_guess, *, edge_only=False):
        """Integrate ODE with trial shooting variables.

        When edge_only=True, only stores the solution at eta_max
        (skips interpolating the full grid). Used for Jacobian
        columns where we only need the boundary state.
        """
        y0 = solver_problem.build_initial_condition(s_guess)
        points = [options.eta_max] if edge_only else eta

        # guard: non-finite y0 (from ill-conditioned Newton steps) causes
        # solve_ivp to raise ValueError before integration starts; catch it
        # and return a sentinel object with non-finite edge values so the
        # Newton iteration can detect the bad step and reduce epsilon
        class _FailedSol:
            success = False
            y = np.full((len(y0), len(points)), np.nan)

        try:
            sol = solve_ivp(
                _timed_ode,
                (0, options.eta_max),
                y0,
                method=options.ode_method,
                t_eval=points,
                rtol=options.tolerance,
                atol=options.tolerance,
            )
        except (ValueError, OverflowError, _SolveTimeout):
            return _FailedSol()

        if not sol.success:
            # fallback: BDF is an implicit method, more stable for stiff ODEs
            try:
                sol = solve_ivp(
                    _timed_ode,
                    (0, options.eta_max),
                    y0,
                    method="BDF",
                    t_eval=points,
                    rtol=options.tolerance,
                    atol=options.tolerance,
                    max_step=options.eta_max / options.bdf_max_step_divisor,
                )
            except (ValueError, OverflowError, _SolveTimeout):
                return _FailedSol()
        return sol

    # --------------------------------------------------
    # Newton-Raphson iteration loop
    #
    # shooting variables s = [f''(0), g_var, ...] are the unknowns
    # each iteration:
    #   1. integrate the ODE from eta=0 to eta=eta_max using current s
    #   2. evaluate the residual F(s) = y(eta_max) - target boundary conditions
    #   3. build the Jacobian J = dF/ds via centered finite differences
    #   4. solve J * delta_s = F for the Newton correction
    #   5. update s <- s - delta_s
    # repeat until ||F|| < tolerance
    # --------------------------------------------------

    # initialize convergence flag
    converged = False
    # flag set when the wall-clock limit aborts the loop before convergence
    timed_out = False
    # initialize iteration counter
    iteration = 0
    # initialize residual vector
    F = np.zeros(solver_problem.n_shooting)
    # convergence history: populated only when options.store_history=True
    history: list[dict] | None = [] if options.store_history else None

    for iteration in range(1, options.max_iterations + 1):  # noqa: B007

        # check wall-clock limit before doing any work this iteration
        # (t_shoot_start / deadline are set above, before the integrator helper)
        if time.monotonic() > deadline:
            timed_out = True
            break  # converged stays False; loop exits cleanly

        # integrate ODE with current shooting variables
        sol = _integrate(s)

        # solution sol.y is arranged as (n_variables, n_eta); we want the edge values at eta_max for the residual -> slice the last column
        y_edge = sol.y[:, -1]

        # guard: if the ODE integration diverged (both LSODA and BDF failed),
        # y_edge may contain inf/nan; clamp to a large finite value so that the
        # residual function receives valid input and the Newton iteration can
        # detect the large residual and take a corrective step rather than
        # propagating non-finite values through the Jacobian build
        y_edge = np.where(np.isfinite(y_edge), y_edge, np.sign(y_edge + 1) * 1.0e30)

        # compute residual vector F(s) = y(eta_max) - target boundary conditions (set in builder.py of the respective model)
        F = solver_problem.residual_function(y_edge)

        # record convergence history if requested
        if history is not None:
            history.append({
                "iteration": iteration,
                "shooting_vars": s.copy(),
                "residual_norm": float(np.linalg.norm(F)),
            })

        # check convergence (tolerance defaults to 1e-8 if it is not user defined, see options.py)
        if np.all(np.abs(F) < options.tolerance):
            converged = True
            break

        # build Jacobian via centered finite differences:
        #   J[:, j] = (F(s + eps_j * e_j) - F(s - eps_j * e_j)) / (2 * eps_j)
        # each column j requires two perturbed integrations (forward and backward)

        # initialize Jacobian matrix
        J = np.zeros((solver_problem.n_shooting, solver_problem.n_shooting))

        # loop over each shooting variable to compute finite difference columns
        for j in range(solver_problem.n_shooting):

            # true copy of shooting variables for forward and backward perturbations (avoid in-place modification issues)
            s_fwd = s.copy()
            s_bwd = s.copy()

            # perturb the j-th shooting variable forward and backward by epsilon[j]
            s_fwd[j] += epsilon[j]
            s_bwd[j] -= epsilon[j]

            # integrate with forward and backward perturbed shooting variables, only storing edge values for residual evaluation
            # edge_only=True speeds up integration by skipping interpolation of the full solution array on the whole eta grid
            sol_fwd = _integrate(s_fwd, edge_only=True)
            sol_bwd = _integrate(s_bwd, edge_only=True)

            # extract edge values from the forward and backward solutions for residual evaluation
            # guard against malformed sol.y (e.g., non-array result from integration failure)
            # use a large-value sentinel so the residual function gets valid (bad) values
            n_state = solver_problem.build_initial_condition(s).shape[0]
            y_raw_fwd = sol_fwd.y
            y_raw_bwd = sol_bwd.y
            if not isinstance(y_raw_fwd, np.ndarray) or y_raw_fwd.ndim != 2:
                y_raw_fwd = np.full((n_state, 1), np.nan)
            if not isinstance(y_raw_bwd, np.ndarray) or y_raw_bwd.ndim != 2:
                y_raw_bwd = np.full((n_state, 1), np.nan)
            y_edge_fwd = y_raw_fwd[:, -1]
            y_edge_bwd = y_raw_bwd[:, -1]

            # clamp non-finite edge values so Jacobian columns stay finite
            y_edge_fwd = np.where(np.isfinite(y_edge_fwd), y_edge_fwd, np.sign(y_edge_fwd + 1) * 1.0e30)
            y_edge_bwd = np.where(np.isfinite(y_edge_bwd), y_edge_bwd, np.sign(y_edge_bwd + 1) * 1.0e30)

            # compute residuals for forward and backward perturbations
            F_fwd = solver_problem.residual_function(y_edge_fwd)
            F_bwd = solver_problem.residual_function(y_edge_bwd)

            # compute the j-th column of the Jacobian using centered finite difference
            J[:, j] = (F_fwd - F_bwd) / (2.0 * epsilon[j])

        # newton-raphson update
        try:
            # solve for the update to shooting variables: J * delta_s = F => delta_s = J^{-1} * F
            delta_s = np.linalg.solve(J, F)

            # proposed next shooting variables
            s_next = s - delta_s

            # guard: if the Newton step produced non-finite shooting variables
            # (can happen when the Jacobian is ill-conditioned and F contains
            # inf/nan from ODE overflow during perturbation integrations),
            # reject the step and reduce epsilon instead of propagating inf
            if not np.all(np.isfinite(s_next)):
                epsilon *= options.relaxation_factor
                continue

            # update shooting variables for the next iteration
            s = s_next

        except np.linalg.LinAlgError:
            # singular Jacobian => cannot solve for update => reduce epsilon and try again
            # a smaller epsilon can resolve near-singular columns caused by
            # perturbations that overshoot into a flat region of the residual
            epsilon *= options.relaxation_factor

            continue

    #--------------------------------------------------
    # final integration with converged shooting variables
    #--------------------------------------------------
    sol = _integrate(s)

    # --------------------------------------------------
    # post-convergence physical plausibility check
    #
    # a genuine asymptotic solution has f''(eta_max) -> 0 as the velocity
    # profile flattens to the freestream value.  a large f''(eta_max) means
    # Newton converged to a spurious oscillatory branch that satisfies the
    # far-field BC at one point by coincidence, not by asymptotic decay.
    # --------------------------------------------------
    if converged and abs(sol.y[2, -1]) >= options.fpp_edge_threshold:
        fpp_edge_val = float(sol.y[2, -1])
        warnings.warn(
            f"Possible spurious solution: f''(eta_max) = {fpp_edge_val:.4f} "
            f"exceeds threshold {options.fpp_edge_threshold}. "
            "Newton may have converged to a non-physical oscillatory branch. "
            "Consider using a better initial guess or beta-continuation.",
            stacklevel=3,
        )
        converged = False

    if converged and abs(sol.y[4, -1]) >= options.taup_edge_threshold:
        taup_edge_val = float(sol.y[4, -1])
        warnings.warn(
            f"Possible spurious solution: tau'(eta_max) = {taup_edge_val:.4f} "
            f"exceeds threshold {options.taup_edge_threshold}. "
            "Temperature profile has not relaxed to the edge condition; "
            "Newton may have converged to a non-physical branch or eta_max is too small.",
            stacklevel=3,
        )
        converged = False

    # -------------------------------------------------
    # package results into ShootingResult dataclass for convenient access by caller
    # -------------------------------------------------
    return ShootingResult(
        converged=converged,
        timed_out=timed_out,
        iterations=iteration,
        eta=eta,
        solution=sol.y,
        shooting_vars=s.copy(),
        residual=F,
        history=history,
    )
