"""Newton-Raphson shooting method

Solves boundary value problems using the shooting method with
finite-difference Jacobian. Model-agnostic -- operates only
through the SolverProblem abstraction.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

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


# --------------------------------------------------
# ShootingResult dataclass
# --------------------------------------------------
@dataclass
class ShootingResult:
    """Result from shooting method"""

    converged: bool
    """Indicates whether the solution converged within the specified tolerance and iteration limit"""
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
                solver_problem.ode_function,
                (0, options.eta_max),
                y0,
                method=options.ode_method,
                t_eval=points,
                rtol=options.tolerance,
                atol=options.tolerance,
            )
        except (ValueError, OverflowError):
            return _FailedSol()

        if not sol.success:
            # fallback: BDF is an implicit method, more stable for stiff ODEs
            try:
                sol = solve_ivp(
                    solver_problem.ode_function,
                    (0, options.eta_max),
                    y0,
                    method="BDF",
                    t_eval=points,
                    rtol=options.tolerance,
                    atol=options.tolerance,
                )
            except (ValueError, OverflowError):
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
    # initialize iteration counter
    iteration = 0
    # initialize residual vector
    F = np.zeros(solver_problem.n_shooting)

    for iteration in range(1, options.max_iterations + 1):  # noqa: B007

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

    # -------------------------------------------------
    # package results into ShootingResult dataclass for convenient access by caller
    # -------------------------------------------------
    return ShootingResult(
        converged=converged,
        iterations=iteration,
        eta=eta,
        solution=sol.y,
        shooting_vars=s.copy(),
        residual=F,
    )
