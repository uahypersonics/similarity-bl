"""Falkner-Skan-Cooke problem builder

Constructs a SolverProblem for the 3D compressible Falkner-Skan-Cooke model
(7-equation ODE, 3 shooting variables).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from simbl.solver.falkner_skan_cooke.initial_guess import build_y0, get_initial_guess
from simbl.solver.falkner_skan_cooke.ode import bl_ode
from simbl.solver.shooting import SolverProblem

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from flow_state.transport import TransportModel

    from simbl.solver.inputs import SimilarityInputs


# --------------------------------------------------
# problem builder
# called from main.py to assemble everything the generic shooting method
# needs: ODE function, initial condition builder, residual, and initial guess
# --------------------------------------------------
def build_solver_problem(
    problem: SimilarityInputs,
    visc_model: TransportModel,
    *,
    initial_fpp: float | None = None,
    initial_wp: float | None = None,
    initial_gvar: float | None = None,
) -> SolverProblem:
    """Build a SolverProblem for the Falkner-Skan-Cooke model

    Parameters
    ----------
    problem : SimilarityInputs
        Physics specification (including sweep_angle).
    visc_model : TransportModel
        Viscosity model instance.
    initial_fpp : float, optional
        Override initial guess for f''(0).
    initial_wp : float, optional
        Override initial guess for w'(0).
    initial_gvar : float, optional
        Override initial guess for g'(0) (isothermal) or g(0) (adiabatic).

    Returns
    -------
    SolverProblem
        Ready-to-solve problem for the generic shooter.
    """
    # --------------------------------------------------
    # resolve initial guess
    # user can partially override initial guess or rely entirely on the lookup table
    # --------------------------------------------------
    if initial_fpp is not None and initial_wp is not None and initial_gvar is not None:
        # if user provided all three shooting variables, use them directly (precedence over lookup table)
        guess = np.array([initial_fpp, initial_wp, initial_gvar])
    else:
        # some lookup logic is needed

        # build key_values dict for lookup table
        key_values = {"mach": problem.mach_edge, "beta": problem.beta, "sweep_angle": problem.sweep_angle}
        if problem.wall_bc == "isothermal":
            key_values["g_wall"] = problem.g_wall

        # get initial guess from lookup table based on problem parameters (wired in initial_guess module)
        guess = get_initial_guess(wall_bc=problem.wall_bc, key_values=key_values).copy()

        # apply any user overrides if provided on top of lookup prediction
        if initial_fpp is not None:
            guess[0] = initial_fpp
        if initial_wp is not None:
            guess[1] = initial_wp
        if initial_gvar is not None:
            guess[2] = initial_gvar

    # --------------------------------------------------
    # ODE setup: bind problem, visc_model, and tan^2(sweep) into a two argument function for solve_ivp
    # --------------------------------------------------
    sweep_rad = np.radians(problem.sweep_angle)
    tan2_sweep = np.tan(sweep_rad) ** 2

    def ode_func(eta: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return bl_ode(eta, y, problem, visc_model, tan2_sweep)

    # --------------------------------------------------
    # initial condition builder (3 shooting variables -> 7-element y0)
    # --------------------------------------------------
    def build_initial_condition(s: NDArray[np.float64]) -> NDArray[np.float64]:
        return build_y0(wall_bc=problem.wall_bc, g_wall=problem.g_wall, shooting_vars=s)

    # --------------------------------------------------
    # residual function: f'(inf) = 1, w(inf) = w_e_normalized, g(inf) = 1
    # for aligned flow (sweep_angle = 0), w_e = 0 so w(inf) = 0
    # for swept flow (sweep_angle > 0), w is normalized by edge crossflow so w(inf) = 1
    # --------------------------------------------------
    w_edge_target = 0.0 if problem.sweep_angle == 0.0 else 1.0

    def residual(y_edge: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.array([
            y_edge[1] - 1.0,         # f'(inf) - 1
            y_edge[3] - w_edge_target,  # w(inf) - w_edge_target
            y_edge[5] - 1.0,         # g(inf) - 1
        ])

    # --------------------------------------------------
    # pack into SolverProblem dataclass (convenience wrapper) for the generic shooting method (defined in shooting.py)
    # --------------------------------------------------
    return SolverProblem(
        ode_function=ode_func,
        build_initial_condition=build_initial_condition,
        residual_function=residual,
        n_shooting=3,
        initial_guess=guess,
    )
