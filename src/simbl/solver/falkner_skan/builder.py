"""Falkner-Skan problem builder

Constructs a SolverProblem for the 2D compressible Falkner-Skan model
(5-equation ODE, 2 shooting variables).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from simbl.solver.falkner_skan.initial_guess import build_y0, get_initial_guess
from simbl.solver.falkner_skan.ode import bl_ode
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
    initial_gvar: float | None = None,
) -> SolverProblem:
    """Build a SolverProblem for the Falkner-Skan model

    Parameters
    ----------
    problem : SimilarityInputs
        Physics specification.
    visc_model : TransportModel
        Viscosity model instance.
    initial_fpp : float, optional
        Override initial guess for f''(0).
    initial_gvar : float, optional
        Override initial guess for g variable -- g(0) for adiabatic, g'(0) for isothermal.

    Returns
    -------
    SolverProblem
        Ready-to-solve problem for the generic shooter.
    """
    # --------------------------------------------------
    # resolve initial guess
    # user can partially override initial guess or rely entirely on the lookup table
    # --------------------------------------------------
    if initial_fpp is not None and initial_gvar is not None:
        # if user provided both shooting variables, use them directly (precedence over lookup table)
        guess = np.array([initial_fpp, initial_gvar])
    else:
        # some lookup logic is needed

        # build key_values dict for lookup table
        key_values = {"mach": problem.mach_edge, "beta": problem.beta}
        if problem.wall_bc == "adiabatic":
            key_values["temp_edge"] = problem.temp_edge
        elif problem.wall_bc == "isothermal":
            key_values["g_wall"] = problem.g_wall

        # get initial guess from lookup table based on problem parameters (wired in initial_guess module)
        guess = get_initial_guess(wall_bc=problem.wall_bc, key_values=key_values).copy()

        # apply any user overrides if provided on top of lookup prediction
        if initial_fpp is not None:
            guess[0] = initial_fpp
        if initial_gvar is not None:
            guess[1] = initial_gvar

        # for adiabatic walls, override the lookup-table g_wall estimate with the
        # analytical recovery temperature formula.  The shipped table was generated at
        # T_edge=300 K; at other T_edge values the IDW prediction is inaccurate because
        # Sutherland viscosity is non-linear in absolute temperature.
        # The recovery formula is a reliable starting point for all Mach / T_edge combos.
        if problem.wall_bc == "adiabatic":
            # recovery factor r ~ Pr^(1/3) for turbulent; use same value for laminar FS
            r = problem.prandtl ** (1.0 / 3.0)
            g_wall_est = 1.0 + r * (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2
            guess[1] = g_wall_est

    # --------------------------------------------------
    # ODE setup: bind problem and visc_model into a two argument function for solve_ivp
    # --------------------------------------------------
    def ode_func(eta: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return bl_ode(eta, y, problem, visc_model)

    # --------------------------------------------------
    # initial condition builder (2 shooting variables -> 5-element y0)
    # --------------------------------------------------
    def build_initial_condition(s: NDArray[np.float64]) -> NDArray[np.float64]:
        return build_y0(wall_bc=problem.wall_bc, g_wall=problem.g_wall, shooting_vars=s)

    # --------------------------------------------------
    # residual function: f'(inf) = 1, g(inf) = 1
    # --------------------------------------------------
    def residual(y_edge: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.array([
            y_edge[1] - 1.0,  # f'(inf) - 1
            y_edge[3] - 1.0,  # g(inf) - 1
        ])

    # --------------------------------------------------
    # pack into SolverProblem dataclass (convenience wrapper) for the generic shooting method (defined in shooting.py)
    # --------------------------------------------------
    return SolverProblem(
        ode_function=ode_func,
        build_initial_condition=build_initial_condition,
        residual_function=residual,
        n_shooting=2,
        initial_guess=guess,
    )
