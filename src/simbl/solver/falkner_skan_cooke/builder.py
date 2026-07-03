"""Falkner-Skan-Cooke problem builder

Constructs a SolverProblem for the 3D compressible Falkner-Skan-Cooke model
(7-equation ODE, 3 shooting variables).

Implements Liu (2021), Phys. Fluids 33, 126109, Eqs. 16-18.
State vector: [f, f', f'', tau, tau', g, g']
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
        Override initial guess for g_cf'(0) (crossflow velocity gradient at wall).
    initial_gvar : float, optional
        Override initial guess for tau'(0) (isothermal) or tau(0) (adiabatic).

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
        if problem.wall_bc == "adiabatic":
            key_values["temp_edge"] = problem.temp_edge
        elif problem.wall_bc == "isothermal":
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

        # for adiabatic walls, override the lookup-table g_wall estimate with the
        # analytical recovery temperature formula.  The shipped table was generated at
        # T_edge=300 K; at other T_edge values the IDW prediction is inaccurate because
        # Sutherland viscosity is non-linear in absolute temperature.
        # The recovery formula is a reliable starting point for all Mach / T_edge combos.
        if problem.wall_bc == "adiabatic":
            # recovery factor r ~ Pr^(1/3) for turbulent; use same value for laminar FS
            r = problem.prandtl ** (1.0 / 3.0)
            g_wall_est = 1.0 + r * (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2
            guess[2] = g_wall_est

    # --------------------------------------------------
    # ODE setup: precompute loop-invariant parameters S and K, then bind
    # everything into a two-argument function for solve_ivp
    # --------------------------------------------------

    # S (Liu Eq. 22) and K (Liu Eq. 21) depend only on problem constants.
    # The solver always receives local edge conditions (chi = xi/xi_ref = 1),
    # so S simplifies to: S = 1 + (gamma-1)/2 * Ma_e^2 * cos^2(Lambda)
    # and K = [1 + (gamma-1)/2 * Ma_e^2] / S
    # Computing them once here avoids repeating trig + arithmetic on every ODE call.
    sweep_rad = np.radians(problem.sweep_angle)
    cos2_lambda = np.cos(sweep_rad) ** 2
    S = 1.0 + (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2 * cos2_lambda
    K = (1.0 + (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2) / S

    def ode_func(eta: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return bl_ode(eta, y, problem, visc_model, S, K)

    # --------------------------------------------------
    # Initial condition builder (3 shooting variables -> 7-element y0)
    # State vector: [f, f', f'', tau, tau', g, g']
    # Note: problem.g_wall = T_wall/T_edge is passed as tau_wall to build_y0
    # --------------------------------------------------
    def build_initial_condition(s: NDArray[np.float64]) -> NDArray[np.float64]:
        return build_y0(wall_bc=problem.wall_bc, tau_wall=problem.g_wall, shooting_vars=s)

    # --------------------------------------------------
    # Residual function for boundary conditions at eta -> infinity:
    #   fp(inf) = 1     (streamwise velocity reaches edge value)
    #   g_cf(inf) = g_cf_edge_target
    #   tau(inf) = 1    (temperature reaches edge value)
    #
    # For aligned flow (sweep_angle = 0), w_e = 0 so g_cf(inf) = 0
    # For swept flow (sweep_angle > 0), g_cf = w/w_e is normalized by edge crossflow so g_cf(inf) = 1
    # --------------------------------------------------
    gcf_edge_target = 0.0 if problem.sweep_angle == 0.0 else 1.0

    def residual(y_edge: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.array([
            y_edge[1] - 1.0,              # f'(inf) - 1
            y_edge[5] - gcf_edge_target,  # g(inf) - g_edge_target
            y_edge[3] - 1.0,              # tau(inf) - 1
        ])

    # --------------------------------------------------
    # boundary condition function for solve_bvp (collocation fallback)
    #
    # state vector: [f, f', f'', tau, tau', g, g']  (indices 0-6)
    # wall: f=0, f'=0, g=0 always; tau or tau' fixed depending on wall_bc
    # edge: f'=1, tau=1, g=gcf_edge_target
    # --------------------------------------------------
    if problem.wall_bc == "isothermal":
        tau_wall_val = problem.g_wall

        def bc_func(ya: NDArray[np.float64], yb: NDArray[np.float64]) -> NDArray[np.float64]:
            return np.array([
                ya[0],                        # f(0) = 0
                ya[1],                        # f'(0) = 0
                ya[3] - tau_wall_val,          # tau(0) = Tw/Te
                ya[5],                         # g(0) = 0
                yb[1] - 1.0,                   # f'(inf) = 1
                yb[3] - 1.0,                   # tau(inf) = 1
                yb[5] - gcf_edge_target,       # g(inf) = gcf_edge_target
            ])

        # shooting unknowns: f''(0), g'(0), tau'(0) — indices 2, 4, 6
        sv_indices = [2, 4, 6]

    else:
        # adiabatic: tau'(0) = 0; tau(0) is the shooting unknown
        def bc_func(ya: NDArray[np.float64], yb: NDArray[np.float64]) -> NDArray[np.float64]:
            return np.array([
                ya[0],                    # f(0) = 0
                ya[1],                    # f'(0) = 0
                ya[4],                    # tau'(0) = 0  (zero heat flux)
                ya[5],                    # g(0) = 0
                yb[1] - 1.0,              # f'(inf) = 1
                yb[3] - 1.0,              # tau(inf) = 1
                yb[5] - gcf_edge_target,  # g(inf) = gcf_edge_target
            ])

        # shooting unknowns: f''(0), g'(0), tau(0) — indices 2, 3, 6
        sv_indices = [2, 3, 6]

    # --------------------------------------------------
    # pack into SolverProblem dataclass (convenience wrapper) for the generic shooting method (defined in shooting.py)
    # --------------------------------------------------
    return SolverProblem(
        ode_function=ode_func,
        build_initial_condition=build_initial_condition,
        residual_function=residual,
        n_shooting=3,
        initial_guess=guess,
        bc_function=bc_func,
        shooting_var_indices=sv_indices,
    )
