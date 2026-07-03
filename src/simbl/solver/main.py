"""Main solver entry point

Dispatches to the appropriate model solver (Falkner-Skan or Falkner-Skan-Cooke),
runs the shooting or the bvp method, saves converged solutions to the lookup table,
and returns the solution.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from flow_state.transport import get_transport_model

from simbl.solver.bvp import bvp_method
from simbl.solver.falkner_skan.solution import FalknerSkanSolution
from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
from simbl.solver.options import SolverOptions
from simbl.solver.shooting import ShootingResult, shooting_method

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.inputs import SimilarityInputs


# --------------------------------------------------
# solve_similarity
# --------------------------------------------------
def solve_similarity(
    problem: SimilarityInputs,
    options: SolverOptions | None = None,
    *,
    initial_fpp: float | None = None,
    initial_gvar: float | None = None,
    initial_wp: float | None = None,
    initial_profile: ShootingResult | None = None,
) -> tuple[FalknerSkanSolution | FalknerSkanCookeSolution, ShootingResult]:
    """Solve the compressible similarity boundary layer equations

    This is the main solver entry point. Dispatches to the appropriate
    solver based on the `equations` argument.

    ``problem`` can be omitted and replaced with keyword arguments that are
    forwarded directly to :class:`~simbl.solver.inputs.SimilarityInputs`.
    ``mach`` is accepted as a shorthand for ``mach_edge``.

    Parameters
    ----------
    problem : SimilarityInputs
        Physics specification (Mach, temperature, wall BC, etc.).
    options : SolverOptions, optional
        Numerical settings. If None, uses defaults. Also carries the `equations` field.
    initial_fpp : float, optional
        Override initial guess for f''(0).
    initial_gvar : float, optional
        Override initial guess for g variable -- g(0) for adiabatic, g'(0) for isothermal.
    initial_wp : float, optional
        Override initial guess for w'(0) (FSC only).
    initial_profile : ShootingResult, optional
        Full solution profile from a nearby converged solve (e.g. the previous
        Mach point).  When ``solver_method="shooting_then_bvp"`` and shooting
        fails, this profile is passed to the BVP fallback as the initial mesh
        guess instead of the failed shooting result.  A physically correct profile
        from an adjacent parameter point gives the collocation solver a far better
        starting state than the diverged shooting trajectory.

    Returns
    -------
    tuple[FalknerSkanSolution | FalknerSkanCookeSolution, ShootingResult]
        Solution profiles and shooting method convergence info.

    Examples
    --------
    Quickest form — pass fields directly:

    >>> from simbl import solve_similarity, SimilarityInputs
    >>> solution, info = solve_similarity(SimilarityInputs(mach_edge=6.0, temp_edge=55.0))

    Explicit ``SimilarityInputs`` object:

    >>> from simbl import solve_similarity, SimilarityInputs
    >>> problem = SimilarityInputs(mach_edge=4.0, temp_edge=220.0, wall_bc="adiabatic")
    >>> solution, info = solve_similarity(problem)

    Falkner-Skan-Cooke with sweep:

    >>> problem = SimilarityInputs(mach_edge=2.0, temp_edge=220.0, sweep_angle=30.0)
    >>> solution, info = solve_similarity(problem, SolverOptions(equations="falkner_skan_cooke"))
    """
    # set default options if not provided
    if options is None:
        options = SolverOptions()

    equations = options.equations

    # catch mismatched equations/sweep_angle before doing any work
    if equations == "falkner_skan" and problem.sweep_angle != 0.0:
        raise ValueError(
            f"sweep_angle={problem.sweep_angle} is set but equations='falkner_skan' (2D). "
            "Use equations='falkner_skan_cooke' for swept-wing problems."
        )
    # FSC at sweep_angle == 0 is allowed: the crossflow equation decouples
    # (tan^2(Lambda) = 0 in the energy equation, w solves an independent
    # linear BVP), and the streamwise/energy quantities recover Falkner-Skan.

    # create viscosity model (shared across all models)
    # viscosity_model_kwargs allows overriding built-in air() preset constants
    visc_model = get_transport_model(problem.viscosity_model, **(problem.viscosity_model_kwargs or {}))

    # --------------------------------------------------
    # dispatch to model-specific builder
    # --------------------------------------------------
    if equations == "falkner_skan":

        # import model-specific builder functions
        from simbl.solver.falkner_skan.builder import build_solver_problem
        from simbl.solver.falkner_skan.solution import build_solution

        # build the model-specific solver problem dataclass, which includes:
        # - the ODE system
        # - the initial condition builder
        # - the residual function
        # - the number of shooting variables
        # - the initial guess
        #
        # the problem is the SimilarityInputs
        # options are the SolverOptions
        # visc_model is the transport model instance (e.g. Sutherland) coming from flow_state
        # initial_fpp and initial_gvar are optional overrides for the initial guess
        solver_problem = build_solver_problem(
            problem, visc_model,
            initial_fpp=initial_fpp,
            initial_gvar=initial_gvar,
        )

    elif equations == "falkner_skan_cooke":

        # import model-specific builder functions
        from simbl.solver.falkner_skan_cooke.builder import build_solver_problem
        from simbl.solver.falkner_skan_cooke.solution import build_solution

        # same as FS above but with w'(0) as additional shooting variable
        solver_problem = build_solver_problem(
            problem, visc_model,
            initial_fpp=initial_fpp,
            initial_wp=initial_wp,
            initial_gvar=initial_gvar,
        )

    else:

        # unknown model => raise error and list valid options
        raise ValueError(
            f"Unknown equations: {equations!r}. "
            "Choose 'falkner_skan' or 'falkner_skan_cooke'."
        )

    # --------------------------------------------------
    # dispatch to solver
    #
    # "bvp"              : scipy.integrate.solve_bvp only (robust, slower)
    # "shooting"         : Newton shooting only (default, fast)
    # --------------------------------------------------

    if options.solver_method == "bvp":
        # BVP solver

        # run shooting first to obtain a physically shaped initial profile
        shooting_result = shooting_method(solver_problem, options)

        # if the shooting converged, use its solution as the initial profile for the BVP solver
        shooting_profile = shooting_result if not np.any(np.isnan(shooting_result.solution)) else None

        # run bvp solver
        result = bvp_method(solver_problem, options, initial_profile=shooting_profile, problem=problem)

    else:

        # shooting solver (default)
        result = shooting_method(solver_problem, options)

        # BVP fallback: if shooting failed and the flag is set, retry with bvp method
        if not result.converged and options.bvp_fallback:
            if initial_profile is not None:
                bvp_initial = initial_profile
            elif not np.any(np.isnan(result.solution)):
                bvp_initial = result
            else:
                bvp_initial = None
            result = bvp_method(solver_problem, options, initial_profile=bvp_initial, problem=problem)

    # --------------------------------------------------
    # build model-specific solution dataclass
    # --------------------------------------------------
    return build_solution(result), result
