"""Main solver entry point

Dispatches to the appropriate model solver (Falkner-Skan or Falkner-Skan-Cooke),
runs the shooting method, saves converged solutions to the lookup table,
and returns the solution.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from typing import TYPE_CHECKING

from flow_state.transport import get_transport_model

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
    problem: SimilarityInputs | None = None,
    options: SolverOptions | None = None,
    *,
    initial_fpp: float | None = None,
    initial_gvar: float | None = None,
    initial_wp: float | None = None,
    **kwargs,
) -> tuple[FalknerSkanSolution | FalknerSkanCookeSolution, ShootingResult]:
    """Solve the compressible similarity boundary layer equations

    This is the main solver entry point. Dispatches to the appropriate
    solver based on the `equations` argument.

    ``problem`` can be omitted and replaced with keyword arguments that are
    forwarded directly to :class:`~simbl.solver.inputs.SimilarityInputs`.
    ``mach`` is accepted as a shorthand for ``mach_edge``.

    Parameters
    ----------
    problem : SimilarityInputs, optional
        Physics specification (Mach, temperature, wall BC, etc.).
        If not provided, ``**kwargs`` are used to construct one.
    options : SolverOptions, optional
        Numerical settings. If None, uses defaults. Also carries the `equations` field.
    initial_fpp : float, optional
        Override initial guess for f''(0).
    initial_gvar : float, optional
        Override initial guess for g variable -- g(0) for adiabatic, g'(0) for isothermal.
    initial_wp : float, optional
        Override initial guess for w'(0) (FSC only).
    **kwargs
        Forwarded to :class:`~simbl.solver.inputs.SimilarityInputs` when ``problem``
        is not provided. ``mach`` is accepted as an alias for ``mach_edge``.

    Returns
    -------
    tuple[FalknerSkanSolution | FalknerSkanCookeSolution, ShootingResult]
        Solution profiles and shooting method convergence info.

    Examples
    --------
    Quickest form — pass fields directly:

    >>> from simbl import solve_similarity
    >>> solution, info = solve_similarity(mach=6.0, temp_edge=55.0)

    Explicit ``SimilarityInputs`` object:

    >>> from simbl import solve_similarity, SimilarityInputs
    >>> problem = SimilarityInputs(mach_edge=4.0, temp_edge=220.0, wall_bc="adiabatic")
    >>> solution, info = solve_similarity(problem)

    Falkner-Skan-Cooke with sweep:

    >>> problem = SimilarityInputs(mach_edge=2.0, temp_edge=220.0, sweep_angle=30.0)
    >>> solution, info = solve_similarity(problem, SolverOptions(equations="falkner_skan_cooke"))
    """
    from simbl.solver.inputs import SimilarityInputs as _SimilarityInputs

    # convenience: build a SimilarityInputs from kwargs if problem not given
    if problem is None:
        if not kwargs:
            raise TypeError("solve_similarity() requires either a 'problem' argument or keyword arguments")
        if "mach" in kwargs:
            kwargs["mach_edge"] = kwargs.pop("mach")
        problem = _SimilarityInputs(**kwargs)

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
    # solve using shooting method in shooting.py
    # result is a ShootingResult dataclass containing:
    # - converged: bool
    # - iterations: int
    # - eta: array of similarity coordinate values
    # - solution: array of solution variables at each eta point
    # - shooting_vars: array of converged shooting variable values at the wall
    # - residual: array of final residual values at the edge
    # --------------------------------------------------
    result = shooting_method(solver_problem, options)

    # --------------------------------------------------
    # save initial guess to lookup table if solution converged
    # this allows building a database of converged solutions for different Mach/temp/wall BCs
    # which can be used as initial guesses for future solves to improve convergence
    # --------------------------------------------------
    if result.converged:

        # load model specific save_converged function (saves to JSON in lookup table directory)
        if equations == "falkner_skan":
            # falkner-skan: initial guess consists of f''(0) and g(0) or g'(0) depending on wall BC type
            from simbl.solver.falkner_skan.initial_guess import save_converged
        else:
            # falkner-skan-cooke: initial guess consists of f''(0), w'(0), and g(0) or g'(0) depending on wall BC type
            from simbl.solver.falkner_skan_cooke.initial_guess import save_converged

        # save the converged solution to the lookup table for future use as an initial guess
        save_converged(result, problem)

    # --------------------------------------------------
    # build model-specific solution dataclass
    # --------------------------------------------------
    return build_solution(result), result
