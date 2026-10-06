"""JSON writer for similarity solutions

Writes similarity profiles and metadata to a JSON file.
Unlike the Tecplot writer, raw similarity variables (f, f', f'', g, g') are
included alongside physical profiles — useful for archiving, post-processing,
and interoperability with other tools.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.falkner_skan.solution import FalknerSkanSolution
    from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.options import SolverOptions
    from simbl.solver.shooting import ShootingResult
    from simbl.solver.solution import SimilaritySolution

from simbl.transform import DEFAULT_ETA2Y_TRANSFORMS


# --------------------------------------------------
# write_json: write similarity profiles and metadata to JSON
# --------------------------------------------------
def _write_json(
    solution: FalknerSkanSolution | FalknerSkanCookeSolution | SimilaritySolution,
    fname: Path,
    problem: SimilarityInputs | None = None,
    solver_options: SolverOptions | None = None,
    shooting_result: ShootingResult | None = None,
) -> None:
    """Write similarity solution to JSON format

    Parameters
    ----------
    solution : FalknerSkanSolution | FalknerSkanCookeSolution
        Similarity solution to write.
    fname : Path
        Output file path.
    problem : SimilarityInputs, optional
        Problem specification for metadata.
    solver_options : SolverOptions, optional
        Numerical solver settings for provenance metadata.
    shooting_result : ShootingResult, optional
        Shooting method convergence info.
    """

    # --------------------------------------------------
    # build metadata: problem inputs, wall values, convergence info
    # mirrors the AUXDATA fields in the Tecplot writer
    # --------------------------------------------------
    metadata: dict = {
        "schema_version": 1,
        "generator": {
            "name": "similarity-bl",
            "version": version("similarity-bl"),
        },
        "generated": datetime.now(UTC).isoformat(),
    }

    # problem inputs (if provided)
    if problem is not None:
        metadata["inputs"] = asdict(problem)

    # numerical settings and transformed-coordinate convention
    if solver_options is not None:
        equations = solver_options.equations
        metadata["solver"] = {
            "equations": equations,
            "configured_method": solver_options.solver_method,
            "ode_method": solver_options.ode_method,
        }
        metadata["profile_transform"] = DEFAULT_ETA2Y_TRANSFORMS[equations]
        metadata["numerics"] = asdict(solver_options)

    # wall values from solution
    metadata["fpp_wall"] = solution.fpp[0]
    metadata["taup_wall"] = solution.taup[0]
    metadata["tau_wall"] = solution.tau[0]

    # write crossflow metadata only when the optional profile pair is populated
    crossflow = getattr(solution, "g", None)
    crossflow_gradient = getattr(solution, "gp", None)
    if crossflow is not None and crossflow_gradient is not None:
        metadata["gp_wall"] = crossflow_gradient[0]

    # convergence info (if provided)
    if shooting_result is not None:
        if "solver" in metadata:
            metadata["solver"]["method"] = shooting_result.method
        metadata["convergence"] = {
            "converged": shooting_result.converged,
            "timed_out": shooting_result.timed_out,
            "iterations": shooting_result.iterations,
            "shooting_variables": shooting_result.shooting_vars.tolist(),
            "residual": shooting_result.residual.tolist(),
            "residual_norm": float(sum(value**2 for value in shooting_result.residual) ** 0.5),
        }

    # bookkeeping
    metadata["n_points"] = len(solution.eta)

    # --------------------------------------------------
    # build profiles: raw similarity variables + crossflow if FSC
    # note: .tolist() converts numpy arrays to plain Python lists (JSON cannot serialize numpy types)
    # --------------------------------------------------
    profiles: dict[str, list[float]] = {
        "eta": solution.eta.tolist(),
        "f": solution.f.tolist(),
        "fp": solution.fp.tolist(),
        "fpp": solution.fpp.tolist(),
        "tau": solution.tau.tolist(),
        "taup": solution.taup.tolist(),
    }

    # crossflow profiles (FSC only)
    if crossflow is not None and crossflow_gradient is not None:
        profiles["g"] = crossflow.tolist()
        profiles["gp"] = crossflow_gradient.tolist()

    # --------------------------------------------------
    # combine metadata and profile dictionaries to output dictionary
    # --------------------------------------------------
    output = {
        "metadata": metadata,
        "profiles": profiles,
    }

    # --------------------------------------------------
    # write output dictionary to JSON file with indentation for readability
    # --------------------------------------------------
    with open(fname, "w") as f:
        json.dump(output, f, indent=2)
