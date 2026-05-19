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
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.falkner_skan.solution import FalknerSkanSolution
    from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.shooting import ShootingResult


# --------------------------------------------------
# write_json: write similarity profiles and metadata to JSON
# --------------------------------------------------
def _write_json(
    solution: FalknerSkanSolution | FalknerSkanCookeSolution,
    fname: Path,
    problem: SimilarityInputs | None = None,
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
    shooting_result : ShootingResult, optional
        Shooting method convergence info.
    """

    # --------------------------------------------------
    # build metadata: problem inputs, wall values, convergence info
    # mirrors the AUXDATA fields in the Tecplot writer
    # --------------------------------------------------
    metadata: dict = {}

    # problem inputs (if provided)
    if problem is not None:
        metadata["mach_edge"] = problem.mach_edge
        metadata["temp_edge"] = problem.temp_edge
        metadata["prandtl"] = problem.prandtl
        metadata["gamma"] = problem.gamma
        metadata["beta"] = problem.beta
        metadata["wall_bc_type"] = problem.wall_bc
        if problem.sweep_angle != 0.0:
            metadata["sweep_angle"] = problem.sweep_angle

    # wall values from solution
    metadata["fpp_wall"] = solution.fpp[0]
    metadata["gp_wall"] = solution.gp[0]
    metadata["g_wall"] = solution.g[0]
    if hasattr(solution, "wp"):
        metadata["wp_wall"] = solution.wp[0]

    # convergence info (if provided)
    if shooting_result is not None:
        metadata["converged"] = shooting_result.converged
        metadata["iterations"] = shooting_result.iterations

    # bookkeeping
    metadata["n_points"] = len(solution.eta)
    metadata["generated"] = datetime.now().isoformat()

    # --------------------------------------------------
    # build profiles: raw similarity variables + crossflow if FSC
    # note: .tolist() converts numpy arrays to plain Python lists (JSON cannot serialize numpy types)
    # --------------------------------------------------
    profiles: dict[str, list[float]] = {
        "eta": solution.eta.tolist(),
        "f": solution.f.tolist(),
        "fp": solution.fp.tolist(),
        "fpp": solution.fpp.tolist(),
        "g": solution.g.tolist(),
        "gp": solution.gp.tolist(),
    }

    # crossflow profiles (FSC only)
    if hasattr(solution, "w"):
        profiles["w"] = solution.w.tolist()
        profiles["wp"] = solution.wp.tolist()

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
