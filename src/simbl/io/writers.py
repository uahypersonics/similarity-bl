"""Writer registry and dispatcher for output formats

Provides a unified write() interface that auto-detects format from
file extension and dispatches to the appropriate writer function.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.falkner_skan.solution import FalknerSkanSolution
    from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.main import ShootingResult

from simbl.io.json_writer import _write_json
from simbl.io.tecplot_writer import _write_tecplot

# --------------------------------------------------
# writer registry: maps file extension --> writer function
# --------------------------------------------------
_WRITERS = {
    "dat": _write_tecplot,
    "json": _write_json,
}


# --------------------------------------------------
# get_supported_formats: list registered file extensions
# --------------------------------------------------
def get_supported_formats() -> list[str]:
    """Return sorted list of supported file extensions"""
    return sorted(_WRITERS.keys())


# --------------------------------------------------
# write: auto-detect format from extension and dispatch to writer
# --------------------------------------------------
def write(
    solution: FalknerSkanSolution | FalknerSkanCookeSolution,
    fname: str | Path,
    *,
    problem: SimilarityInputs | None = None,
    shooting_result: ShootingResult | None = None,
) -> Path:
    """Write solution to file, format auto-detected from extension

    Parameters
    ----------
    solution : FalknerSkanSolution | FalknerSkanCookeSolution
        Solution to write.
    fname : str or Path
        Output file path. Format determined by extension.
    problem : SimilarityInputs, optional
        Problem specification for metadata.
    shooting_result : ShootingResult, optional
        Convergence info for metadata.

    Returns
    -------
    Path
        Path to the written file.

    Raises
    ------
    ValueError
        If the file extension is not supported.
    """

    # convert fname to Path object
    fname = Path(fname)
    # ensure parent directory exists
    fname.parent.mkdir(parents=True, exist_ok=True)

    # extract file extension (lowercase, no dot)
    ext = fname.suffix.lower().lstrip(".")

    # check if extension is supported
    if ext not in _WRITERS:
        supported = ", ".join(f".{e}" for e in get_supported_formats())
        raise ValueError(f"Unknown format '.{ext}'. Supported: {supported}")

    # dispatch to the registered writer function
    writer = _WRITERS[ext]
    writer(solution, fname, problem=problem, shooting_result=shooting_result)

    return fname
