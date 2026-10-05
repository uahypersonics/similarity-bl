"""Reader registry and dispatcher for solution formats."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from simbl.io.json_reader import _read_json
from simbl.solver.solution import SimilaritySolution

# --------------------------------------------------
# reader registry
# --------------------------------------------------
_READERS: dict[str, Callable[[Path], SimilaritySolution]] = {
    "json": _read_json,
}


# --------------------------------------------------
# public API
# --------------------------------------------------
def get_supported_formats() -> list[str]:
    """Return the sorted supported input file extensions."""

    supported_formats = sorted(_READERS.keys())
    return supported_formats


def read(path: str | Path) -> SimilaritySolution:
    """Read a solution file with its format inferred from the extension.

    Args:
        path: Solution file path.

    Returns:
        Validated similarity solution.

    Raises:
        ValueError: If the file extension is unsupported.
    """

    # convert to Path object
    path = Path(path)

    # extract and validate the input format
    extension = path.suffix.lower().lstrip(".")
    if extension not in _READERS:
        supported = ", ".join(f".{item}" for item in get_supported_formats())
        raise ValueError(f"Unknown format '.{extension}'. Supported: {supported}")

    # dispatch to the registered reader
    reader = _READERS[extension]
    solution = reader(path)
    return solution
