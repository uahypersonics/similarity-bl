"""JSON reader for canonical SIMBL similarity solutions."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from simbl.solver.solution import SimilaritySolution


# --------------------------------------------------
# JSON reader
# --------------------------------------------------
def _read_json(path: Path) -> SimilaritySolution:
    """Read a canonical SIMBL JSON solution.

    Args:
        path: SIMBL JSON solution path.

    Returns:
        Validated similarity solution.

    Raises:
        TypeError: If the JSON document structure is invalid.
        ValueError: If required profiles are missing or invalid.
    """

    # read the JSON document
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError("SIMBL JSON root must be an object")

    # extract controlled metadata and profile mappings
    metadata = data.get("metadata", {})
    profiles = data.get("profiles")
    if not isinstance(metadata, dict):
        raise TypeError("SIMBL JSON metadata must be an object")
    if not isinstance(profiles, dict):
        raise TypeError("SIMBL JSON profiles must be an object")

    # require the canonical profile variables produced by SIMBL solvers
    required_names = ("eta", "f", "fp", "fpp", "tau", "taup")
    missing_names = [name for name in required_names if name not in profiles]
    if missing_names:
        missing = ", ".join(missing_names)
        raise ValueError(f"SIMBL JSON is missing profiles: {missing}")

    # convert lists to explicit float64 arrays
    arrays = {name: np.asarray(profiles[name], dtype=np.float64) for name in required_names}
    crossflow = {
        name: np.asarray(profiles[name], dtype=np.float64)
        for name in ("g", "gp")
        if name in profiles
    }
    solution = SimilaritySolution(
        eta=arrays["eta"],
        f=arrays["f"],
        fp=arrays["fp"],
        fpp=arrays["fpp"],
        tau=arrays["tau"],
        taup=arrays["taup"],
        metadata=dict(metadata),
        g=crossflow.get("g"),
        gp=crossflow.get("gp"),
    )

    return solution
