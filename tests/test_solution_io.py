"""Tests for reading and writing canonical SIMBL solutions."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from simbl.io import read, write
from simbl.solver.falkner_skan.solution import FalknerSkanSolution


# --------------------------------------------------
# solution I/O tests
# --------------------------------------------------
def test_json_solution_round_trip_uses_canonical_profiles(tmp_path: Path) -> None:
    """JSON output should round-trip canonical temperature profiles."""

    # build one representative converged solution
    eta = np.array([0.0, 1.0, 2.0, 4.0])
    solution = FalknerSkanSolution(
        eta=eta,
        f=np.array([0.0, 0.2, 1.0, 3.0]),
        fp=np.array([0.0, 0.55, 0.9, 1.0]),
        fpp=np.array([0.55, 0.35, 0.12, 0.0]),
        tau=np.array([1.5, 1.3, 1.1, 1.0]),
        taup=np.array([-0.25, -0.15, -0.05, 0.0]),
    )
    output_path = tmp_path / "simbl.json"

    # write and inspect the canonical JSON field names
    write(solution, output_path)
    output_data = json.loads(output_path.read_text(encoding="utf-8"))
    assert "tau" in output_data["profiles"]
    assert "taup" in output_data["profiles"]
    assert "g" not in output_data["profiles"]
    assert "gp" not in output_data["profiles"]

    # read through the public extension-based dispatcher
    loaded_solution = read(output_path)
    assert np.array_equal(loaded_solution.tau, solution.tau)
    assert np.array_equal(loaded_solution.taup, solution.taup)

    round_trip_path = tmp_path / "simbl_round_trip.json"
    write(loaded_solution, round_trip_path)
    assert round_trip_path.is_file()


def test_read_rejects_unsupported_extension(tmp_path: Path) -> None:
    """Reader dispatch should reject unregistered input formats."""

    input_path = tmp_path / "simbl.dat"

    with pytest.raises(ValueError, match="Unknown format '.dat'. Supported: .json"):
        read(input_path)
