"""Tests for mapping converged similarity profiles onto physical HDF5 grids."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from cfd_io import Dataset, StructuredGrid, read_hdf5, write_hdf5
from flow_state import solve
from flow_state.io import write_json
from typer.testing import CliRunner

from simbl.baseflow import BaseflowConfig, map_baseflow_to_grid
from simbl.cli import cli


# --------------------------------------------------
# fixtures
# --------------------------------------------------
def _write_mapping_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Write one converged profile, edge state, and structured grid."""

    # write a canonical converged similarity profile
    profile_path = tmp_path / "simbl.json"
    profile_data = {
        "metadata": {"beta": 0.0},
        "profiles": {
            "eta": [0.0, 1.0, 2.0, 4.0],
            "f": [0.0, 0.2, 1.0, 3.0],
            "fp": [0.0, 0.55, 0.9, 1.0],
            "fpp": [0.55, 0.35, 0.12, 0.0],
            "tau": [1.5, 1.3, 1.1, 1.0],
            "taup": [-0.25, -0.15, -0.05, 0.0],
        },
    }
    profile_path.write_text(json.dumps(profile_data), encoding="utf-8")

    # write dimensional edge conditions with a transport model
    edge_path = tmp_path / "edge_state.json"
    edge_state = solve(mach=3.0, pres=2500.0, temp=220.0)
    write_json(edge_state, edge_path)

    # write a two-station planar grid extending beyond mapped eta_max
    grid_path = tmp_path / "grid.hdf5"
    x_stations = np.array([0.5, 1.0])
    y_points = np.array([0.0, 0.001, 0.003, 0.02, 0.05])
    grid_x = np.broadcast_to(x_stations[:, None, None], (2, 5, 1)).copy()
    grid_y = np.broadcast_to(y_points[None, :, None], (2, 5, 1)).copy()
    grid_z = np.zeros_like(grid_x)
    grid = StructuredGrid(x=grid_x, y=grid_y, z=grid_z)
    write_hdf5(grid_path, Dataset(grid=grid), dtype="f8")

    return profile_path, edge_path, grid_path


# --------------------------------------------------
# mapping tests
# --------------------------------------------------
def test_map_baseflow_preserves_grid_and_fills_farfield(tmp_path: Path) -> None:
    """Mapped profiles should preserve coordinates and use exact edge fill."""
    profile_path, edge_path, grid_path = _write_mapping_inputs(tmp_path)
    output_path = tmp_path / "baseflow.hdf5"

    config = BaseflowConfig(
        solution_input=profile_path,
        edge_state_input=edge_path,
        grid_input=grid_path,
        output=output_path,
    )
    result_path = map_baseflow_to_grid(config)

    source = read_hdf5(grid_path)
    result = read_hdf5(result_path)
    edge_state = solve(mach=3.0, pres=2500.0, temp=220.0)

    assert result_path == output_path
    np.testing.assert_allclose(result.grid.x, source.grid.x)
    np.testing.assert_allclose(result.grid.y, source.grid.y)
    assert set(result.flow) == {"uvel", "vvel", "wvel", "temp", "pres", "dens", "visc"}

    # verify no-slip wall and exact edge-state fill beyond eta_max
    np.testing.assert_allclose(result.flow["uvel"].data[:, 0, :], 0.0)
    np.testing.assert_allclose(result.flow["uvel"].data[:, -1, :], edge_state.uvel)
    np.testing.assert_allclose(result.flow["temp"].data[:, -1, :], edge_state.temp)
    np.testing.assert_allclose(result.flow["pres"].data, edge_state.pres)
    np.testing.assert_allclose(result.flow["dens"].data[:, -1, :], edge_state.dens)
    np.testing.assert_allclose(result.flow["visc"].data[:, -1, :], edge_state.mu)


def test_baseflow_cli_run_uses_relative_config_paths(tmp_path: Path) -> None:
    """Baseflow CLI should resolve artifacts relative to its TOML file."""
    _write_mapping_inputs(tmp_path)
    config_path = tmp_path / "simbl_baseflow.toml"
    config_path.write_text(
        """\
[baseflow]
solution_input = "simbl.json"
edge_state_input = "edge_state.json"
grid_input = "grid.hdf5"
output = "baseflow.hdf5"
similarity_transform = "levy_lees"
geometry_transform = "none"
farfield_fill = "edge"
endpoint_tolerance = 1.0e-3
dtype = "f8"
""",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["baseflow", "run", str(config_path)])

    assert result.exit_code == 0
    assert (tmp_path / "baseflow.hdf5").is_file()
    assert "Wrote:" in result.output


def test_map_baseflow_rejects_unconverged_endpoint(tmp_path: Path) -> None:
    """Edge filling should require convergence at finite eta_max."""
    profile_path, edge_path, grid_path = _write_mapping_inputs(tmp_path)
    profile_data = json.loads(profile_path.read_text(encoding="utf-8"))
    profile_data["profiles"]["fp"][-1] = 0.95
    profile_path.write_text(json.dumps(profile_data), encoding="utf-8")

    config = BaseflowConfig(
        solution_input=profile_path,
        edge_state_input=edge_path,
        grid_input=grid_path,
        output=tmp_path / "baseflow.hdf5",
    )

    with pytest.raises(ValueError, match="has not converged"):
        map_baseflow_to_grid(config)


def test_map_baseflow_rejects_unimplemented_mangler(tmp_path: Path) -> None:
    """Mangler requests should fail rather than silently use planar mapping."""
    profile_path, edge_path, grid_path = _write_mapping_inputs(tmp_path)
    config = BaseflowConfig(
        solution_input=profile_path,
        edge_state_input=edge_path,
        grid_input=grid_path,
        output=tmp_path / "baseflow.hdf5",
        geometry_transform="mangler",
    )

    with pytest.raises(NotImplementedError, match="mangler"):
        map_baseflow_to_grid(config)
