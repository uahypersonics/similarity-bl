"""Configuration contract for mapping similarity profiles onto physical grids."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# --------------------------------------------------
# baseflow config
# --------------------------------------------------
@dataclass(frozen=True, slots=True)
class BaseflowConfig:
    """Validated inputs for one physical-grid baseflow mapping."""

    solution_input: Path
    edge_state_input: Path
    grid_input: Path
    output: Path
    similarity_transform: str = "levy_lees"
    geometry_transform: str = "none"
    farfield_fill: str = "edge"
    endpoint_tolerance: float = 1.0e-3
    dtype: str = "f8"

    def __post_init__(self) -> None:
        """Validate mapping options and file formats."""

        # validate supported input and output formats
        json_paths = {
            "solution_input": self.solution_input,
            "edge_state_input": self.edge_state_input,
        }
        for field_name, path in json_paths.items():
            if path.suffix.lower() != ".json":
                raise ValueError(f"[baseflow].{field_name} must use a .json file")

        hdf5_paths = {"grid_input": self.grid_input, "output": self.output}
        for field_name, path in hdf5_paths.items():
            if path.suffix.lower() not in {".h5", ".hdf5"}:
                raise ValueError(f"[baseflow].{field_name} must use an HDF5 file")

        # validate transform and interpolation policies
        valid_similarity_transforms = {"levy_lees", "illingworth_stewartson"}
        if self.similarity_transform not in valid_similarity_transforms:
            valid = ", ".join(sorted(valid_similarity_transforms))
            raise ValueError(f"[baseflow].similarity_transform must be one of: {valid}")
        if self.geometry_transform not in {"none", "mangler"}:
            raise ValueError("[baseflow].geometry_transform must be 'none' or 'mangler'")
        if self.farfield_fill != "edge":
            raise ValueError("[baseflow].farfield_fill currently supports only 'edge'")
        if self.endpoint_tolerance <= 0.0:
            raise ValueError("[baseflow].endpoint_tolerance must be greater than zero")
        if self.dtype not in {"f4", "f8"}:
            raise ValueError("[baseflow].dtype must be 'f4' or 'f8'")


# --------------------------------------------------
# config parser
# --------------------------------------------------
def load_baseflow_config(path: str | Path) -> BaseflowConfig:
    """Read and validate a focused ``[baseflow]`` TOML config.

    Args:
        path: Baseflow config path.

    Returns:
        Validated baseflow mapping config.
    """

    # read the focused TOML document
    path = Path(path)
    with path.open("rb") as stream:
        data = tomllib.load(stream)

    section = data.get("baseflow")
    if not isinstance(section, dict):
        raise TypeError("config must contain a [baseflow] table")

    # reject misspelled or unsupported fields
    allowed_keys = {
        "solution_input",
        "edge_state_input",
        "grid_input",
        "output",
        "similarity_transform",
        "geometry_transform",
        "farfield_fill",
        "endpoint_tolerance",
        "dtype",
    }
    unknown_keys = sorted(set(section) - allowed_keys)
    if unknown_keys:
        unknown = ", ".join(unknown_keys)
        raise ValueError(f"unknown [baseflow] fields: {unknown}")

    # resolve artifact paths relative to the config location
    config_directory = path.parent
    config = BaseflowConfig(
        solution_input=_required_path(section, "solution_input", config_directory),
        edge_state_input=_required_path(section, "edge_state_input", config_directory),
        grid_input=_required_path(section, "grid_input", config_directory),
        output=_required_path(section, "output", config_directory),
        similarity_transform=str(section.get("similarity_transform", "levy_lees")),
        geometry_transform=str(section.get("geometry_transform", "none")),
        farfield_fill=str(section.get("farfield_fill", "edge")),
        endpoint_tolerance=float(section.get("endpoint_tolerance", 1.0e-3)),
        dtype=str(section.get("dtype", "f8")),
    )

    return config


def _required_path(section: dict[str, Any], name: str, directory: Path) -> Path:
    """Read one required path relative to the config directory."""
    value = section.get(name)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"[baseflow].{name} must be a nonempty path string")

    path = Path(value)
    if not path.is_absolute():
        path = directory / path

    return path


# --------------------------------------------------
# starter config
# --------------------------------------------------
BASEFLOW_CONFIG_TEMPLATE = """\
# Generated by: simbl baseflow init
# Edit the section below, then run: simbl baseflow run

# --------------------------------------------------
# physical-grid baseflow mapping
# --------------------------------------------------

[baseflow]
# converged SIMBL JSON solution
solution_input = "simbl.json"

# dimensional edge FlowState JSON
edge_state_input = "edge_state.json"

# structured HDF5 grid generated by grid-generator
grid_input = "grid.hdf5"

# HDF5 grid plus mapped dimensional flow fields
output = "baseflow.hdf5"

# similarity-coordinate inversion: levy_lees or illingworth_stewartson
similarity_transform = "levy_lees"

# optional geometry mapping: none or mangler
geometry_transform = "none"

# values outside the converged profile use exact edge conditions
farfield_fill = "edge"

# maximum endpoint residual accepted before edge filling
endpoint_tolerance = 1.0e-3

# HDF5 floating-point storage: f4 or f8
dtype = "f8"
"""
