"""Map converged similarity profiles onto structured physical grids."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

from simbl.baseflow.config import BaseflowConfig
from simbl.io import read
from simbl.solver.solution import SimilaritySolution
from simbl.transform import eta2y


# --------------------------------------------------
# public mapping workflow
# --------------------------------------------------
def map_baseflow_to_grid(config: BaseflowConfig) -> Path:
    """Map one converged similarity profile onto a structured HDF5 grid.

    Args:
        config: Validated baseflow mapping configuration.

    Returns:
        Path to the HDF5 dataset containing the original grid and mapped fields.

    Raises:
        ImportError: If cfd-io is unavailable.
        NotImplementedError: If an unsupported geometry or equation family is requested.
        ValueError: If profile, edge-state, or grid data are inconsistent.
    """

    # reject geometry mappings that are not yet implemented
    if config.geometry_transform == "mangler":
        raise NotImplementedError("geometry_transform='mangler' is not implemented; use 'none'")

    # load optional integration packages only for this workflow
    try:
        from cfd_io import Dataset, Field, StructuredGrid, read_hdf5, write_hdf5
    except ImportError as exc:
        raise ImportError(
            "SIMBL baseflow mapping requires cfd-io; install the cfd-io package"
        ) from exc

    from flow_state.io import read_json
    from flow_state.transport import transport_model_from_spec

    # read and validate the similarity solution
    profile = read(config.solution_input)
    _validate_profile_endpoint(profile, config.endpoint_tolerance)
    if profile.g is not None:
        raise NotImplementedError("Falkner-Skan-Cooke baseflow mapping is not implemented")

    # read dimensional edge conditions required by eta2y and dimensionalization
    edge_state = read_json(config.edge_state_input)
    required_edge_values = {
        "uvel": edge_state.uvel,
        "dens": edge_state.dens,
        "visc": edge_state.mu,
    }
    missing_edge_values = [name for name, value in required_edge_values.items() if value is None]
    if missing_edge_values:
        missing = ", ".join(missing_edge_values)
        raise ValueError(f"edge FlowState is missing required values: {missing}")
    if edge_state.transport_model is None:
        raise ValueError("edge FlowState must contain a transport model")

    uvel_edge = float(edge_state.uvel)
    dens_edge = float(edge_state.dens)
    visc_edge = float(edge_state.mu)
    transport = transport_model_from_spec(edge_state.transport_model)

    # read the structured physical grid
    grid_dataset = read_hdf5(config.grid_input)
    if not isinstance(grid_dataset.grid, StructuredGrid):
        raise TypeError("baseflow mapping requires a structured HDF5 grid")
    grid = grid_dataset.grid
    if grid.x.ndim != 3:
        raise ValueError("structured grid coordinates must have shape (ni, nj, nk)")
    if grid.x.shape[1] < 2:
        raise ValueError("structured grid must contain at least two wall-normal points")

    # allocate dimensional fields on the original grid
    field_data = {
        name: np.empty(grid.shape, dtype=np.float64)
        for name in ("uvel", "vvel", "wvel", "temp", "pres", "dens", "visc")
    }

    # map one similarity profile independently onto each wall-normal grid line
    ni, _, nk = grid.shape
    beta = float(profile.metadata.get("beta", 0.0))
    nu_edge = visc_edge / dens_edge
    for streamwise_index in range(ni):
        for spanwise_index in range(nk):
            x_line = grid.x[streamwise_index, :, spanwise_index]
            y_line = grid.y[streamwise_index, :, spanwise_index]
            x_station = float(x_line[0])

            # current planar mapping requires a constant-x wall-normal grid line
            if not np.allclose(x_line, x_station):
                raise ValueError(
                    "geometry_transform='none' requires constant x along each wall-normal grid line"
                )
            if x_station <= 0.0:
                raise ValueError(
                    "eta2y requires every streamwise station x to be greater than zero"
                )

            # compute physical distance from the local wall point
            wall_distance = np.abs(y_line - y_line[0])
            if not np.all(np.diff(wall_distance) > 0.0):
                raise ValueError(
                    "wall-normal grid distance must be strictly increasing from index j=0"
                )

            # map the converged similarity coordinate to physical distance
            source_y = eta2y(
                eta=profile.eta,
                tau=profile.tau,
                x=x_station,
                dens_edge=dens_edge,
                uvel_edge=uvel_edge,
                visc_edge=visc_edge,
                beta=beta,
                transform=config.similarity_transform,
            )

            # build dimensional source profiles
            source_temp = edge_state.temp * profile.tau
            source_dens = dens_edge / profile.tau
            source_visc = np.asarray(
                [transport.mu(float(temp)) for temp in source_temp],
                dtype=np.float64,
            )
            velocity_scale = np.sqrt(nu_edge * uvel_edge / x_station)
            source_vvel = 0.5 * velocity_scale * (profile.eta * profile.fp - profile.f)

            # interpolate inside eta_max and fill the outer grid with edge values
            field_data["uvel"][streamwise_index, :, spanwise_index] = _interpolate_with_edge_fill(
                source_y,
                uvel_edge * profile.fp,
                wall_distance,
                uvel_edge,
            )
            field_data["vvel"][streamwise_index, :, spanwise_index] = _interpolate_with_edge_fill(
                source_y,
                source_vvel,
                wall_distance,
                float(source_vvel[-1]),
            )
            field_data["wvel"][streamwise_index, :, spanwise_index] = 0.0
            field_data["temp"][streamwise_index, :, spanwise_index] = _interpolate_with_edge_fill(
                source_y,
                source_temp,
                wall_distance,
                edge_state.temp,
            )
            field_data["pres"][streamwise_index, :, spanwise_index] = edge_state.pres
            field_data["dens"][streamwise_index, :, spanwise_index] = _interpolate_with_edge_fill(
                source_y,
                source_dens,
                wall_distance,
                dens_edge,
            )
            field_data["visc"][streamwise_index, :, spanwise_index] = _interpolate_with_edge_fill(
                source_y,
                source_visc,
                wall_distance,
                visc_edge,
            )

    # preserve grid metadata and record the mapping inputs
    attrs = dict(grid_dataset.attrs)
    attrs.update(
        {
            "generated by": "simbl.baseflow",
            "simbl solution input": str(config.solution_input),
            "edge state input": str(config.edge_state_input),
            "similarity transform": config.similarity_transform,
            "geometry transform": config.geometry_transform,
            "simbl eta max": float(profile.eta[-1]),
        }
    )
    output_dataset = Dataset(
        grid=grid,
        flow={name: Field(values) for name, values in field_data.items()},
        attrs=attrs,
    )

    # write the mapped dataset and create requested parent directories
    config.output.parent.mkdir(parents=True, exist_ok=True)
    output_path = write_hdf5(config.output, output_dataset, dtype=config.dtype)

    return output_path


# --------------------------------------------------
# validation and interpolation helpers
# --------------------------------------------------
def _validate_profile_endpoint(
    profile: SimilaritySolution,
    tolerance: float,
) -> None:
    """Require the finite eta domain to have reached uniform edge conditions."""

    endpoint_errors = {
        "fp": abs(float(profile.fp[-1]) - 1.0),
        "tau": abs(float(profile.tau[-1]) - 1.0),
        "fpp": abs(float(profile.fpp[-1])),
        "taup": abs(float(profile.taup[-1])),
    }
    failed = {name: error for name, error in endpoint_errors.items() if error > tolerance}
    if failed:
        details = ", ".join(f"{name} error={error:.3e}" for name, error in failed.items())
        raise ValueError(
            "similarity profile has not converged to edge conditions at eta_max: "
            f"{details}; tolerance={tolerance:.3e}"
        )


def _interpolate_with_edge_fill(
    source_y: NDArray[np.float64],
    source_values: NDArray[np.float64],
    target_y: NDArray[np.float64],
    edge_value: float,
) -> NDArray[np.float64]:
    """Cubic-spline interpolate within eta_max and fill beyond it."""

    # validate the physical source coordinate
    if not np.all(np.diff(source_y) > 0.0):
        raise ValueError("eta2y output must be strictly increasing")

    # evaluate only inside the mapped similarity domain
    result = np.full(target_y.shape, edge_value, dtype=np.float64)
    interior = target_y <= source_y[-1]
    spline = CubicSpline(source_y, source_values, extrapolate=False)
    result[interior] = spline(target_y[interior])

    return result
