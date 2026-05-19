"""Config file operations and solver input conversion.

Three public file operations:
    config_load      - load and validate a TOML config file -> SolverConfig
    config_save      - serialize a SolverConfig back to a TOML file
    config_init      - write a blank template config file

One public conversion:
    config_to_inputs - convert SolverConfig -> SimilarityInputs + SolverOptions (equations lives on SolverOptions)
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
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.options import SolverOptions

import tomllib

from simbl.config.schema import SolverConfig
from simbl.io.toml_writer import _write_toml


# --------------------------------------------------
# config_load: read and validate TOML -> SolverConfig
# --------------------------------------------------
def config_load(fname: str | Path) -> SolverConfig:
    """Load and validate a TOML configuration file

    Parameters
    ----------
    fname : str or Path
        Path to the TOML configuration file.

    Returns
    -------
    SolverConfig
        Validated solver configuration.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If validation fails.
    """

    fname = Path(fname)

    # data is a plain dictionary
    with open(fname, "rb") as f:
        data = tomllib.load(f)

    # validate and convert to SolverConfig
    return SolverConfig(**data)


# --------------------------------------------------
# config_save: serialize SolverConfig -> TOML file
# --------------------------------------------------
def config_save(config: SolverConfig, fname: str | Path) -> None:
    """Write a SolverConfig to a TOML file

    Parameters
    ----------
    config : SolverConfig
        Validated solver configuration to save.
    fname : str or Path
        Output file path.
    """

    # convert fname to Path object
    fname = Path(fname)
    # ensure parent directory exists
    fname.parent.mkdir(parents=True, exist_ok=True)

    # serialize to plain dict, stripping None values (TOML cannot represent None)
    # note: model_dump is a pydantic method that converts the object back to a plain dict
    # exclude_none=True omits any fields with value None (which could not be serialized to TOML)
    data = config.model_dump(exclude_none=True)

    # pydantic serializes enums as their value, but we need the string
    # e.g. WallBCType.ADIABATIC -> "adiabatic"
    if "wall" in data and "type" in data["wall"]:
        data["wall"]["type"] = config.wall.type.value

    _write_toml(data, fname)


# --------------------------------------------------
# config_init: write blank template config file
# --------------------------------------------------
def config_init(fname: str | Path, *, equations: str = "fs", force: bool = False) -> None:
    """Write a template configuration file to disk

    Parameters
    ----------
    fname : str or Path
        Output file path.
    equations : str
        Template type: 'fs' (Falkner-Skan, default) or 'fsc' (Falkner-Skan-Cooke).
    force : bool
        If False (default), raises FileExistsError if file already exists.

    Raises
    ------
    FileExistsError
        If file exists and force=False.
    ValueError
        If equations is not 'fs' or 'fsc'.
    """
    fname = Path(fname)
    if fname.exists() and not force:
        raise FileExistsError(f"File already exists: {fname}. Use force=True to overwrite.")

    # select template based on equations type
    # fs  = 2D Falkner-Skan (flat plate / wedge)
    # fsc = 3D Falkner-Skan-Cooke (swept wing)
    from simbl.config.template import CONFIG_TEMPLATE_FS, CONFIG_TEMPLATE_FSC
    templates = {"fs": CONFIG_TEMPLATE_FS, "fsc": CONFIG_TEMPLATE_FSC}
    if equations not in templates:
        raise ValueError(f"Unknown equations '{equations}'. Choose 'fs' or 'fsc'.")

    fname.write_text(templates[equations])


# --------------------------------------------------
# config_to_inputs: convert SolverConfig -> SimilarityInputs + SolverOptions
# --------------------------------------------------
def config_to_inputs(config: SolverConfig) -> tuple[SimilarityInputs, SolverOptions]:
    """Convert a SolverConfig into SimilarityInputs and SolverOptions

    Parameters
    ----------
    config : SolverConfig
        Validated solver configuration.

    Returns
    -------
    tuple[SimilarityInputs, SolverOptions]
        Problem specification and solver options.
        The equations field on SolverOptions carries config.equations.
    """
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.options import SolverOptions

    # collect non-None Sutherland overrides to forward to get_transport_model()
    visc_kwargs: dict[str, float] = {}
    if config.viscosity.mu_ref is not None:
        visc_kwargs["mu_ref"] = config.viscosity.mu_ref
    if config.viscosity.T_ref is not None:
        visc_kwargs["T_ref"] = config.viscosity.T_ref
    if config.viscosity.S is not None:
        visc_kwargs["S"] = config.viscosity.S

    inputs = SimilarityInputs(
        mach_edge=config.conditions.mach_edge,
        temp_edge=config.conditions.temp_edge,
        wall_bc=config.wall.type.value,
        temp_wall=config.wall.temp_wall,
        prandtl=config.conditions.prandtl,
        gamma=config.conditions.gamma,
        beta=config.conditions.beta,
        sweep_angle=config.conditions.sweep_angle,
        viscosity_model=config.viscosity.model,
        viscosity_model_kwargs=visc_kwargs or None,
    )

    options = SolverOptions(
        eta_max=config.numerics.eta_max,
        n_points=config.numerics.n_points,
        tolerance=config.numerics.tolerance,
        max_iterations=config.numerics.max_iterations,
        ode_method=config.numerics.ode_method,
        equations=config.equations,
    )

    return inputs, options
