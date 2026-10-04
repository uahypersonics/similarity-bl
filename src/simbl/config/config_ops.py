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
def config_init(
    fname: str | Path,
    *,
    equations: str = "fs",
    force: bool = False,
    flow_state: str | Path | None = None,
) -> list[str]:
    """Write a template configuration file to disk

    Parameters
    ----------
    fname : str or Path
        Output file path.
    equations : str
        Template type: 'fs' (Falkner-Skan, default) or 'fsc' (Falkner-Skan-Cooke).
    force : bool
        If False (default), raises FileExistsError if file already exists.
    flow_state : str or Path, optional
        Canonical FlowState JSON used to initialize edge and gas properties.

    Returns
    -------
    list[str]
        Warnings for optional FlowState properties that could not be imported.

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

    # initialize the selected template from an optional canonical FlowState
    template = templates[equations]
    warnings: list[str] = []
    if flow_state is not None:
        template, warnings = _initialize_template_from_flow_state(
            template,
            flow_state,
        )

    # write and validate the generated configuration
    fname.write_text(template, encoding="utf-8")
    config_load(fname)

    return warnings


# --------------------------------------------------
# FlowState template initialization
# --------------------------------------------------
def _initialize_template_from_flow_state(
    template: str,
    flow_state_path: str | Path,
) -> tuple[str, list[str]]:
    """Populate a SIMBL template from a canonical FlowState JSON file.

    Parameters
    ----------
    template : str
        SIMBL TOML template text.
    flow_state_path : str or Path
        Canonical FlowState JSON path.

    Returns
    -------
    tuple[str, list[str]]
        Updated template text and consistency warnings.

    Raises
    ------
    ValueError
        If the FlowState does not contain required edge properties.
    """
    from flow_state.io import read_json

    # read the canonical edge state through flow-state
    flow_state_path = Path(flow_state_path)
    state = read_json(flow_state_path)

    # validate properties required by the similarity equations
    if state.mach is None:
        raise ValueError("FlowState must contain mach for mach_edge")

    # import the required edge and gas properties
    template = _replace_template_value(template, "mach_edge = 6.0", state.mach)
    template = _replace_template_value(template, "temp_edge = 300.0", state.temp)
    template = _replace_template_value(template, "gamma = 1.4", state.gamma)

    # import optional Prandtl data or retain the documented SIMBL default
    warnings: list[str] = []
    if state.pr is None:
        warnings.append(
            "FlowState does not contain a Prandtl number; retaining SIMBL default prandtl = 0.72"
        )
    else:
        template = _replace_template_value(template, "prandtl = 0.72", state.pr)

    # import the transport model when available
    if state.transport_model is None:
        warnings.append(
            "FlowState does not contain a transport model; retaining SIMBL default model = 'sutherland'"
        )
    else:
        transport_spec = state.transport_model
        model_line = f'model = "{transport_spec.model_type}"'
        model_parameters = transport_spec.parameters
        supported_parameters = {"mu_ref", "T_ref", "S"}
        unsupported_parameters = sorted(set(model_parameters) - supported_parameters)

        # only transfer a complete parameter set that SIMBL can represent
        if unsupported_parameters:
            unsupported = ", ".join(unsupported_parameters)
            warnings.append(
                f"SIMBL cannot import parameters for transport model "
                f"'{transport_spec.model_type}': {unsupported}; using its registered default"
            )
        else:
            parameter_lines = [f"{name} = {value}" for name, value in model_parameters.items()]
            if parameter_lines:
                model_line = "\n".join([model_line, *parameter_lines])

        template = template.replace('model = "sutherland"', model_line, 1)

    # record source provenance without creating a live file dependency
    source_comment = f"# initialized from FlowState: {flow_state_path}"
    template_lines = template.splitlines()
    template_lines.insert(1, source_comment)
    template = "\n".join(template_lines) + "\n"

    return template, warnings


def _replace_template_value(template: str, field: str, value: float) -> str:
    """Replace one controlled scalar assignment in a config template."""
    if template.count(field) != 1:
        raise ValueError(f"SIMBL config template does not contain exactly one {field!r}")

    key = field.split("=", 1)[0].strip()
    replacement = f"{key} = {value}"
    updated_template = template.replace(field, replacement, 1)

    return updated_template


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
        prandtl=config.gas.prandtl,
        gamma=config.gas.gamma,
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
        max_solve_time=config.numerics.max_solve_time,
        equations=config.equations,
    )

    return inputs, options
