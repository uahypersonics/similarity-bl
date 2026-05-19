"""Pydantic schemas for TOML configuration validation.

Maps the TOML config format to validated Python objects:

    [conditions]: edge flow conditions (Mach, temperature, beta, gas properties)
    [wall]: wall boundary condition (type, temp_wall)
    [viscosity]: viscosity model selection and parameters
    [numerics]: numerical settings (grid, tolerances, ODE method)
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from flow_state.transport import available_transport_models
from pydantic import BaseModel, Field, field_validator

from simbl.solver.bcs import WallBCType
from simbl.solver.equations import VALID_EQUATIONS


# --------------------------------------------------
# ConditionsConfig: flow conditions at the edge of the boundary layer
# --------------------------------------------------
class ConditionsConfig(BaseModel):
    """Edge flow conditions (Mach, temperature, pressure gradient, gas properties)"""

    # Field(...) with no default means required, omitting these in the TOML raises a ValidationError
    mach_edge: float = Field(..., gt=0, description="Edge Mach number")
    temp_edge: float = Field(..., gt=0, description="Edge temperature [K]")
    # all remaining fields are optional - defaults applied if absent in TOML
    beta: float = Field(default=0.0, description="Falkner-Skan pressure gradient parameter")
    sweep_angle: float = Field(default=0.0, ge=0, lt=90, description="Sweep angle [deg] (0 = 2D Falkner-Skan, >0 = Falkner-Skan-Cooke)")
    gamma: float = Field(default=1.4, gt=1, description="Specific heat ratio")
    prandtl: float = Field(default=0.72, gt=0, description="Prandtl number")


# --------------------------------------------------
# WallConfig: wall boundary condition
# --------------------------------------------------
class WallConfig(BaseModel):
    """Wall boundary condition (adiabatic or isothermal)."""

    # wall BC type: adiabatic or isothermal, default adiabatic if [wall] section is present but type is omitted
    type: WallBCType = Field(default=WallBCType.ADIABATIC, description="Wall BC type")
    # float | None: temp_wall is optional, None is valid for adiabatic, required for isothermal
    # cross-field check enforcing this rule lives in SolverConfig._validate_wall below
    temp_wall: float | None = Field(
        None, gt=0, description="Wall temperature [K] (for isothermal)"
    )


# --------------------------------------------------
# ViscosityConfig: viscosity model selection
# --------------------------------------------------
class ViscosityConfig(BaseModel):
    """Viscosity model selection (model name resolved via flow_state)."""

    model: str = Field(
        default="sutherland", description="Viscosity model (from flow_state)"
    )
    # optional Sutherland overrides - if set, passed directly to get_transport_model()
    # as kwargs, replacing the built-in air() preset values.
    # These match the Sutherland dataclass fields in flow_state:
    # - mu_ref: reference viscosity [Pa s] (air default: 1.716e-5)
    # - T_ref: reference temperature [K] (air default: 273.15)
    # - S: Sutherland constant [K] (air default: 110.4)
    mu_ref: float | None = Field(default=None, gt=0, description="Reference viscosity [Pa s] (Sutherland override)")
    T_ref: float | None = Field(default=None, gt=0, description="Reference temperature [K] (Sutherland override)")
    S: float | None = Field(default=None, gt=0, description="Sutherland constant [K] (Sutherland override)")

    @field_validator("model")
    @classmethod
    def validate_model(cls, v: str) -> str:
        # query flow_state at validation time - rejects unknown model names immediately
        valid = available_transport_models()
        if v not in valid:
            raise ValueError(f"Unknown viscosity model '{v}'. Available: {', '.join(valid)}")
        return v


# --------------------------------------------------
# NumericsConfig: ODE integrator and grid settings
# --------------------------------------------------
class NumericsConfig(BaseModel):
    """Numerical method and grid settings (integrator, grid size, convergence)."""

    eta_max: float = Field(default=15.0, gt=0, description="Maximum similarity coordinate")
    n_points: int = Field(default=500, gt=10, description="Number of grid points")
    tolerance: float = Field(default=1e-8, gt=0, description="Convergence tolerance")
    max_iterations: int = Field(default=100, gt=0, description="Maximum iterations")
    ode_method: str = Field(default="LSODA", description="ODE integrator method")


# --------------------------------------------------
# top-level config
# --------------------------------------------------
class SolverConfig(BaseModel):
    """Unified solver configuration.

        equations = "falkner_skan"

        [conditions]
        mach_edge = 6.0
        temp_edge = 300.0
        beta = 0.0
        gamma = 1.4
        prandtl = 0.72

        [wall]
        type = "isothermal"
        temp_wall = 600.0

        [viscosity]
        model = "sutherland"

        [numerics]
        eta_max = 15.0
        n_points = 500
        tolerance = 1e-8
        max_iterations = 100
        ode_method = "LSODA"
    """

    # equations: which governing ODE system to solve
    # "falkner_skan" (default) or "falkner_skan_cooke" (swept wing)
    # this is the top-level dispatch key passed to solve_similarity(equations=...)
    equations: str = Field(default="falkner_skan", description="Governing equations: falkner_skan or falkner_skan_cooke")
    # conditions is required, no defaults, must be specified in TOML config file
    conditions: ConditionsConfig
    # default_factory=WallConfig constructs a fresh WallConfig() if the [wall] section is absent
    wall: WallConfig = Field(default_factory=WallConfig)
    viscosity: ViscosityConfig = Field(default_factory=ViscosityConfig)
    numerics: NumericsConfig = Field(default_factory=NumericsConfig)

    @field_validator("equations")
    @classmethod
    def _validate_equations(cls, v: str) -> str:
        if v not in VALID_EQUATIONS:
            raise ValueError(f"Unknown equations '{v}'. Choose: {', '.join(sorted(VALID_EQUATIONS))}")
        return v

    @field_validator("wall")
    @classmethod
    def _validate_wall(cls, v: WallConfig) -> WallConfig:
        # cross-field check: isothermal wall requires an explicit wall temperature
        # lives here (on SolverConfig) rather than on WallConfig because it needs
        # both fields of WallConfig visible at once as a complete object
        if v.type == WallBCType.ISOTHERMAL and v.temp_wall is None:
            raise ValueError("temp_wall must be specified for isothermal wall BC")
        return v
