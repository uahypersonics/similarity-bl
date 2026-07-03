"""Physics specifications required to define a similarity boundary layer problem

SimilarityInputs class:
- contains parameters describing the physical problem: edge Mach number, edge temperature, wall boundary condition, gas properties, etc.
- validation logic to ensure inputs are physically meaningful (e.g. positive Mach and temperature)
- immutability to ensure problem definition cannot be accidentally modified after creation
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

from simbl.solver.bcs import VALID_WALL_BC_TYPES


# --------------------------------------------------
# SimilarityInputs dataclass (immutable)
# --------------------------------------------------
@dataclass(frozen=True)
class SimilarityInputs:
    """Inputs for a similarity boundary layer problem independent of numerical solver settings

    Parameters
    ----------
    mach_edge : float
        Edge Mach number
    temp_edge : float
        Edge temperature [K]
    wall_bc : str
        Wall boundary condition type: "isothermal" or "adiabatic".
    temp_wall : float, optional
        Wall temperature [K] for isothermal walls. Default None (adiabatic).
    prandtl : float, optional
        Prandtl number. Default 0.72 (air).
    gamma : float, optional
        Specific heat ratio. Default 1.4 (air).
    beta : float, optional
        Hartree pressure gradient parameter beta = 2m/(m+1). Default 0.0 (flat plate).
    sweep_angle : float, optional
        Sweep angle [deg] for Falkner-Skan-Cooke problems. Default 0.0 (2D).
    viscosity_model : str, optional
        Transport model for dynamic viscosity. Default "sutherland".

    Examples
    --------
    Adiabatic flat plate at Mach 4:

    >>> from simbl import SimilarityInputs
    >>> inputs = SimilarityInputs(
    ...     mach_edge=4.0,
    ...     temp_edge=300.0,
    ...     wall_bc="adiabatic",
    ... )

    Isothermal wall (cold wall, T_wall = 300 K):

    >>> inputs = SimilarityInputs(
    ...     mach_edge=6.0,
    ...     temp_edge=55.0,
    ...     temp_wall=300.0,
    ...     wall_bc="isothermal",
    ... )
    """

    # --------------------------------------------------
    # define parameters with type annotations and default values
    # --------------------------------------------------

    # edge conditions
    mach_edge: float
    temp_edge: float

    # wall boundary condition

    # "isothermal" or "adiabatic"
    wall_bc: str = "adiabatic"
    # wall temperature [K] for isothermal walls; None for adiabatic walls
    temp_wall: float | None = None

    # gas properties
    prandtl: float = 0.72
    gamma: float = 1.4

    # flow configuration
    beta: float = 0.0  # Hartree pressure gradient parameter: beta = 2m/(m+1)
    sweep_angle: float = 0.0  # sweep angle [deg] for FSC problems

    # transport model
    viscosity_model: str = "sutherland"
    # optional kwargs forwarded to get_transport_model() - overrides the built-in air() preset
    # e.g. {"mu_ref": 1.716e-5, "T_ref": 273.15, "S": 110.4} for a custom Sutherland model
    # None means use the default preset (air)
    viscosity_model_kwargs: dict[str, float] | None = None

    # derived quantity (computed in __post_init__, not user-settable)
    g_wall: float | None = field(default=None, init=False)  # T_wall / T_edge

    # --------------------------------------------------
    # post-init validation (runs automatically after dataclass __init__)
    #
    # include sanity checks to ensure inputs are physically meaningful
    # --------------------------------------------------
    def __post_init__(self) -> None:
        """Validate input parameters and compute derived quantities"""

        # --------------------------------------------------
        # hard error checks: raise exceptions for invalid inputs
        # --------------------------------------------------

        # normalize wall_bc to lowercase string
        object.__setattr__(self, "wall_bc", self.wall_bc.strip().lower())
        if self.wall_bc not in VALID_WALL_BC_TYPES:
            raise ValueError(
                f"wall_bc must be one of {VALID_WALL_BC_TYPES}: got {self.wall_bc!r}"
            )

        if self.mach_edge <= 0:
            raise ValueError(f"mach_edge must be positive: {self.mach_edge}")
        if self.temp_edge <= 0:
            raise ValueError(f"temp_edge must be positive: {self.temp_edge}")
        if self.prandtl <= 0:
            raise ValueError(f"prandtl must be positive: {self.prandtl}")
        if self.gamma <= 1:
            raise ValueError(f"gamma must be > 1: {self.gamma}")
        if not 0 <= self.sweep_angle < 90:
            raise ValueError(
                f"sweep_angle must be in [0, 90): {self.sweep_angle}"
            )
        # --------------------------------------------------
        # compute derived quantities
        # --------------------------------------------------

        # g_wall = T_wall / T_edge for isothermal walls
        if self.wall_bc == "isothermal":
            if self.temp_wall is None:
                raise ValueError("temp_wall must be specified for isothermal wall")
            if self.temp_wall <= 0:
                raise ValueError(f"temp_wall must be positive: {self.temp_wall}")
            object.__setattr__(self, "g_wall", self.temp_wall / self.temp_edge)
        else:
            if self.temp_wall is not None:
                raise ValueError("temp_wall should not be specified for adiabatic wall")

        # --------------------------------------------------
        # soft warnings: alert user to potentially problematic inputs but allow code to run
        # --------------------------------------------------
        if self.wall_bc == "isothermal":
            if self.g_wall < 1e-3:
                warnings.warn(
                    f"temp_wall/temp_edge = {self.g_wall:.2e} is extremely low, results may be unreliable",
                    stacklevel=2,
                )
            if self.g_wall > 1000:
                warnings.warn(
                    f"temp_wall/temp_edge = {self.g_wall:.1f} is extremely high, results may be unreliable",
                    stacklevel=2,
                )

