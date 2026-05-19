"""Wall boundary condition types

Defines the valid wall boundary condition types used internally
by the solver. Users specify wall_bc as a plain string
("isothermal" or "adiabatic") on SimilarityInputs.

Isothermal wall
    Known:   g(0) = T_w / T_e
    Unknown: g'(0)  -> shooting variable

Adiabatic wall
    Known:   g'(0) = 0  (zero heat flux)
    Unknown: g(0)   -> shooting variable
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from enum import StrEnum

# valid wall BC type strings (for validation)
VALID_WALL_BC_TYPES = {"isothermal", "adiabatic"}


# --------------------------------------------------
# wall boundary condition type enum (internal use)
# --------------------------------------------------
class WallBCType(StrEnum):
    """Wall boundary condition types (internal)"""

    ISOTHERMAL = "isothermal"
    ADIABATIC = "adiabatic"
