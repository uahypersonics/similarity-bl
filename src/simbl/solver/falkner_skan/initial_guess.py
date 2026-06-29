"""Initial guess handling for the Falkner-Skan (2D) model

Owns:
- FS table configs (adiabatic + isothermal)
- build_y0: shooting variables --> 5-element ODE initial condition
- get_initial_guess: public entry point for predictions
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from simbl.solver.lookup_table import LookupTable

# --------------------------------------------------
# table configurations (Falkner-Skan, 2 shooting variables)
#
# g_wall is a value (unknown) for adiabatic, but a key (prescribed) for isothermal
# adiabatic:  key = (mach, beta),          values = (fpp_wall, g_wall)
# isothermal: key = (mach, beta, g_wall),  values = (fpp_wall, gp_wall)
# --------------------------------------------------
_TABLE_CONFIGS: dict[str, dict] = {
    "adiabatic": {
        "key_fields": ["mach", "beta", "temp_edge"],
        "value_fields": ["fpp_wall", "g_wall"],
        "default_values": [0.55, 3.0],
        "fname": "lookup_fs_adiabatic.json",
    },
    "isothermal": {
        "key_fields": ["mach", "beta", "g_wall"],
        "value_fields": ["fpp_wall", "gp_wall"],
        "default_values": [0.50, 1.0],
        "fname": "lookup_fs_isothermal.json",
    },
}


# --------------------------------------------------
# build_y0: shooting variables --> ODE initial condition [f, f', f'', g, g']
# --------------------------------------------------
def build_y0(
    wall_bc: str,
    g_wall: float | None,
    shooting_vars: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Build 5-element initial condition vector for FS ODE integration

    Parameters
    ----------
    wall_bc : str
        Wall boundary condition type: "isothermal" or "adiabatic"
    g_wall : float or None
        Normalized wall temperature T_w / T_e (isothermal only)
    shooting_vars : NDArray[np.float64]
        Shooting variables [f''(0), g_var]

    Returns
    -------
    NDArray[np.float64]
        Initial conditions [f(0), f'(0), f''(0), g(0), g'(0)]
    """
    fpp_0 = shooting_vars[0]
    g_var = shooting_vars[1]

    if wall_bc == "isothermal":
        # wall temperature prescribed; g'(0) is the shooting variable
        g_0 = g_wall
        gp_0 = g_var
    else:
        # adiabatic: zero heat flux g'(0) = 0; g(0) is the shooting variable
        g_0 = g_var
        gp_0 = 0.0

    return np.array([0.0, 0.0, fpp_0, g_0, gp_0])


# --------------------------------------------------
# get_initial_guess: public entry point for predictions
# --------------------------------------------------
def get_initial_guess(
    wall_bc: str,
    key_values: dict[str, float],
) -> NDArray[np.float64]:
    """Get initial guess for FS shooting method

    Parameters
    ----------
    wall_bc : str
        Wall boundary condition type: "isothermal" or "adiabatic"
    key_values : dict[str, float]
        Key field values for lookup (must contain all key_fields for this wall BC)

    Returns
    -------
    NDArray[np.float64]
        Predicted shooting variable values.
    """

    # initialize LookupTable (defined in lookup_table.py) for this wall BC type
    # ** expands the config dict keys into keyword arguments for the LookupTable constructor
    lookup = LookupTable(**_TABLE_CONFIGS[wall_bc])

    # predict shooting variables for this query point
    return lookup.predict(key_values)
