"""Initial guess handling for the Falkner-Skan-Cooke (3D) model

Owns:
- FSC table configs (adiabatic + isothermal)
- build_y0: shooting variables --> 7-element ODE initial condition
- get_initial_guess: public entry point for predictions
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from simbl.solver.lookup_table import LookupTable

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.shooting import ShootingResult


# --------------------------------------------------
# table configurations (Falkner-Skan-Cooke, 3 shooting variables)
#
# g_wall (= tau_wall = T_wall/T_edge) is a value (unknown) for adiabatic, but a key (prescribed) for isothermal
# Note: Using 'g_wall' name for compatibility with SimilarityInputs.g_wall attribute
# adiabatic:  key = (mach, beta, sweep_angle),          values = (fpp_wall, gcfp_wall, g_wall)
# isothermal: key = (mach, beta, sweep_angle, g_wall),  values = (fpp_wall, gcfp_wall, gp_wall)
# --------------------------------------------------
_TABLE_CONFIGS: dict[str, dict] = {
    "adiabatic": {
        "key_fields": ["mach", "beta", "sweep_angle"],
        "value_fields": ["fpp_wall", "gcfp_wall", "g_wall"],
        "default_values": [0.55, 0.55, 3.0],
        "fname": "lookup_fsc_adiabatic.json",
    },
    "isothermal": {
        "key_fields": ["mach", "beta", "sweep_angle", "g_wall"],
        "value_fields": ["fpp_wall", "gcfp_wall", "gp_wall"],
        "default_values": [0.50, 0.50, 1.0],
        "fname": "lookup_fsc_isothermal.json",
    },
}


# --------------------------------------------------
# build_y0: shooting variables --> ODE initial condition [f, f', f'', tau, tau', g, g']
# --------------------------------------------------
def build_y0(
    wall_bc: str,
    tau_wall: float | None,
    shooting_vars: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Build 7-element initial condition vector for FSC ODE integration

    Parameters
    ----------
    wall_bc : str
        Wall boundary condition type: "isothermal" or "adiabatic".
    tau_wall : float or None
        Normalized wall temperature T_w / T_e (isothermal only).
    shooting_vars : NDArray[np.float64]
        Shooting variables [f''(0), g_cf'(0), tau_var].

    Returns
    -------
    NDArray[np.float64]
        Initial conditions [f(0), f'(0), f''(0), tau(0), tau'(0), g(0), g'(0)].
    """
    fpp_0 = shooting_vars[0]
    gcfp_0 = shooting_vars[1]
    tau_var = shooting_vars[2]

    if wall_bc == "isothermal":
        # Wall temperature prescribed; tau'(0) is the shooting variable
        tau_0 = tau_wall
        taup_0 = tau_var
    else:
        # Adiabatic: zero heat flux tau'(0) = 0; tau(0) is the shooting variable
        tau_0 = tau_var
        taup_0 = 0.0

    return np.array([0.0, 0.0, fpp_0, tau_0, taup_0, 0.0, gcfp_0])


# --------------------------------------------------
# get_initial_guess: public entry point for predictions
# --------------------------------------------------
def get_initial_guess(
    wall_bc: str,
    key_values: dict[str, float],
) -> NDArray[np.float64]:
    """Get initial guess for FSC shooting method

    Parameters
    ----------
    wall_bc : str
        Wall boundary condition type: "isothermal" or "adiabatic".
    key_values : dict[str, float]
        Key field values for lookup (must contain all key_fields for this wall BC).

    Returns
    -------
    NDArray[np.float64]
        Predicted shooting variable values.
    """

    # initialize LookupTable (defined in lookup_table.py) for this wall BC type
    # ** expands the config dict keys into keyword arguments for the LookupTable constructor
    lookup = LookupTable(**_TABLE_CONFIGS[wall_bc])

    # for sweep_angle = 0, bootstrap from the Falkner-Skan lookup table
    # (FSC at lambda=0 reduces to FS in the streamwise/energy fields, and
    # the FSC adiabatic table has no entries at sweep_angle=0 yet, so
    # extrapolation in the sweep dimension is unreliable). The crossflow
    # shooting variable gcf_p(0) gets a sensible default since crossflow decouples.
    if key_values.get("sweep_angle", None) == 0.0:
        from simbl.solver.falkner_skan.initial_guess import (
            _TABLE_CONFIGS as _FS_TABLE_CONFIGS,
        )
        fs_lookup = LookupTable(**_FS_TABLE_CONFIGS[wall_bc])
        # FS keys are a subset of FSC keys (no sweep_angle)
        fs_key_values = {k: v for k, v in key_values.items() if k != "sweep_angle"}
        # FS predict returns [fpp_wall, g_wall] (adiabatic) or [fpp_wall, gp_wall] (isothermal)
        # where g_wall = tau_wall = T_wall/T_edge
        fs_guess = fs_lookup.predict(fs_key_values)
        # build the 3-element FSC guess: [fpp_wall, gcfp_wall, g_var/gp_var]
        # default gcfp_wall = 0.5 (consistent with low-sweep typical values)
        return np.array([fs_guess[0], 0.5, fs_guess[1]])

    # predict shooting variables for this query point
    return lookup.predict(key_values)
