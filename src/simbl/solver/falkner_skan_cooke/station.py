"""Multi-station edge condition utilities for Falkner-Skan-Cooke similarity

The FSC similarity solver (bl_ode) always receives *local* edge conditions:
    - mach_edge  = Ma_e at the current station
    - temp_edge  = T_e at the current station

When marching along a swept surface at multiple x-stations, the edge
temperature changes with station. This module computes the local T_e from
the reference-station quantities following Liu (2021), Phys. Fluids 33,
126109, Eqs. 20-23.

Typical multi-station workflow
------------------------------
1. Pick a reference station (xi_ref) where you know all edge conditions.
2. For each downstream station xi, call compute_local_temp_edge() to get T_e.
3. Pass the result as temp_edge into SimilarityInputs for that station.

Example
-------
>>> from simbl.solver.falkner_skan_cooke.station import compute_local_temp_edge
>>> from simbl import SimilarityInputs, SolverOptions, solve_similarity
>>>
>>> # Reference conditions (e.g. attachment line)
>>> mach_ref  = 3.0
>>> temp_ref  = 280.0  # K
>>> gamma     = 1.4
>>> beta      = 0.0    # flat plate
>>> sweep_deg = 45.0
>>>
>>> # Solve at a downstream station where xi / xi_ref = 1.5
>>> xi_ratio = 1.5
>>> temp_local = compute_local_temp_edge(
...     mach_ref=mach_ref,
...     temp_ref=temp_ref,
...     xi_over_xi_ref=xi_ratio,
...     beta=beta,
...     gamma=gamma,
...     sweep_angle_deg=sweep_deg,
... )
>>>
>>> inputs = SimilarityInputs(
...     mach_edge=mach_ref,   # reference Mach (see note in function docstring)
...     temp_edge=temp_local, # local temperature computed above
...     wall_bc="adiabatic",
...     sweep_angle=sweep_deg,
...     beta=beta,
... )

Note on mach_edge
-----------------
Liu's formulation uses Ma_{e,ref} (the reference-station Mach) throughout
the ODE coefficients -- the local variation in edge velocity is absorbed into
the chi = (xi/xi_ref)^m scaling that feeds S and T_e. Because the solver
always operates at chi=1, the correct mach_edge to pass into SimilarityInputs
is still Ma_{e,ref} for all stations.

If your inviscid solution gives you an explicitly varying Ma_e(xi), use that
value directly instead of Ma_{e,ref}.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import math


# --------------------------------------------------
# public API
# --------------------------------------------------
def compute_local_temp_edge(
    mach_ref: float,
    temp_ref: float,
    xi_over_xi_ref: float,
    beta: float,
    gamma: float,
    sweep_angle_deg: float,
) -> float:
    """Compute local edge temperature at station xi from reference conditions.

    Implements Liu (2021), Phys. Fluids 33, 126109, Eqs. 20-23:

        chi = (xi / xi_ref)^m                    [Eq. 23]
        S   = 1 + (gamma-1)/2 * chi^2 * Ma_ref^2 * cos^2(Lambda)  [Eq. 22]
        S_ref = S evaluated at chi=1 (= denominator of K)
        T_e = T_{e,ref} * S_ref / S              [Eq. 20]

    Parameters
    ----------
    mach_ref : float
        Edge Mach number at the reference station (Ma_{e,ref}).
    temp_ref : float
        Edge temperature at the reference station [K].
    xi_over_xi_ref : float
        Station ratio xi / xi_ref. Use 1.0 to recover T_{e,ref}.
    beta : float
        Hartree pressure gradient parameter beta = 2m/(m+1).
        beta=0 is a flat plate; beta=1 is stagnation.
    gamma : float
        Specific heat ratio (1.4 for air).
    sweep_angle_deg : float
        Sweep angle [degrees].

    Returns
    -------
    float
        Local edge temperature T_e [K] at station xi.
        Pass this as temp_edge into SimilarityInputs.

    Raises
    ------
    ValueError
        If mach_ref, temp_ref, or xi_over_xi_ref are non-positive.
    """
    # --------------------------------------------------
    # validate inputs
    # --------------------------------------------------
    if mach_ref <= 0:
        raise ValueError(f"mach_ref must be positive: {mach_ref}")
    if temp_ref <= 0:
        raise ValueError(f"temp_ref must be positive: {temp_ref}")
    if xi_over_xi_ref <= 0:
        raise ValueError(f"xi_over_xi_ref must be positive: {xi_over_xi_ref}")

    # --------------------------------------------------
    # compute Falkner-Skan exponent m from Hartree beta
    # beta = 2m/(m+1)  =>  m = beta / (2 - beta)
    # stagnation point: beta -> 2 gives m -> infinity (handled as large value)
    # --------------------------------------------------
    if abs(beta - 2.0) > 1e-10:
        m = beta / (2.0 - beta)
    else:
        # stagnation point limit: m -> infinity, chi^(2m) -> large; T_e -> T_ref / S
        m = 1e10

    # --------------------------------------------------
    # compute chi = (xi / xi_ref)^m  [Liu Eq. 23]
    # for flat plate (m=0): chi = 1 regardless of xi_over_xi_ref
    # --------------------------------------------------
    chi = xi_over_xi_ref**m

    # --------------------------------------------------
    # build the compressibility group for cos^2(Lambda)
    # --------------------------------------------------
    sweep_rad = math.radians(sweep_angle_deg)
    cos2_lambda = math.cos(sweep_rad) ** 2

    # --------------------------------------------------
    # S at reference station (chi=1): S_ref = 1 + (gamma-1)/2 * Ma_ref^2 * cos^2(Lambda)
    # S at local station:             S     = 1 + (gamma-1)/2 * chi^2 * Ma_ref^2 * cos^2(Lambda)
    # [Liu Eq. 22]
    # --------------------------------------------------
    compressibility_group = (gamma - 1.0) / 2.0 * mach_ref**2 * cos2_lambda
    S_ref = 1.0 + compressibility_group
    S_local = 1.0 + chi**2 * compressibility_group

    # --------------------------------------------------
    # local edge temperature from Liu Eq. 20: T_e = T_{e,ref} * S_ref / S_local
    # --------------------------------------------------
    temp_local = temp_ref * S_ref / S_local

    return temp_local
