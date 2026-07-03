"""Initial profile builders for BVP solver initialization

Provides physics-informed starting profiles for scipy.integrate.solve_bvp.

Design
------
Two-layer design:

1. Shape functions on normalized domain in [0, 1] -> [0, 1].
   These are pure math -- no physics, no state.  The caller scales t from
   the physical eta coordinate using a boundary layer thickness estimate delta:

       t = clip(eta / delta, 0, 1)

   For eta > delta the profile is already at its edge values (fp=1, tau=1),
   giving a physically correct flat freestream tail.

2. Profile builder using the Crocco-Busemann temperature relation:

       tau(eta) = tau_w + (1 - tau_w)*fp + (tau_aw - 1)*fp*(1 - fp)

   where
       tau_aw = 1 + r*(gamma-1)/2*M^2: adiabatic wall temperature  (tau = T/T_e)
       tau_w  = Tw/Te: for isothermal; tau_aw for adiabatic
       r      = Pr^(1/3): recovery factor (Pohlhausen approximation)
       fp     = u/u_e: velocity profile from the chosen shape function

   Boundary condition check:
       tau(fp=0) = tau_w: wall BC satisfied by construction
       tau(fp=1) = 1: edge BC satisfied by construction

Available shape functions
-------------------------
    "erf"        : error function - closest match to Blasius flat-plate profile (default)
    "cubic"      : cubic - smooth, symmetric, zero-derivative at both ends
    "exponential": concave shape  - steeper near-wall gradient, fast rise

Usage
-----
Called from bvp.py when no initial_profile is available:

    from simbl.solver.profiles import build_initial_profile
    y_init = build_initial_profile(problem, solver_problem, eta_init)
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy.special import erf

if TYPE_CHECKING:
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.shooting import SolverProblem


# --------------------------------------------------
# shape functions: normalized profile on t in [0, 1]
#
# each function maps t -> fp in [0, 1] with fp(0)=0, fp(1)=1
# --------------------------------------------------

def _shape_erf(t: NDArray, c: float = 1.13, smooth_endpoints: bool = True) -> NDArray:
    """Error function shape.

    Calibrated constant c=1.13 minimises RMS error against tabulated Blasius
    values over t in [0, 1].  Increase c for a steeper near-wall gradient.

    When smooth_endpoints=True (default), t is pre-mapped through sin(pi/2*t)
    before being passed to erf.  This forces df'/dt = 0 at t=1 only (smooth
    transition into the flat freestream), while keeping a nonzero slope at
    t=0 (physically correct wall shear).  The Blasius-like interior shape
    is preserved.
    """
    # sin(pi/2*t) maps [0,1]->[0,1], zero slope at t=1, nonzero at t=0
    s = np.sin(np.pi / 2.0 * t) if smooth_endpoints else t
    # normalize so that the profile reaches 1 at the end of the domain
    return erf(c * s) / float(erf(c))


def _shape_cubic(t: NDArray) -> NDArray:
    """Cubic curve.
    """

    # cubic polynomial with zero derivative at both ends
    return 3.0 * t**2 - 2.0 * t**3


def _shape_exponential(t: NDArray, c: float = 3.0, smooth_endpoints: bool = True) -> NDArray:
    """Exponential curve.

    When smooth_endpoints=True (default), applies the same sin(pi/2*t) pre-map
    as _shape_erf: zero slope at t=1 only, nonzero at t=0.
    """
    # sin(pi/2*t): zero slope at t=1, nonzero at t=0
    s = np.sin(np.pi / 2.0 * t) if smooth_endpoints else t
    # normalize so that the profile reaches 1 at the end of the domain
    den = 1.0 - np.exp(-c)
    return (1.0 - np.exp(-c * s)) / den


# --------------------------------------------------
# registry: map string name -> shape function
# --------------------------------------------------
_SHAPE_REGISTRY: dict[str, object] = {
    "erf":         _shape_erf,
    "cubic":       _shape_cubic,
    "exponential": _shape_exponential,
}

VALID_SHAPES = set(_SHAPE_REGISTRY)


# --------------------------------------------------
# crocco_tau: Crocco-Busemann temperature approximation
# --------------------------------------------------
def crocco_tau(
    fp: NDArray[np.float64],
    tau_w: float,
    mach: float,
    gamma: float = 1.4,
    prandtl: float = 0.72,
) -> NDArray[np.float64]:
    """Crocco-Busemann temperature profile as a function of velocity

    Computes the approximate temperature ratio tau = T/Te as a function of
    the normalized velocity fp = u/ue using the Crocco-Busemann relation:

        tau(fp) = tau_w + (1 - tau_w)*fp + (tau_aw - 1)*fp*(1 - fp)

    where tau_aw = 1 + r*(gamma-1)/2*M^2 is the adiabatic wall temperature
    and r = Pr^(1/3) is the Pohlhausen recovery factor.

    Boundary conditions are satisfied by construction:
        tau(fp=0) = tau_w  (wall)
        tau(fp=1) = 1      (edge)

    Parameters
    ----------
    fp : NDArray
        Normalized velocity profile u/ue.  Values should lie in [0, 1].
    tau_w : float
        Wall temperature ratio Tw/Te.
    mach : float
        Edge Mach number.
    gamma : float, optional
        Specific heat ratio.  Default 1.4 (air).
    prandtl : float, optional
        Prandtl number.  Default 0.72 (air).

    Returns
    -------
    NDArray
        Temperature ratio T/Te, same shape as fp.
    """
    # recovery factor: Pohlhausen approximation r = Pr^(1/3)
    r = prandtl ** (1.0 / 3.0)
    # adiabatic wall temperature (non-dimensional: T/Te)
    tau_aw = 1.0 + r * (gamma - 1.0) / 2.0 * mach**2
    return tau_w + (1.0 - tau_w) * fp + (tau_aw - 1.0) * fp * (1.0 - fp)


# --------------------------------------------------
# build_initial_profile: public entry point
# --------------------------------------------------
def build_initial_profile(
    problem: SimilarityInputs,
    solver_problem: SolverProblem,
    eta: NDArray[np.float64],
    shape: str = "erf",
    delta_fraction: float = 0.8,
) -> NDArray[np.float64]:
    """Build a physics-informed BVP initial profile using the Crocco-Busemann relation

    Constructs a starting profile for scipy.integrate.solve_bvp that satisfies
    the wall and edge boundary conditions exactly and encodes the compressible
    temperature distribution via the Crocco-Busemann approximation.

    Parameters
    ----------
    problem : SimilarityInputs
        Physics specification.  Provides Mach, gamma, Pr, wall BC, and Tw/Te.
    solver_problem : SolverProblem
        Solver problem definition.  Used to determine state vector size and
        to extract the wall boundary condition initial state.
    eta : NDArray
        Similarity coordinate array for the initial mesh (output of clustered_mesh).
    shape : str, optional
        Profile shape function.  One of "erf" (default), "cubic", "exponential".
        "erf" best matches the Blasius flat-plate profile.
    delta_fraction : float, optional
        Fraction of eta_max used as the boundary layer thickness estimate delta.
        The velocity and temperature profiles reach their edge values at eta=delta
        and remain flat for eta > delta.  Default 0.8 (80% of domain).

    Returns
    -------
    NDArray
        Initial profile array, shape (n_state, len(eta)).  Ready to pass to
        scipy.integrate.solve_bvp as the ``y`` argument.
    """

    # validate shape (is it in registry?)
    if shape not in _SHAPE_REGISTRY:
        raise ValueError(
            f"shape must be one of {VALID_SHAPES}: got {shape!r}"
        )

    # assign shape function from registry based on validated shape name
    shape_func = _SHAPE_REGISTRY[shape]

    # determine state vector size from wall initial condition
    y_wall = solver_problem.build_initial_condition(solver_problem.initial_guess)

    # build_initial_condition maps shooting vars -> full y0 vector
    n_state = len(y_wall)
    n_pts   = len(eta)

    # boundary layer thickness estimate (based on delta_fraction of eta_max)
    eta_max = float(eta[-1])
    delta   = delta_fraction * eta_max

    # normalized coordinate t in [0, 1]
    t = np.clip(eta / delta, 0.0, 1.0)

    # velocity profile: fp(eta) from shape function
    fp = shape_func(t)

    # stream function f(eta): integrate fp(eta)
    # use cumulative trapezoidal integration on the (possibly non-uniform) eta grid
    f = np.zeros(n_pts)
    f[1:] = np.cumsum(0.5 * (fp[:-1] + fp[1:]) * np.diff(eta))

    # velocity gradient fpp(eta): differentiate fp(eta)
    fpp = np.zeros(n_pts)
    # centered finite difference for interior points
    fpp[1:-1] = (fp[2:] - fp[:-2]) / (eta[2:] - eta[:-2])
    # forward difference for first point
    fpp[0]    = (fp[1] - fp[0]) / (eta[1] - eta[0])
    # backward difference for last point
    fpp[-1]   = (fp[-1] - fp[-2]) / (eta[-1] - eta[-2])

    # Crocco-Busemann temperature profile

    # wall temperature: prescribed for isothermal, recovery temperature for adiabatic
    if problem.wall_bc == "isothermal":
        # g_wall = Tw/Te (computed in SimilarityInputs.__post_init__)
        tau_w = float(problem.g_wall)
    else:
        # adiabatic: wall temperature equals the recovery temperature
        r = problem.prandtl ** (1.0 / 3.0)
        tau_w = 1.0 + r * (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2

    # use the public crocco_tau function
    tau = crocco_tau(fp, tau_w, problem.mach_edge, problem.gamma, problem.prandtl)

    # temperature gradient taup(eta) = dtau/dfp * fpp
    r        = problem.prandtl ** (1.0 / 3.0)
    tau_aw   = 1.0 + r * (problem.gamma - 1.0) / 2.0 * problem.mach_edge**2
    dtau_dfp = (1.0 - tau_w) + (tau_aw - 1.0) * (1.0 - 2.0 * fp)
    taup     = dtau_dfp * fpp

    # assemble state vector
    # FS  (n_state=5): [f, fp, fpp, tau, taup]
    # FSC (n_state=7): [f, fp, fpp, tau, taup, g, gp]
    y = np.zeros((n_state, n_pts))
    y[0, :] = f
    y[1, :] = fp
    y[2, :] = fpp
    y[3, :] = tau
    y[4, :] = taup

    if n_state >= 7:
        # Falkner-Skan-Cooke: crossflow g and gp
        # use same shape as streamwise velocity (rough but reasonable starting point)
        # g(0)=0 (wall), g(inf)=gcf_edge_target (edge)
        # gcf_edge_target is 0 for aligned flow, 1 for swept -- approximate with fp
        y[5, :] = fp
        y[6, :] = fpp

    return y
