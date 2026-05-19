"""ODE system for compressible Falkner-Skan-Cooke boundary layer equations

Three-dimensional extension of the similarity boundary layer equations
with crossflow, following Liu (2021), Phys. Fluids 33, 126109.

The system of 7 equations in transformed variables:
    f'   = fp
    fp'  = fpp
    fpp' = (g/mu)[-f*fpp + fpp*gp*(mu/g - dmu/dT)/g + beta*fp^2 - beta*g]
    w'   = wp
    wp'  = (g/mu)[-f*wp + wp*gp*(mu/g - dmu/dT)/g]
    g'   = gp
    gp'  = (g/mu)[gp^2*(mu/g - dmu/dT)/g - Pr*f*gp
                   - Pr*mu*(gamma-1)*M^2*(fpp^2 + tan^2(Lambda)*wp^2)/g]

Where:
    f   = streamwise stream function
    w   = crossflow velocity ratio w/w_e
    g   = T/T_e (normalized temperature)
    eta = similarity coordinate
    M   = edge Mach number (streamwise)
    Lambda = sweep angle
    Pr  = Prandtl number
    gamma = specific heat ratio
    beta = pressure gradient parameter
    mu  = viscosity ratio mu/mu_e
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from flow_state.transport import TransportModel

    from simbl.solver.inputs import SimilarityInputs


# --------------------------------------------------
# ODE system for Falkner-Skan-Cooke similarity equations
# --------------------------------------------------
def bl_ode(
    eta: float,
    y: NDArray[np.float64],
    problem: SimilarityInputs,
    visc_model: TransportModel,
    tan2_sweep: float,
) -> NDArray[np.float64]:
    """Compressible Falkner-Skan-Cooke ODE system (7 equations)

    Parameters
    ----------
    eta : float
        Similarity coordinate (independent variable).
    y : NDArray[np.float64]
        State vector [f, f', f'', w, w', g, g'].
    problem : SimilarityInputs
        Physics specification (Mach, Pr, gamma, beta, temp_edge).
    visc_model : TransportModel
        Viscosity model instance (from flow_state).
    tan2_sweep : float
        Crossflow dissipation parameter tan^2(Lambda).

    Returns
    -------
    NDArray[np.float64]
        Derivatives [f', f'', f''', w', w'', g', g''].
    """
    # unpack state vector
    f, fp, fpp, _w, wp, g, gp = y

    # compute dimensional temperature
    # clamp to a small positive value to survive Newton-Raphson trial steps;
    # intermediate shooting variables can produce unphysical g during iteration
    # and the viscosity model raises ValueError on non-positive temperatures
    temp = max(g * problem.temp_edge, 1.0)

    # viscosity from flow_state transport model
    mu_val = visc_model.mu(temp)
    mu_edge = visc_model.mu(problem.temp_edge)
    dmudt_val = visc_model.dmudt(temp)

    # normalized viscosity: mu = mu(T) / mu(T_e)
    mu = mu_val / mu_edge
    # dmu_dT = T_e / mu_e * (d mu/d T)
    dmu_dT = problem.temp_edge / mu_edge * dmudt_val

    # coefficient from expanding d/d(eta) of the Chapman-Rubesin parameter C = mu/g
    # C' = d/deta (mu/g) = (mu/g - dmu/dT) / g * g'
    # C'/g' = (mu/g - dmu_dT) / g
    # this term appears repeatedly in the equations, so we compute it once here
    cr_deriv = (mu / g - dmu_dT) / g

    # initialize derivatives array
    dy = np.zeros(7)

    # --------------------------------------------------
    # streamwise momentum (identical to 2D)
    # --------------------------------------------------
    dy[0] = fp
    dy[1] = fpp
    dy[2] = (g / mu) * (
        -f * fpp
        + fpp * gp * cr_deriv
        + problem.beta * fp**2
        - problem.beta * g
    )

    # --------------------------------------------------
    # crossflow: (N w')' + f w' = 0
    # --------------------------------------------------
    dy[3] = wp
    dy[4] = (g / mu) * (
        -f * wp
        + wp * gp * cr_deriv
    )

    # --------------------------------------------------
    # energy (with crossflow dissipation)
    # --------------------------------------------------
    dy[5] = gp

    # clip gp before squaring to prevent overflow during Newton-Raphson trial
    # steps; unphysical shooting guesses can produce very large gp values that
    # cause gp**2 * cr_deriv to overflow to inf, stalling the ODE integrator
    gp_safe = np.clip(gp, -500.0, 500.0)

    dy[6] = (g / mu) * (
        gp_safe**2 * cr_deriv
        - problem.prandtl * f * gp
        - problem.prandtl * mu * (problem.gamma - 1) * problem.mach_edge**2 * (fpp**2 + tan2_sweep * wp**2) / g
    )

    # replace any remaining inf/nan with a large finite value so that the ODE
    # integrator receives a valid derivative and can step toward convergence
    # rather than signalling failure and triggering the slow BDF fallback path
    dy = np.where(np.isfinite(dy), dy, np.sign(dy) * 1.0e30)

    return dy
