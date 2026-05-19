"""ODE system for compressible Falkner-Skan similarity equations

The system of 5 equations in transformed variables:
    f'   = fp
    fp'  = fpp
    fpp' = (g/mu)[-f*fpp + fpp*gp*(mu/g - dmu/dT)/g + beta*fp^2 - beta*g]
    g'   = gp
    gp'  = (g/mu)[gp^2*(mu/g - dmu/dT)/g - Pr*f*gp - Pr*mu*(gamma-1)*M^2*fpp^2/g]

Where:
    f   = stream function
    g   = T/T_e (normalized temperature)
    eta = similarity coordinate
    M   = edge Mach number
    Pr  = Prandtl number
    gamma = specific heat ratio
    beta = pressure gradient parameter
    mu  = viscosity ratio (mu/mu_e)
"""""

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
# ODE system for Falkner-Skan similarity equations
# --------------------------------------------------
def bl_ode(
    eta: float,
    y: NDArray[np.float64],
    problem: SimilarityInputs,
    visc_model: TransportModel,
) -> NDArray[np.float64]:
    """Compressible boundary layer ODE system (5 equations)

    Parameters
    ----------
    eta : float
        Similarity coordinate (independent variable).
    y : NDArray[np.float64]
        State vector [f, f', f'', g, g'].
    problem : SimilarityInputs
        Physics specification (Mach, Pr, gamma, beta, temp_edge).
    visc_model : TransportModel
        Viscosity model instance (from flow_state).

    Returns
    -------
    NDArray[np.float64]
        Derivatives [f', f'', f''', g', g''].
    """

    # unpack state vector
    f, fp, fpp, g, gp = y

    # compute dimensional temperature
    temp = g * problem.temp_edge

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
    dy = np.zeros(5)

    # --------------------------------------------------
    # streamwise momentum
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
    # energy
    # --------------------------------------------------
    dy[3] = gp
    dy[4] = (g / mu) * (
        gp**2 * cr_deriv
        - problem.prandtl * f * gp
        - problem.prandtl * mu * (problem.gamma - 1) * problem.mach_edge**2 * fpp**2 / g
    )

    return dy
