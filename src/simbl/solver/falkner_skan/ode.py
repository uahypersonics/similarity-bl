"""ODE system for compressible Falkner-Skan similarity equations

The system of 5 equations in transformed variables:

    state vector:
    y0 = f
    y1 = f'
    y2 = f''
    y3 = tau
    y4 = tau'

    y0' = f' = y1
    y1' = f'' = y2
    y2' = f''' =
    y3' = tau' = y4
    y4' = tau'' =

Where:
    f    = stream function
    tau  = T/T_e (normalized temperature)
    eta  = similarity coordinate
    M    = edge Mach number
    Pr   = Prandtl number
    gamma = specific heat ratio
    beta = pressure gradient parameter
    mu   = viscosity ratio (mu/mu_e)
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
        State vector [f, f', f'', tau, tau'].
    problem : SimilarityInputs
        Physics specification (Mach, Pr, gamma, beta, temp_edge).
    visc_model : TransportModel
        Viscosity model instance (from flow_state).

    Returns
    -------
    NDArray[np.float64]
        Derivatives [f', f'', f''', tau', tau''].
    """

    # unpack state vector
    f, f_p, f_pp, tau, tau_p = y

    # compute dimensional temperature (needed for viscosity model)
    # clamp temperature to a small positive value to survive Newton-Raphson trial steps
    temp = max(tau * problem.temp_edge, 1.0)

    # viscosity from flow_state transport model
    mu_val = visc_model.mu(temp)
    mu_edge = visc_model.mu(problem.temp_edge)
    dmudt_val = visc_model.dmudt(temp)

    # normalized viscosity: mu = mu(T) / mu(T_e)
    mu = mu_val / mu_edge

    # Chapman-Rubesin factor: C = rho*mu/(rho_e*mu_e) = mu/tau for perfect gas
    C = mu / tau

    # mu'/mu = (T_e/mu_local)*(dmu/dT)*tau'  (chain rule, T = T_e*tau)
    # precompute T_e/mu_local * dmu/dT (reused in both momentum and energy)
    visc_term = problem.temp_edge * dmudt_val / mu_val

    # initialize derivatives array
    dy = np.zeros(5)

    # --------------------------------------------------
    # streamwise momentum
    # f''' = (tau'/tau)*f'' - (T_e/mu)*(dmu/dT)*tau'*f'' - f*f''/C - beta_H*(tau-f'^2)/C
    # --------------------------------------------------
    dy[0] = f_p
    dy[1] = f_pp
    dy[2] = (
        tau_p / tau * f_pp
        - visc_term * tau_p * f_pp
        - (f * f_pp + problem.beta * (tau - f_p**2)) / C
    )

    # --------------------------------------------------
    # energy
    # tau'' = tau'^2/tau - (T_e/mu)*(dmu/dT)*tau'^2
    #       - Pr*f*tau'/C
    #       - Pr*(gamma-1)*M_e^2*[f''^2 - beta_H*tau*f'/C]
    # --------------------------------------------------
    dy[3] = tau_p
    dy[4] = (
        tau_p**2 / tau
        - visc_term * tau_p**2
        - problem.prandtl * f * tau_p / C
        - problem.prandtl * (problem.gamma - 1) * problem.mach_edge**2
        * (f_pp**2 - problem.beta * tau * f_p / C)
    )

    return dy
