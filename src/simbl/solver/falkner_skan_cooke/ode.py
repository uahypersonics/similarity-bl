"""ODE system for compressible Falkner-Skan-Cooke boundary layer equations

Three-dimensional extension of the similarity boundary layer equations
with crossflow, following Liu (2021), Phys. Fluids 33, 126109.

Implementation of Liu Eqs. 16-18 using state vector:
    y = [f, f', f'', tau, tau', g, g']

Where:
    f    = streamwise stream function
    f'   = u/u_e (streamwise velocity ratio)
    f''  = streamwise velocity gradient
    tau  = T/T_e (normalized temperature)
    tau' = temperature gradient
    g    = w/w_e (crossflow velocity ratio)
    g'   = crossflow velocity gradient
    eta   = similarity coordinate

Edge conditions:
    M   = edge Mach number (streamwise)
    Lambda = sweep angle
    Pr  = Prandtl number
    gamma = specific heat ratio
    beta = 2m/(m+1) (pressure gradient parameter)
    mu  = viscosity ratio mu/mu_e

Note: Variable naming follows Liu (2021) convention to avoid confusion:
    - g is crossflow velocity (not to be confused with temperature)
    - tau is temperature (not to be confused with shear stress)

    - C = rho*mu / (rho_e*mu_e) is the Chapman-Rubesin factor
    - perfect gas: p = rho R T and p_e = rho_e R T_e
    - dp/dy = 0 implies p = p_e everywhere in the boundary layer
    - Therefore: rho/rho_e = 1/tau for perfect gas at constant pressure
    - C = mu/mu_e * rho/rho_e = mu/tau  (where mu = mu_local/mu_e)
    - mu'/mu = (T_e/mu_local)(dmu/dT) * tau'  (chain rule via T = T_e*tau)
    - precomputed as visc_term = dmu_dT/mu = T_e/mu_local * dmu/dT
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
    S: float,
    K: float,
) -> NDArray[np.float64]:
    """Compressible Falkner-Skan-Cooke ODE system (7 equations)

    Implements Liu (2021), Phys. Fluids 33, 126109, Eqs. 16-18.

    S and K are loop-invariant (depend only on problem constants, not on eta
    or the state vector). They are precomputed once in build_solver_problem
    and passed in here so this function does not recompute them on every call.

    Parameters
    ----------
    eta : float
        Similarity coordinate (independent variable).
    y : NDArray[np.float64]
        State vector [f, f', f'', tau, tau', g, g'].
    problem : SimilarityInputs
        Physics specification (Mach, Pr, gamma, beta, temp_edge).
    visc_model : TransportModel
        Viscosity model instance (from flow_state).
    S : float
        Stagnation enthalpy parameter (Liu Eq. 22, at chi=1 for local station).
        S = 1 + (gamma-1)/2 * Ma_e^2 * cos^2(Lambda)
    K : float
        Crossflow recovery factor (Liu Eq. 21).
        K = [1 + (gamma-1)/2 * Ma_e^2] / S

    Returns
    -------
    NDArray[np.float64]
        Derivatives [f', f'', f''', g', g'', tau', tau''].
    """

    # unpack state vector
    f, f_p, f_pp, tau, tau_p, g, g_p = y

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
    # precompute T_e/mu_local * dmu/dT (reused in all three equations)
    visc_term = problem.temp_edge * dmudt_val / mu_val

    # Initialize derivatives array
    dy = np.zeros(7)

    # --------------------------------------------------
    # x-momentum (Liu Eq. 16)
    # (C f'')' + f f'' + beta_H*(tau - f'^2) = 0
    #
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
    # z-momentum / crossflow (Liu Eq. 17)
    # (C g')' + f g' = 0
    #
    # g'' = (tau'/tau)*g' - (T_e/mu)*(dmu/dT)*tau'*g' - f*g'/C
    # --------------------------------------------------
    dy[5] = g_p
    dy[6] = (
        tau_p / tau * g_p
        - visc_term * tau_p * g_p
        - f * g_p / C
    )

    # --------------------------------------------------
    # energy (Liu Eq. 18)
    # tau'' is computed last — uses f''' (dy[2]) and g'' (dy[6]) already computed above
    #
    # tau'' = tau'^2/tau - (T_e/mu)*(dmu/dT)*tau'^2
    #       - Pr*f*tau'/C
    #       - 2Pr(S-1)*[(T_e/mu*dmu/dT - 1/tau)*tau'*f'*f'' + f''^2 + f'*f''']
    #       - 2Pr(K-1)S*[(T_e/mu*dmu/dT - 1/tau)*tau'*g*g' + g'^2 + g*g'']
    #       - 2Pr*f/C*[(S-1)*f'*f'' + (K-1)*S*g*g']
    # --------------------------------------------------
    f_ppp = dy[2]
    g_pp = dy[6]

    dy[3] = tau_p
    dy[4] = (
        tau_p**2 / tau
        - visc_term * tau_p**2
        - problem.prandtl * f * tau_p / C
        - 2.0 * problem.prandtl * (S - 1.0) * (
            (visc_term - 1.0 / tau) * tau_p * f_p * f_pp
            + f_pp**2
            + f_p * f_ppp
        )
        - 2.0 * problem.prandtl * (K - 1.0) * S * (
            (visc_term - 1.0 / tau) * tau_p * g * g_p
            + g_p**2
            + g * g_pp
        )
        - 2.0 * problem.prandtl * f / C * (
            (S - 1.0) * f_p * f_pp + (K - 1.0) * S * g * g_p
        )
    )

    return dy
