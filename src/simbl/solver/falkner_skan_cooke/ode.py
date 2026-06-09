"""ODE system for compressible Falkner-Skan-Cooke boundary layer equations

Three-dimensional extension of the similarity boundary layer equations
with crossflow, following Liu (2021), Phys. Fluids 33, 126109.

Implementation of Liu Eqs. 16-18 using state vector:
    y = [f, fp, fpp, g_cf, gcf_p, tau, tau_p]

Where:
    f     = streamwise stream function
    fp    = f' = u/u_e (streamwise velocity ratio)
    fpp   = f'' (streamwise velocity gradient)
    g_cf  = w/w_e (crossflow velocity ratio)
    gcf_p = g_cf' (crossflow velocity gradient)
    tau   = T/T_e (normalized temperature)
    tau_p = tau' (temperature gradient)
    eta   = similarity coordinate

Edge conditions:
    M   = edge Mach number (streamwise)
    Lambda = sweep angle
    Pr  = Prandtl number
    gamma = specific heat ratio
    beta = 2m/(m+1) (pressure gradient parameter)
    mu  = viscosity ratio mu/mu_e

Note: Variable naming follows Liu (2021) convention to avoid confusion:
    - g_cf is crossflow velocity (not to be confused with temperature)
    - tau is temperature (not to be confused with shear stress)

    - N = rho*mu / (rho_e*mu_e) is the Chapman-Rubesin parameter
    - perfect gas: p = rho R T and p_e = rho_e R T_e
    - dp/dy = 0 implies p = p_e everywhere in the boundary layer
    - Therefore: p/p_e = 1 = (rho R T) / (rho_e R T_e) = (rho/rho_e) * (T/T_e) = (rho/rho_e) * tau
    - rho/rho_e = 1/tau for perfect gas at constant pressure
    - N = mu/mu_e * rho/rho_e = mu_ratio / tau
    - N' = d(N)/deta = d/deta(mu_ratio/tau) = (mu_ratio' * tau - mu_ratio * tau') / tau^2
    - N' = mu_ratio'/tau - mu_ratio * tau' / tau^2
    - mu_ratio' = d(mu_ratio)/deta = d/deta(mu/mu_e) = (dmu/dT) * (dT/deta) / mu_e
    - tau = T/T_e => dT/deta = T_e * tau' => mu_ratio' = (dmu/dT) * (T_e * tau') / mu_e
    - N' = (dmu/dT) * (T_e * tau') / mu_e * (1/tau) - mu_ratio * tau' / tau^2
    - N' = tau'/tau * [ (dmu/dT) * (T_e / mu_e) - mu_ratio / tau ]
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
) -> NDArray[np.float64]:
    """Compressible Falkner-Skan-Cooke ODE system (7 equations)

    Implements Liu (2021), Phys. Fluids 33, 126109, Eqs. 16-18.

    Parameters
    ----------
    eta : float
        Similarity coordinate (independent variable).
    y : NDArray[np.float64]
        State vector [f, fp, fpp, g_cf, gcf_p, tau, tau_p].
    problem : SimilarityInputs
        Physics specification (Mach, Pr, gamma, beta, temp_edge).
    visc_model : TransportModel
        Viscosity model instance (from flow_state).

    Returns
    -------
    NDArray[np.float64]
        Derivatives [f', f'', f''', g_cf', g_cf'', tau', tau''].
    """
    # --------------------------------------------------
    # Unpack state vector
    # --------------------------------------------------
    f, fp, fpp, g_cf, gcf_p, tau, tau_p = y

    # --------------------------------------------------
    # Compute viscosity properties
    # --------------------------------------------------
    # Clamp temperature to a small positive value to survive Newton-Raphson trial steps;
    # intermediate shooting variables can produce unphysical tau during iteration
    # and the viscosity model raises ValueError on non-positive temperatures
    temp = max(tau * problem.temp_edge, 1.0)

    # Viscosity from flow_state transport model
    mu_val = visc_model.mu(temp)
    mu_edge = visc_model.mu(problem.temp_edge)
    dmudt_val = visc_model.dmudt(temp)

    # Normalized viscosity: mu_ratio = mu(T) / mu(T_e)
    mu_ratio = mu_val / mu_edge
    # dmu_dT = T_e / mu_e * (d mu/d T)
    dmu_dT = problem.temp_edge / mu_edge * dmudt_val

    # --------------------------------------------------
    # Chapman-Rubesin parameter: N = rho*mu / (rho_e*mu_e)
    # For perfect gas at constant pressure: rho/rho_e = 1/tau
    # Therefore: N = mu_ratio / tau
    # --------------------------------------------------
    N = mu_ratio / tau

    # --------------------------------------------------
    # Derivative of the Chapman-Rubesin parameter: N'
    # N = mu_ratio / tau depends on eta only through tau, so by the chain /
    # quotient rule: N' = tau' (dmu_dT - mu_ratio/tau) / tau
    # Computed here and reused in every equation below (both momentum
    # equations and the energy equation) so the algebra lives in one place.
    # --------------------------------------------------
    N_prime = tau_p * (dmu_dT - mu_ratio / tau) / tau

    # Initialize derivatives array
    dy = np.zeros(7)

    # --------------------------------------------------
    # Compute Liu parameters S, K, and related quantities
    # --------------------------------------------------
    # mach_e_ref = edge Mach number (streamwise component)
    mach_e_ref = problem.mach_edge

    # Sweep angle Lambda
    sweep_rad = np.radians(problem.sweep_angle)
    cos_lambda = np.cos(sweep_rad)
    cos2_lambda = cos_lambda**2

    # Hartree pressure gradient parameter beta = 2m/(m+1)
    # Solve for Falkner-Skan exponent m: m = beta/(2-beta) for beta ≠ 2
    hartree_beta = problem.beta
    if abs(hartree_beta - 2.0) > 1e-10:
        m = hartree_beta / (2.0 - hartree_beta)
    else:
        m = 1e10  # large value for beta = 2 (stagnation point limit)

    # Liu Eq. 23: v = (xi/xi_ref)^m
    # For flat plate (m=0), v = 1 regardless of xi/xi_ref
    # For pressure gradient flows, v varies with streamwise position
    xi_over_xi_ref = problem.xi_over_xi_ref
    v = xi_over_xi_ref**m

    # Liu Eq. 22: S = 1 + (gamma-1)/2 * v^2 * mach_e_ref^2 * cos^2(Lambda)
    # S is the stagnation enthalpy parameter (NOT absorbed into beta)
    S = 1.0 + (problem.gamma - 1.0) / 2.0 * v**2 * mach_e_ref**2 * cos2_lambda

    # Liu Eq. 21: K = [1 + (gamma-1)/2 * mach_e_ref^2] /
    #                  [1 + (gamma-1)/2 * mach_e_ref^2 * cos^2(Lambda)]
    # K is the crossflow recovery factor
    numerator = 1.0 + (problem.gamma - 1.0) / 2.0 * mach_e_ref**2
    denominator = 1.0 + (problem.gamma - 1.0) / 2.0 * mach_e_ref**2 * cos2_lambda
    K = numerator / denominator

    # --------------------------------------------------
    # Liu Eq. 16: Streamwise momentum
    # (N f'')' + f f'' = (beta / S) * (fp^2 - tau)
    #
    # where:
    #   beta = 2m/(m+1) is the Hartree pressure gradient parameter
    #   S is the stagnation enthalpy parameter (Liu Eq. 22)
    #   Liu Eq. 16 RHS is 2m/[(m+1) S] = beta / S (S appears in the denominator)
    #
    # Expanding (N f'')' = N f''' + N' f''
    # where N' = tau' (dmu_dT - mu_ratio/tau) / tau
    #
    # Therefore:
    # N f''' = -f f'' + (beta / S) * (fp^2 - tau) - N' f''
    # f''' = (tau/mu_ratio) [-f f'' + (beta / S) * (fp^2 - tau)
    #                        + f'' tau' (mu_ratio/tau - dmu_dT) / tau]
    # --------------------------------------------------
    dy[0] = fp
    dy[1] = fpp

    # Using N' (computed once above): the inline group
    # fpp * tau' (mu_ratio/tau - dmu_dT)/tau equals -N' fpp,
    # and tau/mu_ratio = 1/N. So Liu Eq. 16 becomes:
    # f''' = (1/N) [-f f'' + (beta/S)(fp^2 - tau) - N' f'']

    dy[2] = (1.0 / N) * (
        -f * fpp
        + hartree_beta / S * (fp**2 - tau)
        - N_prime * fpp
    )

    # --------------------------------------------------
    # Liu Eq. 17: Crossflow momentum
    # (N g_cf')' + f g_cf' = 0
    #
    # Expanding (N g_cf')' = N g_cf'' + N' g_cf'
    # N g_cf'' = (N g_cf')' - N' g_cf'
    #          = -f g_cf' - g_cf' tau' (dmu_dT - mu_ratio/tau)/tau
    # g_cf'' = (1/N) [-f g_cf' - g_cf' tau' (dmu_dT - mu_ratio/tau)/tau]
    # --------------------------------------------------
    dy[3] = gcf_p

    # Same structure as Liu Eq. 16 but with no pressure-gradient term:
    # g_cf'' = (1/N) [-f g_cf' - N' g_cf']
    dy[4] = (1.0 / N) * (
        -f * gcf_p
        - N_prime * gcf_p
    )

    # --------------------------------------------------
    # Liu Eq. 18: Energy equation (EXACT FORM)
    # (N/Pr tau')' + (S - 1) [N (fp^2)']' + (K - 1) [S N (g_cf^2)']'
    # + f [tau' + (S - 1)(fp^2)' + (K - 1) S (g_cf^2)'] = 0
    #
    # where:
    #   N = mu_ratio / tau (Chapman-Rubesin parameter)
    #   K = [1 + (gamma-1)/2 * mach_e_ref^2] / [1 + (gamma-1)/2 * mach_e_ref^2 * cos^2(Lambda)]  [Liu Eq. 21]
    #   S = 1 + (gamma-1)/2 * v^2 * mach_e_ref^2 * cos^2(Lambda)  [Liu Eq. 22]
    #   v = (xi/xi_ref)^m, for ZPG (m=0): v = 1  [Liu Eq. 23]
    #
    # Derivation of first-order form for tau'':
    #
    # Term 1: (N/Pr tau')' = (N/Pr) tau'' + (N/Pr)' tau'
    #         = (N/Pr) tau'' + (N'/Pr) tau'
    #
    # Term 2: (S - 1) [N (fp^2)']'
    #         (fp^2)' = 2 fp fpp
    #         [N (fp^2)']' = [N * 2 fp fpp]'
    #                      = N' * 2 fp fpp + N * 2 (fpp^2 + fp fppp)
    #         Note: fppp is available from Liu Eq. 16 (momentum equation)
    #
    # Term 3: (K - 1) [S N (g_cf^2)']'
    #         (g_cf^2)' = 2 g_cf gcf_p
    #         [S N (g_cf^2)']' = [S N * 2 g_cf gcf_p]'
    #                          = S N' * 2 g_cf gcf_p + S N * 2 (gcf_p^2 + g_cf gcf_pp)
    #         Note: gcf_pp is available from Liu Eq. 17 (crossflow momentum equation)
    #
    # Term 4: f [tau' + (S - 1)(fp^2)' + (K - 1) S (g_cf^2)']
    #         = f tau' + f (S - 1) * 2 fp fpp + f (K - 1) S * 2 g_cf gcf_p
    #
    # Solving for tau'':
    # (N/Pr) tau'' = - (N'/Pr) tau'
    #                - (S - 1) [N' * 2 fp fpp + N * 2 (fpp^2 + fp fppp)]
    #                - (K - 1) S [N' * 2 g_cf gcf_p + N * 2 (gcf_p^2 + g_cf gcf_pp)]
    #                - f tau' - f (S - 1) * 2 fp fpp - f (K - 1) S * 2 g_cf gcf_p
    #
    # tau'' = (Pr/N) * [all RHS terms]
    # --------------------------------------------------

    # Note: S and K are already computed above (before momentum equations)
    # K = [1 + (gamma-1)/2 * mach_e_ref^2] / [1 + (gamma-1)/2 * mach_e_ref^2 * cos^2(Lambda)]  [Liu Eq. 21]
    # S = 1 + (gamma-1)/2 * v^2 * mach_e_ref^2 * cos^2(Lambda)  [Liu Eq. 22]

    # N_prime was computed once near the top of this function (reused here).

    # fppp and gcf_pp are exactly the momentum derivatives already solved
    # above (Liu Eqs. 16 and 17). Reuse them directly instead of recomputing
    # so the two copies can never drift out of sync.
    fppp = dy[2]
    gcf_pp = dy[4]

    # Compute derivatives of kinetic energy terms
    fp2_prime = 2.0 * fp * fpp
    gcf2_prime = 2.0 * g_cf * gcf_p

    dy[5] = tau_p

    # Term 1: -(N'/Pr) tau'
    term1 = -(N_prime / problem.prandtl) * tau_p

    # Term 2: -(S - 1) [N' * 2 fp fpp + N * 2 (fpp^2 + fp fppp)]
    term2_bracket = N_prime * fp2_prime + N * 2.0 * (fpp**2 + fp * fppp)
    term2 = -(S - 1.0) * term2_bracket

    # Term 3: -(K - 1) S [N' * 2 g_cf gcf_p + N * 2 (gcf_p^2 + g_cf gcf_pp)]
    term3_bracket = S * (N_prime * gcf2_prime + N * 2.0 * (gcf_p**2 + g_cf * gcf_pp))
    term3 = -(K - 1.0) * term3_bracket

    # Term 4: -f [tau' + (S - 1)(fp^2)' + (K - 1) S (g_cf^2)']
    term4 = -f * (tau_p + (S - 1.0) * fp2_prime + (K - 1.0) * S * gcf2_prime)

    # Solve for tau'': tau'' = (Pr/N) * [sum of all terms]
    dy[6] = (problem.prandtl / N) * (term1 + term2 + term3 + term4)

    return dy
