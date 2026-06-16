# Falkner-Skan-Cooke Equations

The Falkner-Skan-Cooke (FSC) equations obtained by reducing the [quasi-3D compressible boundary layer equations](../boundary_layer_equations/quasi_3d_equations.md) to an ODE system [@liu_2021]. 

??? abstract "Inherited assumptions"
    - Steady, laminar flow with three velocity components $(u, v, w)$
    - Homogeneous flow in the spanwise direction ($\partial/\partial z = 0$)
    - Calorically perfect gas ($\gamma$, $c_p$, $\mathrm{Pr}$ constant)
    - Dynamic viscosity via Sutherland's law
    - Thin boundary layer ($\delta \ll L$)

## Similarity Ansatz

The edge velocity follows a power law in the transformed coordinate:

$$U_e(\tilde{\xi}) = C_1\cdot\tilde{\xi}^m$$

The Falkner-Skan similarity variable and stream function are:

$$\eta = \bar{\eta}\sqrt{\frac{m+1}{2}\frac{U_e}{\nu_{e0}\tilde{\xi}}},
\qquad
\psi = f(\eta)\sqrt{\frac{2\nu_{e0}U_e\tilde{\xi}}{m+1}}$$

where $\nu_{e0} = \mu_{e0}/\rho_{e0}$ is the kinematic viscosity at the reference point.
The normalized velocity components are [@liu_2021]:

$$f'(\eta) = \frac{U}{U_e} = \frac{u}{u_e}, \qquad g(\eta) = \frac{w}{w_e}, \qquad \tau(\eta) = \frac{T}{T_e}$$

## Illingworth-Stewartson Transformation

The IS transformation removes density from the governing equations [@liu_2021]:

$$d\tilde{\xi} = \frac{\mu_e a_e \rho_e}{\mu_{e0} a_{e0} \rho_{e0}}\,dx,
\qquad
d\bar{\eta} = \frac{a_e \rho}{a_{e0} \rho_{e0}}\,dy$$

where $a$ is the local sound speed and subscript $e0$ denotes conditions at the
reference edge point. The stream function $\psi$ in the transformed space is defined as:

$$\frac{\partial\psi}{\partial\bar{\eta}} = U, \qquad \frac{\partial\psi}{\partial\tilde{\xi}} = -V$$

## Definitions

**Hartree parameter**:

$$\beta_H = \frac{2m}{m+1}$$

**Chapman-Rubesin factor**:

$$
C \equiv \frac{\rho\mu}{\rho_e\mu_e}
$$

**Compressibility parameters**:

$$
K = \frac{1 + \dfrac{\gamma-1}{2}M_e^2}{1 + \dfrac{\gamma-1}{2}M_e^2\cos^2\!\Lambda},
\qquad
S = 1 + \frac{\gamma-1}{2}M_e^2\cos^2\!\Lambda
$$

where $\Lambda$ is the local swept angle and $M_e$ the local edge Mach number.

!!! warning "Deviation from Liu (2021)"
    [@liu_2021] defines $K$ and $S$ in terms of a reference Mach number $Ma_{e,\mathrm{ref}}$ and a
    streamwise parameter $\chi = (\tilde{\xi}/\tilde{\xi}_\mathrm{ref})^m$ that tracks the variation
    of edge conditions along the surface. Here we adopt a **locally self-similar** formulation, which
    sets $\chi = 1$ and consequently $Ma_{e,\mathrm{ref}} = M_e$. The parameters $\beta_H$, $K$, $S$,
    and $T_e$ are then known inputs evaluated from the local edge conditions at each station.


## ODE System

**x-momentum** [@liu_2021]:

$$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

??? note "Compare to Falkner-Skan"
    The [Falkner-Skan x-momentum](../falkner_skan/index.md#ode-system) is:

    $$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

**z-momentum** (crossflow) [@liu_2021]:

$$(Cg')' + fg' = 0$$

**Energy** [@liu_2021]:

$$\begin{aligned}
&\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
+ (S-1)\left(C(f'^2)'\right)'
+ (K-1)S\left(C(g^2)'\right)' \\
&\quad + f\!\left[\tau' + (S-1)(f'^2)' + (K-1)S(g^2)'\right] = 0
\end{aligned}$$

## Boundary Conditions

**Wall** ($\eta = 0$):

$$f = 0, \qquad f' = 0, \qquad g = 0$$

- **Isothermal**: $\tau(0) = T_w/T_e$ (prescribed)
- **Adiabatic**: $\tau'(0) = 0$

**Edge** ($\eta \to \infty$):

$$f' = 1, \qquad g = 1, \qquad \tau = 1$$

