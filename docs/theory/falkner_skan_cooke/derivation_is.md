# Derivation (Illingworth-Stewartson)

!!! note
    This derivation uses the Illingworth-Stewartson transformation following [@liu_2021].
    For the equivalent derivation using the Lévy-Lees transformation,
    see [Derivation (Lévy-Lees)](derivation.md).

Starting from the [quasi-3D compressible BL equations](../boundary_layer_equations/quasi_3d_equations.md):

$$
\frac{\partial(\rho u)}{\partial x} + \frac{\partial(\rho v)}{\partial y} = 0
$$

$$
\rho\!\left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}\right)
= -\frac{\partial p}{\partial x} + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
$$

$$
\rho\!\left(u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y}\right)
= \frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y}\right)
$$

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= u\frac{\partial p}{\partial x}
+ \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
+ \mu\!\left[\left(\frac{\partial u}{\partial y}\right)^{\!2}
+ \left(\frac{\partial w}{\partial y}\right)^{\!2}\right]
$$

The z-momentum equation carries **no** pressure-gradient term (dropped under the quasi-3D
assumption $\partial p/\partial z = 0$), and the energy equation carries an **extra**
crossflow dissipation term $\mu(\partial w/\partial y)^2$ relative to the 2D case.

### Similarity Ansatz

**Edge velocity power law** (chordwise component):

$$u_e(x) = C x^m$$

The spanwise edge velocity $w_e$ is **constant** along the chord (Cooke's independence
principle: the spanwise momentum equation decouples and $w_e$ is a prescribed constant).

**Dimensionless stream function** and normalized dependent variables [@liu_2021]:

$$\psi = \sqrt{2\xi}\,f(\eta), \qquad
f'(\eta) = \frac{u}{u_e}, \qquad
g(\eta) = \frac{w}{w_e}, \qquad
\tau(\eta) = \frac{T}{T_e}$$

### Outer flow

The outer inviscid flow satisfies the Euler x-momentum equation (see
[quasi-3D BL equations](../boundary_layer_equations/quasi_3d_equations.md)):

$$-\frac{dp}{dx} = \rho_e u_e \frac{du_e}{dx}$$

For the power-law edge velocity $u_e = Cx^m$:

$$
\frac{du_e}{dx} = \frac{d}{dx}\left(Cx^{m}\right)
= m C x^{m-1} = \frac{m}{x} \underbrace{C x^m}_{u_e} = \frac{m u_e}{x}
$$

so the pressure gradient becomes:

$$-\frac{dp}{dx} = \rho_e u_e^2\frac{m}{x}$$

### Definitions

**Transformed coordinates** (IS physical-to-transformed map). The IS transformation references freestream constants $\rho_\infty$, $\mu_\infty$, so $\xi$ requires no assumption on the streamwise variation of edge quantities:

$$\bar{x} = \int_0^x \frac{\rho_e \mu_e}{\rho_\infty \mu_\infty}\,dx', \qquad
\bar{y} = \frac{1}{\rho_\infty}\int_0^y \rho\,dy'$$

In the transformed space $(\bar{x}, \bar{y})$, the continuity equation takes the
incompressible form (unit density). The stream function $\psi$ is therefore defined
**without** a density weighting:

$$\frac{\partial\psi}{\partial \bar{y}} = u, \qquad
\frac{\partial\psi}{\partial \bar{x}} = -\bar{v}$$

The similarity coordinates follow from $\bar{x}$ and $\bar{y}$:

$$\xi = \int_0^{\bar{x}} u_e\,d\bar{x}' = \int_0^x \frac{\rho_e \mu_e u_e}{\rho_\infty\mu_\infty}\,dx',
\qquad
\eta = \frac{u_e}{\sqrt{2\xi}}\,\bar{y} = \frac{u_e}{\rho_\infty\sqrt{2\xi}}\int_0^y\rho\,dy'$$

**Hartree parameter**:

$$\beta_H = \frac{2m}{m+1}$$

**Chapman-Rubesin factor** (IS, referenced to freestream):

$$C = \frac{\rho\mu}{\rho_\infty\mu_\infty}$$

For a uniform freestream ($\rho_e = \rho_\infty$, $\mu_e = \mu_\infty$) this coincides with
the edge-referenced form $C = \rho\mu/(\rho_e\mu_e)$ used in the
[ODE system](index.md#ode-system) and [reduction of order](reduction_of_order.md).

**Temperature ratio**:

$$\tau = \frac{T}{T_e}$$

**Ideal gas law**:

From the ideal gas law $p = \rho R T$, the density is $\rho = p/(RT)$, and at the boundary
layer edge $\rho_e = p_e/(RT_e)$. From the
[y-momentum boundary layer equation](../boundary_layer_equations/quasi_3d_equations.md)

$$
\frac{\partial p}{\partial y} = 0 \rightarrow p = p_e
$$

Therefore:

$$
\frac{\rho}{\rho_e}
= \frac{p/(RT)}{p_e/(RT_e)}
= \frac{T_e}{T}
= \frac{1}{\tau}
$$

**Sweep decomposition and edge Mach number**:

The total edge velocity $Q_e$ splits into chordwise and spanwise components through the
local swept angle $\Lambda$:

$$u_e = Q_e\cos\Lambda, \qquad w_e = Q_e\sin\Lambda, \qquad
M_e = \frac{Q_e}{\sqrt{\gamma R T_e}}$$

With $c_p = \gamma R/(\gamma-1)$ it follows that

$$
\frac{u_e^2}{c_p T_e} = (\gamma-1)M_e^2\cos^2\!\Lambda,
\qquad
\frac{w_e^2}{c_p T_e} = (\gamma-1)M_e^2\sin^2\!\Lambda
$$

**Compressibility parameters**:

$$
K = \frac{1 + \dfrac{\gamma-1}{2}M_e^2}{1 + \dfrac{\gamma-1}{2}M_e^2\cos^2\!\Lambda},
\qquad
S = 1 + \frac{\gamma-1}{2}M_e^2\cos^2\!\Lambda
$$

where $\Lambda$ is the local swept angle and $M_e$ the local edge Mach number. These
combine with the relations above into the compact identities used throughout the energy
assembly:

$$
2(S-1) = (\gamma-1)M_e^2\cos^2\!\Lambda = \frac{u_e^2}{c_p T_e},
\qquad
2(K-1)S = (\gamma-1)M_e^2\sin^2\!\Lambda = \frac{w_e^2}{c_p T_e}
$$

!!! warning "Deviation from Liu (2021)"
    [@liu_2021] defines $K$ and $S$ in terms of a reference Mach number $Ma_{e,\mathrm{ref}}$ and a
    streamwise parameter $\chi = (\tilde{\xi}/\tilde{\xi}_\mathrm{ref})^m$ that tracks the variation
    of edge conditions along the surface. Here we adopt a **locally self-similar** formulation, which
    sets $\chi = 1$ and consequently $Ma_{e,\mathrm{ref}} = M_e$. The parameters $\beta_H$, $K$, $S$,
    and $T_e$ are then known inputs evaluated from the local edge conditions at each station.

**Prandtl number**:

$$\mathrm{Pr} = \frac{\mu c_p}{k}$$

so that $k = \mu c_p/\mathrm{Pr}$.

### Partial derivatives of the similarity coordinates

By the Leibniz integral rule:

$$\frac{\partial\xi}{\partial x} = \frac{\rho_e\mu_e u_e}{\rho_\infty\mu_\infty}, \qquad
\frac{\partial\eta}{\partial y} = \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}$$

For the power-law $u_e = Cx^m$ with $\rho_e\mu_e = \text{const}$ (uniform freestream, $\rho_e = \rho_\infty$, $\mu_e = \mu_\infty$):

$$
\xi = \frac{\rho_\infty\mu_\infty}{\rho_\infty\mu_\infty}\int_0^x Cx'^m\,dx'
    = \frac{\overbrace{Cx^m}^{u_e}\,x}{m+1}
    = \frac{\rho_\infty\mu_\infty u_e x}{m+1}
$$

so that $m/x = \beta_H\rho_\infty\mu_\infty u_e/(2\xi)$. The $\partial\eta/\partial x$ term cancels in the convective and energy terms.

### Transformation operators

$$\frac{\partial F}{\partial y}\bigg|_x
= \frac{\partial F}{\partial \xi}\bigg|_\eta \underbrace{\frac{\partial \xi}{\partial y}\bigg|_x}_{=\,0}
+ \frac{\partial F}{\partial \eta}\bigg|_\xi \frac{\partial \eta}{\partial y}\bigg|_x
= \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\frac{\partial F}{\partial \eta}\bigg|_\xi$$

$$\frac{\partial F}{\partial x}\bigg|_y
= \frac{\partial F}{\partial \xi}\bigg|_\eta \frac{\partial \xi}{\partial x}\bigg|_y
+ \frac{\partial F}{\partial \eta}\bigg|_\xi \frac{\partial \eta}{\partial x}\bigg|_y$$

### Streamwise velocity

From $\partial\psi/\partial\bar{y} = u$ it follows that $u = u_e f'(\eta)$. The wall-normal derivative:

$$
\frac{\partial \overbrace{u}^{u_e f'(\eta)}}{\partial y}
= \frac{\partial \left(u_e f'(\eta)\right)}{\partial\eta} \overbrace{\frac{\partial\eta}{\partial y}}^{\frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}}
= \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\frac{\partial(u_e f'(\eta))}{\partial\eta}
= \frac{\rho u_e^2}{\rho_\infty\sqrt{2\xi}}\,f''(\eta)
$$

The streamwise derivative:

$$
\frac{\partial \overbrace{u}^{u_e f'(\eta)}}{\partial x}\bigg|_y
= \overbrace{\frac{du_e}{dx}}^{mu_e/x} f'(\eta)
+ u_e\,f''(\eta)\,\frac{\partial\eta}{\partial x}\bigg|_y
= \frac{m u_e}{x}\,f'(\eta) + u_e\,f''(\eta)\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

### Continuity

Satisfied identically by construction of $\psi$.

### x-momentum

Using the expressions derived above:

**$\rho u\frac{\partial u}{\partial x}$ term.**

$$\rho u\frac{\partial u}{\partial x} = \rho u_e f' \left(\frac{m u_e}{x}f' + u_e f'' \frac{\partial\eta}{\partial x}\right)
= \frac{\rho u_e^2}{x}\!\left(m f'^2 + x f'' f'\frac{\partial\eta}{\partial x}\right)$$

**$\rho v \frac{\partial u}{\partial y}$ term.** From $\partial\psi/\partial\bar{x} = -\bar{v}$, applying the product rule to $\psi = \sqrt{2\xi}\,f(\eta)$ with $\rho_e = \rho_\infty$, $\mu_e = \mu_\infty$:

$$
\rho v = -\frac{\rho_\infty\mu_\infty u_e}{\sqrt{2\xi}}\,f - \sqrt{2\xi}\,\rho_\infty\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

Multiplying by $\partial u/\partial y = \rho u_e^2 f''/(\rho_\infty\sqrt{2\xi})$:

$$
\rho v\,\frac{\partial u}{\partial y}
= \left(-\frac{\rho_\infty\mu_\infty u_e}{\sqrt{2\xi}}\,f - \sqrt{2\xi}\,\rho_\infty\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y\right)
  \frac{\rho u_e^2 f''}{\rho_\infty\sqrt{2\xi}}
= -\frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,ff''
  - \rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**Combined convective term.** Adding both terms, the $\partial\eta/\partial x$ pieces cancel:

$$
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \frac{m\rho u_e^2}{x}\,f'^2
+ \cancel{\rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}}
- \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,ff''
- \cancel{\rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}}
= \frac{m\rho u_e^2}{x}\,f'^2 - \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,ff''
$$

Recast $\frac{m}{x}$ using $\xi = \rho_\infty\mu_\infty u_e x/(m+1)$:

$$
\frac{m}{\underbrace{x}_{\frac{\xi(m+1)}{\rho_\infty\mu_\infty u_e}}}
= \frac{m\rho_\infty\mu_\infty u_e}{\xi(m+1)}
= \overbrace{\frac{m}{(m+1)}}^{\beta_H/2} \frac{\rho_\infty\mu_\infty u_e}{\xi}
= \frac{\beta_H\rho_\infty\mu_\infty u_e}{2\xi}
$$

Substituting:

$$
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \overbrace{\frac{m}{x}}^{\frac{\beta_H\rho_\infty\mu_\infty u_e}{2\xi}} \rho u_e^2 f'^2
  - \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,ff''
= \frac{\beta_H\rho\rho_\infty\mu_\infty u_e^3}{2\xi} f'^2 - \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,ff''
$$

$$\boxed{
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\!\left(\beta_H f'^2 - ff''\right)
}$$

**Pressure term.**

$$
-\frac{dp}{dx}
= \rho_e u_e^2\frac{m}{x}
= \rho_e u_e^2 \frac{\beta_H \rho_\infty\mu_\infty u_e}{2\xi}
= \overbrace{\rho_e}^{\rho\tau} \frac{\rho_\infty\mu_\infty u_e^3}{2\xi}\beta_H
= \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,\tau\beta_H
$$

$$\boxed{-\frac{dp}{dx} = \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,\tau\beta_H}$$

**Viscous term.**

$$
\mu\frac{\partial \overbrace{u}^{u_e f'}}{\partial y}
= \overbrace{\mu}^{C\rho_\infty\mu_\infty/\rho} \cdot \overbrace{\frac{\partial u}{\partial y}}^{\rho u_e^2 f''/(\rho_\infty\sqrt{2\xi})}
= \frac{C\mu_\infty u_e^2}{\sqrt{2\xi}}\,f''
$$

Applying the outer $\partial/\partial y = (\rho u_e/(\rho_\infty\sqrt{2\xi}))\,\partial/\partial\eta$:

$$
\frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
= \overbrace{\frac{\partial}{\partial y}}^{\frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\partial/\partial\eta}
  \left(\overbrace{\mu\frac{\partial u}{\partial y}}^{C\mu_\infty u_e^2 f''/\sqrt{2\xi}}\right)
= \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\cdot\frac{\rho_\infty\mu_\infty u_e^2}{\sqrt{2\xi}}\,(Cf'')'
$$

$$\boxed{\frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
= \frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}\,(Cf'')'}$$

**Assembly.**

The x-momentum equation becomes:

$$
\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}}\!\left(\beta_H f'^2 - ff''\right)
= \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}}\beta_H\tau
+ \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^3}{2\xi}}\,(Cf'')'
$$

Dividing through by $\rho\rho_\infty\mu_\infty u_e^3/(2\xi)$ and rearranging:

!!! info ""
    $$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

??? note "Compare to Falkner-Skan"
    The x-momentum equation is identical to the
    [Falkner-Skan x-momentum](../falkner_skan/derivation_is.md#x-momentum) —
    the crossflow $w$ does not appear.

### z-momentum

The crossflow momentum equation has no pressure-gradient term:

$$
\rho\!\left(u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y}\right)
= \frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y}\right)
$$

With $w = w_e\,g(\eta)$ and $w_e = \text{const}$, the crossflow is a function of $\eta$ only,
so there is **no** $m u_e/x$ contribution from the streamwise derivative:

$$
\frac{\partial \overbrace{w}^{w_e g(\eta)}}{\partial y}
= w_e\,g'(\eta)\,\overbrace{\frac{\partial\eta}{\partial y}}^{\frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}}
= \frac{\rho u_e w_e}{\rho_\infty\sqrt{2\xi}}\,g'(\eta),
\qquad
\frac{\partial w}{\partial x}\bigg|_y
= w_e\,g'(\eta)\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**$\rho u\frac{\partial w}{\partial x}$ term.**

$$
\rho u\frac{\partial w}{\partial x}
= \rho u_e f'\cdot w_e\,g'\,\frac{\partial\eta}{\partial x}
= \rho u_e w_e\,f'g'\,\frac{\partial\eta}{\partial x}
$$

**$\rho v\frac{\partial w}{\partial y}$ term.** Reusing $\rho v$ from the x-momentum section and multiplying by $\partial w/\partial y = \rho u_e w_e g'/(\rho_\infty\sqrt{2\xi})$:

$$
\rho v\,\frac{\partial w}{\partial y}
= \left(-\frac{\rho_\infty\mu_\infty u_e}{\sqrt{2\xi}}\,f - \sqrt{2\xi}\,\rho_\infty\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y\right)
  \frac{\rho u_e w_e\,g'}{\rho_\infty\sqrt{2\xi}}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}\,fg'
  - \rho u_e w_e\,f'g'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**Combined convective term.** Adding both terms, the $\partial\eta/\partial x$ pieces cancel exactly as in the x-momentum case:

$$
\rho u\,\frac{\partial w}{\partial x} + \rho v\,\frac{\partial w}{\partial y}
= \cancel{\rho u_e w_e\,f'g'\,\frac{\partial\eta}{\partial x}}
- \frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}\,fg'
- \cancel{\rho u_e w_e\,f'g'\,\frac{\partial\eta}{\partial x}}
$$

$$\boxed{
\rho u\,\frac{\partial w}{\partial x} + \rho v\,\frac{\partial w}{\partial y}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}\,fg'
}$$

**Viscous term.**

$$
\mu\frac{\partial \overbrace{w}^{w_e g}}{\partial y}
= \overbrace{\mu}^{C\rho_\infty\mu_\infty/\rho} \cdot \overbrace{\frac{\partial w}{\partial y}}^{\rho u_e w_e g'/(\rho_\infty\sqrt{2\xi})}
= \frac{C\mu_\infty u_e w_e}{\sqrt{2\xi}}\,g'
$$

Applying the outer $\partial/\partial y = (\rho u_e/(\rho_\infty\sqrt{2\xi}))\,\partial/\partial\eta$:

$$
\frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y}\right)
= \overbrace{\frac{\partial}{\partial y}}^{\frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\partial/\partial\eta}
  \left(\overbrace{\mu\frac{\partial w}{\partial y}}^{C\mu_\infty u_e w_e g'/\sqrt{2\xi}}\right)
= \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\cdot\frac{\rho_\infty\mu_\infty u_e w_e}{\sqrt{2\xi}}\,(Cg')'
$$

$$\boxed{\frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y}\right)
= \frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}\,(Cg')'}$$

**Assembly.** Every term carries the common factor $\rho\rho_\infty\mu_\infty u_e^2 w_e/(2\xi)$:

$$
-\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}}\,fg'
= \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 w_e}{2\xi}}\,(Cg')'
$$

Dividing through by $\rho\rho_\infty\mu_\infty u_e^2 w_e/(2\xi)$ and rearranging:

!!! info ""
    $$(Cg')' + fg' = 0$$

### Energy

The energy equation carries an extra crossflow dissipation term:

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= u\frac{\partial p}{\partial x}
+ \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
+ \mu\!\left[\left(\frac{\partial u}{\partial y}\right)^{\!2}
+ \left(\frac{\partial w}{\partial y}\right)^{\!2}\right]
$$

Note $T = T_e\,\tau(\eta) \rightarrow$ function of $\eta$ only

$$
\frac{\partial \overbrace{T}^{T_e\,\tau(\eta)}}{\partial y}
= T_e\,\overbrace{\frac{\partial\tau}{\partial\eta}}^{\tau'}\,
  \overbrace{\frac{\partial\eta}{\partial y}}^{\rho u_e/(\rho_\infty\sqrt{2\xi})}
= \frac{\rho u_e T_e}{\rho_\infty\sqrt{2\xi}}\,\tau'
$$

$$
\frac{\partial \overbrace{T}^{T_e\,\tau(\eta)}}{\partial x}\bigg|_y
= T_e\,\overbrace{\frac{\partial\tau}{\partial\eta}}^{\tau'}\,\frac{\partial\eta}{\partial x}\bigg|_y
= T_e\,\tau'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**Convective term.** Streamwise and wall-normal contributions combined (the $\partial\eta/\partial x$ cross terms cancel):

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= \cancel{\rho c_p u_e T_e f'\tau'\frac{\partial\eta}{\partial x}}
  - \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,f\tau'
  - \cancel{\rho c_p u_e T_e f'\tau'\frac{\partial\eta}{\partial x}}
$$

$$\boxed{\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,f\tau'}$$

**Diffusion term.**

$$
k\frac{\partial T}{\partial y}
= \overbrace{k}^{C\rho_\infty\mu_\infty c_p/(\rho\,\mathrm{Pr})} \cdot
  \overbrace{\frac{\partial T}{\partial y}}^{\rho u_e T_e\tau'/(\rho_\infty\sqrt{2\xi})}
= \frac{C\mu_\infty u_e c_p T_e}{\mathrm{Pr}\sqrt{2\xi}}\,\tau'
$$

Applying the outer $\partial/\partial y = (\rho u_e/(\rho_\infty\sqrt{2\xi}))\,\partial/\partial\eta$:

$$
\frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
= \overbrace{\frac{\partial}{\partial y}}^{\frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\partial/\partial\eta}
  \left(\overbrace{k\frac{\partial T}{\partial y}}^{C\mu_\infty u_e c_p T_e\tau'/(\mathrm{Pr}\sqrt{2\xi})}\right)
= \frac{\rho u_e}{\rho_\infty\sqrt{2\xi}}\cdot\frac{\rho_\infty\mu_\infty u_e c_p T_e}{\mathrm{Pr}\sqrt{2\xi}}\,\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
$$

$$\boxed{\frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\left(\frac{C}{\mathrm{Pr}}\tau'\right)'}$$

**Pressure work term.** Using $u_e^2/(c_p T_e) = 2(S-1)$:

$$
u\frac{dp}{dx}
= \overbrace{u_e f'}^{u}\cdot\left(-\rho_e u_e^2\overbrace{\frac{m}{x}}^{\beta_H\rho_\infty\mu_\infty u_e/(2\xi)}\right)
= -\frac{\beta_H\rho_e\rho_\infty\mu_\infty u_e^4 f'}{2\xi}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{\beta_H\tau u_e^2 f'}{c_p T_e}}^{2(S-1)\beta_H\tau f'}
$$

$$\boxed{u\frac{dp}{dx}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,2(S-1)\,\beta_H\tau f'}$$

**Streamwise dissipation term.** Using $u_e^2/(c_p T_e) = 2(S-1)$:

$$
\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \overbrace{\frac{C\rho_\infty\mu_\infty}{\rho}}^{\mu}
  \left(\overbrace{\frac{\rho u_e^2 f''}{\rho_\infty\sqrt{2\xi}}}^{\partial u/\partial y}\right)^{\!2}
= \frac{C\rho\mu_\infty u_e^4 f''^2}{\rho_\infty\cdot 2\xi}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{C u_e^2 f''^2}{c_p T_e}}^{2(S-1)\,Cf''^2}
$$

$$\boxed{\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,2(S-1)\,Cf''^2}$$

**Crossflow dissipation term.** Using $w_e^2/(c_p T_e) = 2(K-1)S$:

$$
\mu\!\left(\frac{\partial w}{\partial y}\right)^{\!2}
= \overbrace{\frac{C\rho_\infty\mu_\infty}{\rho}}^{\mu}
  \left(\overbrace{\frac{\rho u_e w_e g'}{\rho_\infty\sqrt{2\xi}}}^{\partial w/\partial y}\right)^{\!2}
= \frac{C\rho\mu_\infty u_e^2 w_e^2 g'^2}{\rho_\infty\cdot 2\xi}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{C w_e^2 g'^2}{c_p T_e}}^{2(K-1)S\,Cg'^2}
$$

$$\boxed{\mu\!\left(\frac{\partial w}{\partial y}\right)^{\!2}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,2(K-1)S\,Cg'^2}$$

**Assembly.** Every term carries the common factor $\rho\rho_\infty\mu_\infty u_e^2 c_p T_e/(2\xi)$:

$$
\begin{aligned}
-\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,f\tau'
&= -\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,2(S-1)\beta_H\tau f'
+ \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\left(\frac{C}{\mathrm{Pr}}\tau'\right)' \\
&\quad + \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,2(S-1)Cf''^2
+ \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,2(K-1)S\,Cg'^2
\end{aligned}
$$

Dividing through and rearranging gives the directly-derived **static-temperature** energy equation:

!!! info ""
    $$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau'
    + 2(S-1)\!\left[Cf''^2 - \beta_H\tau f'\right]
    + 2(K-1)S\,Cg'^2 = 0$$

This is the FSC analog of the
[Falkner-Skan energy equation](../falkner_skan/index.md#ode-system): setting $w_e = 0$ and
$\Lambda = 0$ makes $2(S-1) \to (\gamma-1)M_e^2$ and drops the crossflow term, recovering
$(C/\mathrm{Pr}\,\tau')' + f\tau' + (\gamma-1)M_e^2[Cf''^2 - \beta_H\tau f'] = 0$.

!!! danger "Discrepancy with the stated ODE system"
    The [FSC ODE system](index.md#ode-system) and
    [reduction of order](reduction_of_order.md) list the energy equation in
    **total-temperature (stagnation-enthalpy)** form:

    $$\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
    + (S-1)\left(C(f'^2)'\right)'
    + (K-1)S\left(C(g^2)'\right)'
    + f\!\left[\tau' + (S-1)(f'^2)' + (K-1)S(g^2)'\right] = 0$$

    Since $(S-1)f'^2 = u^2/(2c_pT_e)$ and $(K-1)Sg^2 = w^2/(2c_pT_e)$, the convected
    bracket $\tau + (S-1)f'^2 + (K-1)Sg^2 = T_0/T_e$ is the total-temperature ratio, so
    this is a stagnation-enthalpy reorganization of the energy equation.

    Expanding that form and substituting the x-momentum ODE
    $(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$ and the z-momentum ODE
    $(Cg')' + fg' = 0$ reduces it to

    $$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau'
    + 2(S-1)\!\left[Cf''^2 - \beta_H\tau f'\right]
    + 2(K-1)S\,Cg'^2
    \;+\; 2(S-1)\beta_H f'^3 = 0$$

    i.e. it differs from the directly-derived static-temperature result above by the term

    $$2(S-1)\beta_H f'^3 = (\gamma-1)M_e^2\cos^2\!\Lambda\;\beta_H f'^3,$$

    which does **not** vanish on solutions (it is not a multiple of either momentum
    residual). The boxed static-temperature equation is the form that follows directly from
    the quasi-3D energy equation and is consistent with the convention used for the
    [Falkner-Skan energy equation](../falkner_skan/index.md#ode-system). The
    stagnation-enthalpy form in `index.md` carries this extra term; per the task this is
    reported rather than absorbed by forcing the algebra.

### Boundary Conditions

**Wall** ($\eta = 0$):

$$f = 0, \qquad f' = 0, \qquad g = 0$$

- **Isothermal**: $\tau(0) = T_w/T_e$ (prescribed)
- **Adiabatic**: $\tau'(0) = 0$

**Edge** ($\eta \to \infty$):

$$f' = 1, \qquad g = 1, \qquad \tau = 1$$
