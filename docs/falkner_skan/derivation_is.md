# Derivation (Illingworth-Stewartson)

!!! note
    This derivation uses the Illingworth-Stewartson (IS) transformation.
    For the equivalent derivation using the Levy-Lees transformation,
    see [Derivation (Levy-Lees)](derivation.md).

Starting from the steady [2D compressible BL equations](../compressible_boundary_layer_equations/2d_equations.md):

$$
\frac{\partial(\rho u)}{\partial x} + \frac{\partial(\rho v)}{\partial y} = 0
$$

$$
\rho\!\left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}\right)
= -\frac{dp}{dx} + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
$$

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= u\frac{dp}{dx} + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
+ \mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
$$

### Similarity Ansatz

**Edge velocity power law**:

$$u_e(x) = C x^m$$

**Dimensionless stream function**:

$$\psi = \sqrt{2\xi}\,f(\eta)$$

### Outer flow

The outer inviscid flow satisfies the Euler x-momentum equation (see
[2D BL equations](../compressible_boundary_layer_equations/2d_equations.md)):

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

The dimensionless stream function is:

$$\psi = \sqrt{2\xi}\,f(\eta), \qquad f' = \frac{u}{u_e}$$

**Hartree parameter**:

$$\beta_H = \frac{2m}{m+1}$$

**Chapman-Rubesin factor** (IS, referenced to freestream):

$$C = \frac{\rho\mu}{\rho_\infty\mu_\infty}$$

**Temperature ratio**:

$$\tau = \frac{T}{T_e}$$

**Ideal gas law**:

From the ideal gas law $p = \rho R T$, the density is $\rho = p/(RT)$, and at the boundary
layer edge $\rho_e = p_e/(RT_e)$. From the
[y-momentum boundary layer equation](../compressible_boundary_layer_equations/2d_equations.md)

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

**Edge Mach number**:

$$
M_e = \frac{u_e}{\sqrt{\gamma R T_e}}
$$

With $c_p = \gamma R/(\gamma-1)$ it follows that

$$
\frac{u_e^2}{c_p T_e} = (\gamma-1)M_e^2
$$

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

### Energy

The energy equation is:

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= u\frac{dp}{dx}
+ \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
+ \mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
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

**Pressure work term.**

$$
u\frac{dp}{dx}
= \overbrace{u_e f'}^{u}\cdot\left(-\rho_e u_e^2\overbrace{\frac{m}{x}}^{\beta_H\rho_\infty\mu_\infty u_e/(2\xi)}\right)
= -\frac{\beta_H\rho_e\rho_\infty\mu_\infty u_e^4 f'}{2\xi}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{\beta_H\tau u_e^2 f'}{c_p T_e}}^{(\gamma-1)M_e^2\beta_H\tau f'}
$$

$$\boxed{u\frac{dp}{dx}
= -\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,(\gamma-1)M_e^2\,\beta_H\tau f'}$$

**Dissipation term.**

$$
\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \overbrace{\frac{C\rho_\infty\mu_\infty}{\rho}}^{\mu}
  \left(\overbrace{\frac{\rho u_e^2 f''}{\rho_\infty\sqrt{2\xi}}}^{\partial u/\partial y}\right)^{\!2}
= \frac{C\rho\mu_\infty u_e^4 f''^2}{\rho_\infty\cdot 2\xi}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{C u_e^2 f''^2}{c_p T_e}}^{(\gamma-1)M_e^2\,Cf''^2}
$$

$$\boxed{\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}\,(\gamma-1)M_e^2\,Cf''^2}$$

**Assembly.** Every term carries $\rho\rho_\infty\mu_\infty u_e^2 c_p T_e/(2\xi)$:

$$
-\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,f\tau'
= -\cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,(\gamma-1)M_e^2\beta_H\tau f'
+ \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
+ \cancel{\frac{\rho\rho_\infty\mu_\infty u_e^2 c_p T_e}{2\xi}}\,(\gamma-1)M_e^2 Cf''^2
$$

Dividing through and rearranging:

!!! info ""
    $$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau'
    + (\gamma-1)M_e^2\!\left[Cf''^2 - \beta_H\tau f'\right] = 0$$

---

