# Derivation

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

**Stream function** (compressible form):

$$\rho u = \frac{\partial\psi}{\partial y}, \qquad \rho v = -\frac{\partial\psi}{\partial x}$$

**Levy-Lees similarity coordinates**:

$$\xi = \int_0^x \rho_e\mu_e u_e\,dx', \qquad \eta = \frac{u_e}{\sqrt{2\xi}}\int_0^y\rho\,dy'$$

For the power-law $u_e = Cx^m$, the edge quantities $\rho_e$ and $\mu_e$ are constant
(isentropic edge, uniform composition), so $\rho_e\mu_e$ can be taken out of the integral:

$$
\xi = \rho_e\mu_e\int_0^x C x'^m\,dx'
    = \rho_e\mu_e\,C\,\frac{x^{m+1}}{m+1}
    = \frac{\rho_e\mu_e\,\overbrace{Cx^m}^{u_e}\,x}{m+1}
    = \frac{\rho_e\mu_e u_e x}{m+1}
$$

**Hartree parameter**:

$$\beta_H = \frac{2m}{m+1}$$

**Chapman-Rubesin factor**:

$$C = \frac{\rho\mu}{\rho_e\mu_e}$$

**temperature ratio**:

$$
\tau = \frac{T}{T_e}
$$

**Ideal gas law**:

From the ideal gas law $p = \rho R T$, the density is $\rho = p/(RT)$, and at the boundary layer edge $\rho_e = p_e/(RT_e)$.

From the [y-momentum boundary layer equation](../compressible_boundary_layer_equations/2d_equations.md) 

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

To transform any $\partial/\partial x\big|_y$ or $\partial/\partial y\big|_x$ term, we need
the partial derivatives of $\xi$ and $\eta$ with respect to $x$ and $y$. Both follow from the
**Leibniz integral rule**:

$$\frac{d}{dx}\!\left(\int_{a(x)}^{b(x)} f(x,t)\,dt\right)
= f(x,b(x))\frac{d b}{dx} - f(x,a(x))\frac{d a}{dx}
+ \int_{a(x)}^{b(x)}\frac{\partial f}{\partial x}\,dt$$

For $\xi$: the lower limit is constant ($a=0$), the upper limit is $b=x$, and the integrand
$\rho_e\mu_e u_e$ depends only on the dummy variable $x'$ (not on $x$ explicitly):

$$
\frac{\partial\xi}{\partial x} = \frac{\partial}{\partial x}\left(\int_0^x \rho_e\mu_e u_e\,dx'\right) = \rho_e\mu_e u_e
$$

For $\eta$: differentiating with respect to $y$, the lower limit is constant ($a=0$), the upper
limit is $b=y$, and the integrand $\rho$ depends only on $y'$:

$$
\frac{\partial\eta}{\partial y} = \frac{u_e}{\sqrt{2\xi}}\frac{\partial}{\partial y}\left(\int_0^y \rho\,dy'\right) = \frac{\rho u_e}{\sqrt{2\xi}}
$$

$$\frac{\partial\eta}{\partial x}\bigg|_y = \frac{\partial}{\partial x}\!\left(\frac{u_e}{\sqrt{2\xi}}\int_0^y\rho\,dy'\right)$$

### Transformation operators

The change of variables is $(x, y) \to (\xi, \eta)$. For any function $F(\xi, \eta)$, the
chain rule gives:

$$\frac{\partial F}{\partial y}\bigg|_x
= \frac{\partial F}{\partial \xi}\bigg|_\eta \underbrace{\frac{\partial \xi}{\partial y}\bigg|_x}_{=\,0}
+ \frac{\partial F}{\partial \eta}\bigg|_\xi \frac{\partial \eta}{\partial y}\bigg|_x
= \frac{\rho u_e}{\sqrt{2\xi}}\frac{\partial F}{\partial \eta}\bigg|_\xi$$

The $\partial\xi/\partial y\big|_x = 0$ because $\xi = \int_0^x \rho_e\mu_e u_e\,dx'$ contains
no $y$-dependence.

$$\frac{\partial F}{\partial x}\bigg|_y
= \frac{\partial F}{\partial \xi}\bigg|_\eta \frac{\partial \xi}{\partial x}\bigg|_y
+ \frac{\partial F}{\partial \eta}\bigg|_\xi \frac{\partial \eta}{\partial x}\bigg|_y$$

### Streamwise velocity

From $\rho u = \partial\psi/\partial y$ and the $\partial/\partial y$ operator:

$$
\rho u = \frac{\partial\psi}{\partial y}
= \frac{\partial\overbrace{\psi}^{\sqrt{2\xi} f(\eta)}}{\partial \eta} \overbrace{\frac{\partial \eta}{\partial y}}^{\frac{\rho u_e}{\sqrt{2\xi}}}
= \sqrt{2\xi} \overbrace{\frac{\partial f(\eta)}{\partial \eta}}^{f'(\eta)} \frac{\rho u_e}{\sqrt{2\xi}}
= \rho\,u_e f'(\eta)
\qquad\Longrightarrow\qquad u = u_e f'(\eta)
$$

The streamwise derivative of $u$ follows the chain rule on $u = u_e(\xi)\,f'(\eta)$.
Since $u_e$ depends on $x$ only through $\xi$, the two $u_e$ chain-rule factors collapse
immediately: $(\partial u_e/\partial\xi)\,(\partial\xi/\partial x) = du_e/dx$:

$$
\frac{\partial \overbrace{u}^{u_e f'(\eta)}}{\partial x}\bigg|_y
= \overbrace{\frac{du_e}{dx}}^{mu_e/x} f'(\eta)
+ u_e\,\overbrace{\frac{\partial f'(\eta)}{\partial \eta}\bigg|_\xi}^{f''(\eta)}\,
  \frac{\partial \eta}{\partial x}\bigg|_y
= \frac{m u_e}{x}\,f'(\eta) + u_e\,f''(\eta)\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

The wall-normal derivative of $u$:

$$
\frac{\partial \overbrace{u}^{u_e f'(\eta)}}{\partial y} 
= \frac{\partial \left(u_e f'(\eta)\right)}{\partial\eta} \overbrace{\frac{\partial\eta}{\partial y}}^{\frac{\rho u_e}{\sqrt{2 \xi}}}
= \frac{\rho u_e}{\sqrt{2\xi}}\frac{\partial(u_e f'(\eta))}{\partial\eta}
= \frac{\rho u_e^2}{\sqrt{2\xi}}f''(\eta)
$$

### Continuity

$$
\frac{\partial(\overbrace{\rho u}^{\frac{\partial\psi}{\partial y}})}{\partial x} 
+ \frac{\partial(\overbrace{\rho v}^{-\frac{\partial\psi}{\partial x}})}{\partial y}
= 
\frac{\partial^2\psi}{\partial x\,\partial y} 
- \frac{\partial^2\psi}{\partial y\,\partial x} = 0
$$

Continuity is satisfied identically.

### x-momentum

Using the expression derived above the transformed terms are:

**$\rho u\frac{\partial u}{\partial x}$ term.**

$$\rho u\frac{\partial u}{\partial x} = \rho u_e f' \left(\frac{m u_e}{x}f' + u_e f'' \frac{\partial\eta}{\partial x}\right)
= \frac{\rho u_e^2}{x}\!\left(m f'^2 + x f'' f'\frac{\partial\eta}{\partial x}\right)$$

**$\rho v \frac{\partial u}{\partial y}$ term.** Expand $\rho v = -\partial\psi/\partial x$
by applying the product rule to $\psi = \sqrt{2\xi}\,f(\eta)$:

$$
\rho v = -\frac{\partial \overbrace{\psi}^{\sqrt{2\xi}\,f(\eta)}}{\partial x}\bigg|_y
= -\left[
    \overbrace{\frac{\partial\sqrt{2\xi}}{\partial x}}^{\rho_e\mu_e u_e/\sqrt{2\xi}}\,f
    + \sqrt{2\xi}\;\overbrace{\frac{\partial f(\eta)}{\partial x}\bigg|_y}^{f'\,\partial\eta/\partial x\big|_y}
  \right]
= -\frac{\rho_e\mu_e u_e}{\sqrt{2\xi}}\,f - \sqrt{2\xi}\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

Multiplying by $\partial u/\partial y = \rho u_e^2 f''/\sqrt{2\xi}$:

$$
\rho v\,\frac{\partial u}{\partial y}
= \left(-\frac{\rho_e\mu_e u_e}{\sqrt{2\xi}}\,f - \sqrt{2\xi}\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y\right)
  \frac{\rho u_e^2}{\sqrt{2\xi}}\,f''
= -\frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,ff''
  - \rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**Combined convective term.** Adding both terms the $\partial\eta/\partial x$ pieces cancel:

$$
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \frac{m\rho u_e^2}{x}\,f'^2
+ \cancel{\rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}}
- \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,ff''
- \cancel{\rho u_e^2\,f'f''\,\frac{\partial\eta}{\partial x}}
= \frac{m\rho u_e^2}{x}\,f'^2 - \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,ff''
$$

Recast $\frac{m}{x}$:

$$
\frac{m}{\underbrace{x}_{\frac{\xi (m+1)}{\rho_e \mu_e u_e}}} 
= \frac{m\rho_e\mu_e u_e}{\xi(m+1)} 
= \overbrace{\frac{m}{(m+1)}}^{\beta_H/2} \frac{\rho_e \mu_e u_e}{\xi}
= \frac{\beta_H\rho_e\mu_e u_e}{2\xi}
$$

Substituting:

$$
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \overbrace{\frac{m}{x}}^{\frac{\beta_H\rho_e\mu_e u_e}{2\xi}} \rho u_e^2 f'^2
  - \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,ff''
= \frac{\beta_H\rho\rho_e\mu_e u_e^3}{2\xi} f'^2 - \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,ff''
$$

$$\boxed{
\rho u\,\frac{\partial u}{\partial x} + \rho v\,\frac{\partial u}{\partial y}
= \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\!\left(\beta_H f'^2 - ff''\right)
}$$

**Pressure term.** 

$$
-\frac{dp}{dx} 
= \rho_e u_e^2\frac{m}{x} 
= \rho_e u_e^2 \frac{\beta_H \rho_e \mu_e u_e}{2 \xi}
= \overbrace{\rho_e}^{\rho \tau} \frac{\rho_e \mu_e u_e^3}{2\xi} \beta_H
= \frac{\rho\rho_e\mu_e u_e^3}{2\xi} \tau \beta_H
$$

$$\boxed{-\frac{dp}{dx} = \frac{\rho\rho_e\mu_e u_e^3}{2\xi} \tau \beta_H}$$

**Viscous term.** 

$$
\mu\frac{\partial \overbrace{u}^{u_e f'}}{\partial y}
= \overbrace{\mu}^{C\rho_e\mu_e/\rho} \cdot \overbrace{\frac{\partial u}{\partial y}}^{\rho u_e^2 f''/\sqrt{2\xi}}
= \frac{C\rho_e\mu_e u_e^2}{\sqrt{2\xi}}\,f''
$$

Now apply the outer $\partial/\partial y = (\rho u_e/\sqrt{2\xi})\,\partial/\partial\eta$,
noting that $\rho_e\mu_e u_e^2/\sqrt{2\xi}$ does not depend on $\eta$:

$$
\frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
= \overbrace{\frac{\partial}{\partial y}}^{\frac{\rho u_e}{\sqrt{2\xi}}\partial/\partial\eta}
  \left(\overbrace{\mu\frac{\partial u}{\partial y}}^{C\rho_e\mu_e u_e^2 f''/\sqrt{2\xi}}\right)
= \frac{\rho u_e}{\sqrt{2\xi}}\cdot\frac{\rho_e\mu_e u_e^2}{\sqrt{2\xi}}\,(Cf'')'
$$

$$\boxed{\frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
= \frac{\rho\rho_e\mu_e u_e^3}{2\xi}\,(Cf'')'}$$

**Assembly.** 

The x-momentum equation $\rho u\,\partial u/\partial x + \rho v\,\partial u/\partial y = -dp/dx + \partial(\mu\,\partial u/\partial y)/\partial y$ becomes:

$$
\cancel{\frac{\rho\rho_e\mu_e u_e^3}{2\xi}}\!\left(\beta_H f'^2 - ff''\right)
= \cancel{\frac{\rho\rho_e\mu_e u_e^3}{2\xi}}\beta_H\tau
+ \cancel{\frac{\rho\rho_e\mu_e u_e^3}{2\xi}}\,(Cf'')'
$$

Dividing through by $\rho\rho_e\mu_e u_e^3/(2\xi)$ and rearranging:

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
  \overbrace{\frac{\partial\eta}{\partial y}}^{\rho u_e/\sqrt{2\xi}}
= \frac{\rho u_e T_e}{\sqrt{2\xi}}\,\tau'
$$

$$
\frac{\partial \overbrace{T}^{T_e\,\tau(\eta)}}{\partial x}\bigg|_y
= T_e\,\overbrace{\frac{\partial\tau}{\partial\eta}}^{\tau'}\,\frac{\partial\eta}{\partial x}\bigg|_y
= T_e\,\tau'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

**Convective term.** 

*Streamwise*: 

$$
\rho c_p u\,\partial T/\partial x = \rho c_p u_e f' \cdot T_e\tau'\,\partial\eta/\partial x\big|_y
$$

*Wall-normal*: using $\rho v$ from the x-momentum section:

$$
c_p(\rho v)\frac{\partial T}{\partial y}
= c_p\!\left(-\frac{\rho_e\mu_e u_e}{\sqrt{2\xi}}\,f
    - \sqrt{2\xi}\,f'\,\frac{\partial\eta}{\partial x}\bigg|_y\right)
  \frac{\rho u_e T_e\,\tau'}{\sqrt{2\xi}}
= -\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\,f\tau'
  - \rho u_e c_p T_e f'\tau'\,\frac{\partial\eta}{\partial x}\bigg|_y
$$

Combining the above terms

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= \cancel{\rho c_p u_e T_e f'\tau'\frac{\partial\eta}{\partial x}}
  - \frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\,f\tau'
  - \cancel{\rho c_p u_e T_e f'\tau'\frac{\partial\eta}{\partial x}}
$$

$$\boxed{\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= -\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\,f\tau'}$$

**Diffusion term.** 

$$
k\frac{\partial T}{\partial y}
= \overbrace{k}^{C\rho_e\mu_e c_p/(\rho\,\mathrm{Pr})} \cdot
  \overbrace{\frac{\partial T}{\partial y}}^{\rho u_e T_e\tau'/\sqrt{2\xi}}
= \frac{C\rho_e\mu_e u_e c_p T_e}{\mathrm{Pr}\sqrt{2\xi}}\,\tau'
$$

Applying the outer $\partial/\partial y = (\rho u_e/\sqrt{2\xi})\,\partial/\partial\eta$:

$$
\frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
= \overbrace{\frac{\partial}{\partial y}}^{\frac{\rho u_e}{\sqrt{2\xi}}\partial/\partial\eta}
  \left(\overbrace{k\frac{\partial T}{\partial y}}^{C\rho_e\mu_e u_e c_p T_e\tau'/(\mathrm{Pr}\sqrt{2\xi})}\right)
= \frac{\rho u_e}{\sqrt{2\xi}}\cdot\frac{\rho_e\mu_e u_e c_p T_e}{\mathrm{Pr}\sqrt{2\xi}}\,(C\tau')'
$$

$$\boxed{\frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
= \frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\left(\frac{C}{\mathrm{Pr}}\tau'\right)'}$$

**Pressure work term.**

$$
u\frac{dp}{dx}
= \overbrace{u_e f'}^{u}\cdot\left(-\rho_e u_e^2\overbrace{\frac{m}{x}}^{\beta_H\rho_e\mu_e u_e/(2\xi)}\right)
= -\frac{\beta_H\rho_e^2\mu_e u_e^4 f'}{2\xi}
= -\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{\beta_H\tau u_e^2 f'}{c_p T_e}}^{(\gamma-1)M_e^2\beta_H\tau f'}
$$

$$\boxed{u\frac{dp}{dx}
= -\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\,(\gamma-1)M_e^2\,\beta_H\tau f'}$$

**Dissipation term.**

$$
\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \overbrace{\frac{C\rho_e\mu_e}{\rho}}^{\mu}
  \left(\overbrace{\frac{\rho u_e^2 f''}{\sqrt{2\xi}}}^{\partial u/\partial y}\right)^{\!2}
= \frac{C\rho\rho_e\mu_e u_e^4 f''^2}{2\xi}
= \frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\cdot
  \overbrace{\frac{C u_e^2 f''^2}{c_p T_e}}^{(\gamma-1)M_e^2\,Cf''^2}
$$

$$\boxed{\mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
= \frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}\,(\gamma-1)M_e^2\,Cf''^2}$$

**Assembly.** Every term carries $\rho\rho_e\mu_e u_e^2 c_p T_e/(2\xi)$:

$$
-\cancel{\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}}\,f\tau'
= -\cancel{\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}}\,(\gamma-1)M_e^2\beta_H\tau f'
+ \cancel{\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}}\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
+ \cancel{\frac{\rho\rho_e\mu_e u_e^2 c_p T_e}{2\xi}}\,(\gamma-1)M_e^2 Cf''^2
$$

Dividing through and rearranging:

!!! info ""
    $$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau'
    + (\gamma-1)M_e^2\!\left[Cf''^2 - \beta_H\tau f'\right] = 0$$




