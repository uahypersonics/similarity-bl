## Assumptions

- Steady, axisymmetric, laminar flow (no azimuthal dependence, $\partial/\partial\phi = 0$)
- Body of revolution with local surface radius $r_0(x)$
- Calorically perfect gas (constant $\gamma$, $c_p$, $\mathrm{Pr}$)
- Temperature-dependent viscosity (Sutherland or power law)
- Thin boundary layer ($\delta \ll r_0$, $\delta \ll L$)

## Coordinates

Body-fitted coordinates are used throughout:

| Symbol | Meaning |
|---|---|
| $x$ | Streamwise arc length along the body surface |
| $y$ | Wall-normal distance from the body surface |
| $r_0(x)$ | Perpendicular distance from the symmetry axis to the body surface |
| $\theta(x)$ | Local surface half-angle |

Within the boundary layer the radial distance from the symmetry axis is

$$
r = r_0(x) + y\cos\theta(x)
$$

Under the thin BL assumption $\delta \ll r_0$, so $r \approx r_0(x)$ throughout
the layer. This approximation is what allows the momentum and energy equations
to retain their 2D form (see below).

## Governing Equations

The compressible laminar boundary layer equations for an axisymmetric body are
(see for example [^schlichting_2017]):

**Continuity**

$$
\frac{\partial(\rho u r_0)}{\partial x}
+ \frac{\partial(\rho v r_0)}{\partial y} = 0
$$

**x-momentum**

$$
\rho\!\left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}\right)
= -\frac{dp}{dx}
+ \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}\right)
$$

**y-momentum**

$$
\frac{\partial p}{\partial y} = 0
$$

**Energy**

$$
\rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
= u\frac{dp}{dx}
+ \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
+ \mu\!\left(\frac{\partial u}{\partial y}\right)^{\!2}
$$

**Perfect gas state equation**

$$
p = \rho R T
$$

From the y-momentum equation, pressure is uniform across the layer and imposed
by the outer inviscid flow.

## Outer Inviscid Flow

Identical to the 2D case. Pressure is set by the edge conditions via the
Euler x-momentum equation:

$$
-\frac{dp}{dx} = \rho_e u_e \frac{du_e}{dx}
$$

Azimuthal symmetry eliminates any swirl contribution.

??? details "Derivation from the Navier-Stokes equations"

    ### Axisymmetric Compressible Navier-Stokes Equations

    Start from the compressible NS for a Newtonian fluid with Stokes' hypothesis
    ($\lambda = -2\mu/3$) in cylindrical coordinates $(x, r)$, retaining
    azimuthal symmetry ($\partial/\partial\phi = 0$, $w = 0$).

    **Continuity**

    $$
    \frac{\partial\rho}{\partial t}
    + \frac{\partial(\rho u)}{\partial x}
    + \frac{1}{r}\frac{\partial(r\rho v)}{\partial r} = 0
    $$

    **x-momentum**

    $$
    \rho\!\left(\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x}
      + v\frac{\partial u}{\partial r}\right)
    = -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\!\left(\frac{\partial v}{\partial r} + \frac{v}{r}\right)\right)
    + \frac{1}{r}\frac{\partial}{\partial r}\!\left(r\mu\!\left(\frac{\partial u}{\partial r}
      + \frac{\partial v}{\partial x}\right)\right)
    $$

    **r-momentum**

    $$
    \rho\!\left(\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x}
      + v\frac{\partial v}{\partial r}\right)
    = -\frac{\partial p}{\partial r}
    + \frac{\partial}{\partial x}\!\left(\mu\!\left(\frac{\partial v}{\partial x}
      + \frac{\partial u}{\partial r}\right)\right)
    + \frac{1}{r}\frac{\partial(r\tau_{rr})}{\partial r}
    - \frac{\tau_{\phi\phi}}{r}
    $$

    where

    $$
    \tau_{rr} = \frac{4\mu}{3}\frac{\partial v}{\partial r}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{v}{r}\right),
    \qquad
    \tau_{\phi\phi} = \frac{4\mu}{3}\frac{v}{r}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial r}\right)
    $$

    **Energy** (in terms of temperature)

    $$
    \rho c_p\!\left(\frac{\partial T}{\partial t} + u\frac{\partial T}{\partial x}
      + v\frac{\partial T}{\partial r}\right)
    = \left(\frac{\partial p}{\partial t} + u\frac{\partial p}{\partial x}
      + v\frac{\partial p}{\partial r}\right)
    + \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{1}{r}\frac{\partial}{\partial r}\!\left(kr\frac{\partial T}{\partial r}\right)
    + \Phi
    $$

    where

    $$
    \Phi = \mu\!\left[
        2\!\left(\frac{\partial u}{\partial x}\right)^{\!2}
        + 2\!\left(\frac{\partial v}{\partial r}\right)^{\!2}
        + 2\!\left(\frac{v}{r}\right)^{\!2}
        + \!\left(\frac{\partial u}{\partial r} + \frac{\partial v}{\partial x}\right)^{\!2}
        - \frac{2}{3}\!\left(\frac{\partial u}{\partial x}
          + \frac{\partial v}{\partial r} + \frac{v}{r}\right)^{\!2}
    \right]
    $$

    ### Steady Flow

    For steady flow, all $\partial/\partial t$ terms vanish.
    The continuity, momentum, and energy equations reduce to:

    $$
    \frac{\partial(\rho u)}{\partial x}
    + \frac{1}{r}\frac{\partial(r\rho v)}{\partial r} = 0
    $$

    $$
    \rho\!\left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial r}\right)
    = -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\!\left(\frac{\partial v}{\partial r} + \frac{v}{r}\right)\right)
    + \frac{1}{r}\frac{\partial}{\partial r}\!\left(r\mu\!\left(\frac{\partial u}{\partial r}
      + \frac{\partial v}{\partial x}\right)\right)
    $$

    $$
    \rho\!\left(u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial r}\right)
    = -\frac{\partial p}{\partial r}
    + \frac{\partial}{\partial x}\!\left(\mu\!\left(\frac{\partial v}{\partial x}
      + \frac{\partial u}{\partial r}\right)\right)
    + \frac{1}{r}\frac{\partial(r\tau_{rr})}{\partial r}
    - \frac{\tau_{\phi\phi}}{r}
    $$

    $$
    \rho c_p\!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial r}\right)
    = u\frac{\partial p}{\partial x} + v\frac{\partial p}{\partial r}
    + \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{1}{r}\frac{\partial}{\partial r}\!\left(kr\frac{\partial T}{\partial r}\right)
    + \Phi
    $$

    ### Boundary Layer Scaling

    The steady NS above simplifies when the boundary layer is thin relative to the
    streamwise length scale. Introduce a reference length $L$ and free-stream
    speed $U_\infty$, and let $\delta \ll L$ be the BL thickness.

    For BL analysis on a body of revolution, it is convenient to work in
    surface-aligned coordinates: $s$ (streamwise arc length along the body) and
    $y$ (wall-normal distance from the surface). The local body surface radius
    $r_0(s)$ can be any smooth function — no specific body geometry is assumed
    at this stage. In surface-aligned coordinates the steady continuity transforms
    to[^schlichting_2017]

    $$
    \frac{\partial(\rho u\, r_0)}{\partial s}
    + \frac{\partial(\rho v\, r_0)}{\partial y} = 0
    $$

    The $r_0(s)$ factor arises from the azimuthal metric of the cylindrical
    geometry; it is $O(L)$ and is **not** removed at any stage. The curvature
    corrections to the momentum and energy equations carry a factor $y/R(s)$
    (where $R$ is the surface radius of curvature) and are $O(\delta/L) = O(\varepsilon)$.
    They will be dropped together with the streamwise viscous terms at Step 5.

    **Step 1: get $V$ scaling from continuity.**
    Both terms must be the same order:

    $$
    \frac{U_\infty}{L} \sim \frac{V}{\delta}
    \qquad\Longrightarrow\qquad
    V \sim \frac{\delta}{L}\,U_\infty
    $$

    **Step 2: get $\delta$ by balancing convection with wall-normal viscosity in x-momentum.**

    $$
    \underbrace{\rho\,\frac{U_\infty^2}{L}}_{\text{inertia}}
    \sim \underbrace{\frac{\mu\,U_\infty}{\delta^2}}_{\text{viscosity}}
    \qquad\Longrightarrow\qquad
    \delta \sim \frac{L}{\sqrt{Re_L}}, \qquad Re_L = \frac{\rho U_\infty L}{\mu}
    $$

    **Step 3: define dimensionless variables.**

    $$
    \begin{aligned}
    s^* &= \frac{s}{L}, &\quad
    y^* &= \frac{y}{\delta}, &\quad
    u^* &= \frac{u}{U_\infty}, \\[6pt]
    v^* &= \frac{v}{\varepsilon U_\infty}, &\quad
    p^* &= \frac{p}{\rho_\infty U_\infty^2}, &\quad
    \rho^* &= \frac{\rho}{\rho_\infty}, \\[6pt]
    \mu^* &= \frac{\mu}{\mu_\infty}, &\quad
    T^* &= \frac{T}{T_\infty}, &\quad
    k^* &= \frac{k}{k_\infty}
    \end{aligned}
    $$

    where $\varepsilon = \delta/L = Re_L^{-1/2} \ll 1$.
    The body radius scales as $r_0^* = r_0/L$.

    **Step 4: Substitute**

    Using $\partial/\partial s = (1/L)\,\partial/\partial s^*$ and
    $\partial/\partial y = (1/\delta)\,\partial/\partial y^*$, substitute into each
    equation and divide out the common dimensional factor.

    **Continuity**

    $$
    \cancel{\frac{\rho_\infty U_\infty}{L}}\,
    \frac{\partial(\rho^* u^* r_0^*)}{\partial s^*}
    + \cancel{\frac{\rho_\infty U_\infty}{L}}\,
    \frac{\partial(\rho^* v^* r_0^*)}{\partial y^*} = 0
    $$

    $$
    \frac{\partial(\rho^* u^* r_0^*)}{\partial s^*}
    + \frac{\partial(\rho^* v^* r_0^*)}{\partial y^*} = 0
    $$

    Note the $r_0^*$ factor is carried through exactly. Both terms are
    $O(1)$ — the body radius does not produce an additional small parameter.

    **x-momentum**

    The cylindrical NS introduces curvature corrections of order $y/R(s)$
    in the inertia and viscous terms. In dimensionless form these carry a
    factor $\varepsilon^2$ (since $y/R \sim \delta/L = \varepsilon$, and the
    leading viscous term already contributes a $1/\varepsilon^2$ factor that
    cancels the $\varepsilon^2$ of the correction — net: $O(\varepsilon^2)$
    relative to the retained term). All other structure is identical to
    the [2D case](2d_equations.md):

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial s^*}
      + v^*\frac{\partial u^*}{\partial y^*}\right)
    &= -\frac{\partial p^*}{\partial s^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right) \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial s^*}\!\left(\frac{4\mu^*}{3}\frac{\partial u^*}{\partial s^*}
          - \frac{2\mu^*}{3}\frac{\partial v^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial s^*}\right)
        + \text{(curvature)}
    \right]
    \end{aligned}
    $$

    **Wall-normal momentum**

    The same order-of-magnitude argument as the 2D y-momentum applies. All inertia
    and viscous terms are $O(\varepsilon)$ or smaller relative to $\partial p^*/\partial y^*$:

    $$
    \varepsilon^2\,\rho^*\!\left(u^*\frac{\partial v^*}{\partial s^*}
      + v^*\frac{\partial v^*}{\partial y^*}\right)
    = -\frac{\partial p^*}{\partial y^*}
    + \varepsilon^2[\cdots]
    $$

    **Energy**

    The azimuthal terms $(v/r)^2$ in $\Phi^*$ and the cylindrical conduction
    term $k\partial T/\partial r / r$ both carry a factor $\varepsilon^2$ in
    dimensionless form (same argument as the curvature terms above). All other
    structure matches the 2D case:

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial s^*}
      + v^*\frac{\partial T^*}{\partial y^*}\right)
    &= \mathrm{Ec}\,u^*\frac{\partial p^*}{\partial s^*}
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2} \\[6pt]
    &+ \varepsilon^2\!\left[\frac{1}{\mathrm{Pr}}\frac{\partial}{\partial s^*}\!\left(k^*\frac{\partial T^*}{\partial s^*}\right)
      + \text{(curvature \& azimuthal terms)}\right]
    \end{aligned}
    $$

    where $\mathrm{Pr} = \mu_\infty c_p / k_\infty$ and $\mathrm{Ec} = U_\infty^2/(c_p T_\infty)$.

    **Step 5: drop $\mathcal{O}(\varepsilon^2)$ terms.**

    Setting $\varepsilon \to 0$:

    $$
    \frac{\partial(\rho^* u^* r_0^*)}{\partial s^*}
    + \frac{\partial(\rho^* v^* r_0^*)}{\partial y^*} = 0
    $$

    $$
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial s^*} + v^*\frac{\partial u^*}{\partial y^*}\right)
    = -\frac{\partial p^*}{\partial s^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right)
    $$

    $$\frac{\partial p^*}{\partial y^*} = 0$$

    $$
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial s^*} + v^*\frac{\partial T^*}{\partial y^*}\right)
    = \mathrm{Ec}\,u^*\frac{\partial p^*}{\partial s^*}
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2}
    $$

    **Step 6: re-dimensionalize.**

    Reversing the substitutions and applying $\partial p^*/\partial y^* = 0
    \Rightarrow \partial p/\partial y = 0$ recovers the dimensional BL equations
    shown at the top of this page. The x-momentum and energy equations are
    identical to the [2D case](2d_equations.md); the only difference is that
    the continuity equation retains the body radius $r_0(s)$.

[^schlichting_2017]: Schlichting, H. & Gersten, K. (2017). *Boundary Layer Theory*, 9th ed. Springer. DOI: [10.1007/978-3-662-52919-5](https://doi.org/10.1007/978-3-662-52919-5)
