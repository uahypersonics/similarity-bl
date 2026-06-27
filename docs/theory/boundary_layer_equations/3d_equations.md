## Assumptions

- Steady, three-dimensional, laminar flow
- Cartesian coordinates: $x$ streamwise, $z$ crossflow (spanwise), $y$ wall-normal
- Calorically perfect gas (constant $\gamma$, $c_p$, $\mathrm{Pr}$)
- Temperature-dependent viscosity (Sutherland or power law)
- Thin boundary layer ($\delta \ll L$)

## Governing Equations

The three-dimensional compressible laminar boundary layer equations are
(see [^liu_2021]):

**Continuity**

$$
\frac{\partial (\rho u)}{\partial x}
+ \frac{\partial (\rho v)}{\partial y}
+ \frac{\partial (\rho w)}{\partial z} = 0
$$

**x-momentum**

$$
\rho \left( u \frac{\partial u}{\partial x}
           + v \frac{\partial u}{\partial y}
           + w \frac{\partial u}{\partial z} \right)
= -\frac{\partial p}{\partial x}
+ \frac{\partial}{\partial y}\!\left( \mu \frac{\partial u}{\partial y} \right)
$$

**z-momentum**

$$
\rho \left( u \frac{\partial w}{\partial x}
           + v \frac{\partial w}{\partial y}
           + w \frac{\partial w}{\partial z} \right)
= -\frac{\partial p}{\partial z}
+ \frac{\partial}{\partial y}\!\left( \mu \frac{\partial w}{\partial y} \right)
$$

**y-momentum**

$$
\frac{\partial p}{\partial y} = 0
$$

**Energy**

$$
\rho c_p \left( u \frac{\partial T}{\partial x}
               + v \frac{\partial T}{\partial y}
               + w \frac{\partial T}{\partial z} \right)
= \left( u \frac{\partial p}{\partial x} + w \frac{\partial p}{\partial z} \right)
+ \frac{\partial}{\partial y}\!\left( k \frac{\partial T}{\partial y} \right)
+ \mu \left[ \left(\frac{\partial u}{\partial y}\right)^{\!2}
           + \left(\frac{\partial w}{\partial y}\right)^{\!2} \right]
$$

**Perfect gas state equation**

$$
p = \rho R T
$$

From the y-momentum equation it follows that pressure is uniform across the layer and imposed entirely by the outer inviscid flow (see derivation below).

## Outer Inviscid Flow

The boundary layer edge conditions are imposed by the inviscid outer flow:

$$
-\frac{\partial p}{\partial x} = \rho_e u_e \frac{\partial u_e}{\partial x}
+ \rho_e w_e \frac{\partial u_e}{\partial z}
$$

$$
-\frac{\partial p}{\partial z} = \rho_e u_e \frac{\partial w_e}{\partial x}
+ \rho_e w_e \frac{\partial w_e}{\partial z}
$$

??? details "Derivation from the Navier-Stokes equations"

    ### Compressible Navier-Stokes Equations (3D)

    Start with the compressible NS equations for a Newtonian fluid with Stokes' hypothesis ($\lambda = -2\mu/3$).

    **Continuity**

    $$
    \frac{\partial \rho}{\partial t}
    + \frac{\partial (\rho u)}{\partial x}
    + \frac{\partial (\rho v)}{\partial y}
    + \frac{\partial (\rho w)}{\partial z} = 0
    $$

    **x-momentum**

    $$
    \begin{aligned}
    \rho \!\left(\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} + w\frac{\partial u}{\partial z}\right)
    &= -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\!\left(\frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\right)\right) \\
    &+ \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y} + \mu\frac{\partial v}{\partial x}\right)
    + \frac{\partial}{\partial z}\!\left(\mu\frac{\partial u}{\partial z} + \mu\frac{\partial w}{\partial x}\right)
    \end{aligned}
    $$

    **z-momentum**

    $$
    \begin{aligned}
    \rho \!\left(\frac{\partial w}{\partial t} + u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y} + w\frac{\partial w}{\partial z}\right)
    &= -\frac{\partial p}{\partial z}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial w}{\partial x} + \mu\frac{\partial u}{\partial z}\right)
    + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y} + \mu\frac{\partial v}{\partial z}\right) \\
    &+ \frac{\partial}{\partial z}\!\left(\frac{4\mu}{3}\frac{\partial w}{\partial z}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)\right)
    \end{aligned}
    $$

    **y-momentum**

    $$
    \begin{aligned}
    \rho \!\left(\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} + w\frac{\partial v}{\partial z}\right)
    &= -\frac{\partial p}{\partial y}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial v}{\partial x} + \mu\frac{\partial u}{\partial y}\right) \\
    &+ \frac{\partial}{\partial y}\!\left(\frac{4\mu}{3}\frac{\partial v}{\partial y}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial w}{\partial z}\right)\right)
    + \frac{\partial}{\partial z}\!\left(\mu\frac{\partial v}{\partial z} + \mu\frac{\partial w}{\partial y}\right)
    \end{aligned}
    $$

    **Energy** (in terms of temperature)

    $$
    \begin{aligned}
    \rho c_p \!\left(\frac{\partial T}{\partial t} + u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y} + w\frac{\partial T}{\partial z}\right)
    &= \left(\frac{\partial p}{\partial t} + u\frac{\partial p}{\partial x} + v\frac{\partial p}{\partial y} + w\frac{\partial p}{\partial z}\right) \\
    &+ \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
    + \frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right)
    + \Phi
    \end{aligned}
    $$

    where

    $$
    \begin{aligned}
    \Phi = \mu\Bigl[
    &2\!\left(\frac{\partial u}{\partial x}\right)^2
    + 2\!\left(\frac{\partial v}{\partial y}\right)^2
    + 2\!\left(\frac{\partial w}{\partial z}\right)^2 \\
    &+ \!\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2
    + \!\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right)^2
    + \!\left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right)^2 \\
    &- \frac{2}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\right)^2
    \Bigr]
    \end{aligned}
    $$

    ### Steady Flow

    For steady flow, all $\partial/\partial t$ terms vanish.
    The continuity, momentum, and energy equations reduce to:

    $$
    \frac{\partial (\rho u)}{\partial x}
    + \frac{\partial (\rho v)}{\partial y}
    + \frac{\partial (\rho w)}{\partial z} = 0
    $$

    $$
    \begin{aligned}
    \rho \!\left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} + w\frac{\partial u}{\partial z}\right)
    &= -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\!\left(\frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\right)\right) \\
    &+ \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y} + \mu\frac{\partial v}{\partial x}\right)
    + \frac{\partial}{\partial z}\!\left(\mu\frac{\partial u}{\partial z} + \mu\frac{\partial w}{\partial x}\right)
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho \!\left(u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y} + w\frac{\partial w}{\partial z}\right)
    &= -\frac{\partial p}{\partial z}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial w}{\partial x} + \mu\frac{\partial u}{\partial z}\right)
    + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial w}{\partial y} + \mu\frac{\partial v}{\partial z}\right) \\
    &+ \frac{\partial}{\partial z}\!\left(\frac{4\mu}{3}\frac{\partial w}{\partial z}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)\right)
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho \!\left(u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} + w\frac{\partial v}{\partial z}\right)
    &= -\frac{\partial p}{\partial y}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial v}{\partial x} + \mu\frac{\partial u}{\partial y}\right) \\
    &+ \frac{\partial}{\partial y}\!\left(\frac{4\mu}{3}\frac{\partial v}{\partial y}
      - \frac{2\mu}{3}\!\left(\frac{\partial u}{\partial x} + \frac{\partial w}{\partial z}\right)\right)
    + \frac{\partial}{\partial z}\!\left(\mu\frac{\partial v}{\partial z} + \mu\frac{\partial w}{\partial y}\right)
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho c_p \!\left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y} + w\frac{\partial T}{\partial z}\right)
    &= u\frac{\partial p}{\partial x} + v\frac{\partial p}{\partial y} + w\frac{\partial p}{\partial z} \\
    &+ \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
    + \frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right)
    + \Phi
    \end{aligned}
    $$

    ### Boundary Layer Scaling

    The 3D scaling follows the same argument as the 2D case. Both $x$ and $z$ are streamwise-like coordinates ($\sim L$), so $z$-gradients scale identically to $x$-gradients.

    **Step 1: get $V$ scaling from continuity.**

    All three terms must balance:

    $$
    \frac{U_\infty}{L} \sim \frac{V}{\delta} \sim \frac{U_\infty}{L}
    \qquad\Longrightarrow\qquad
    V \sim \frac{\delta}{L}\,U_\infty
    $$

    **Step 2: get $\delta$ by balancing convection with wall-normal viscosity in x-momentum.**

    $$
    \underbrace{\rho\,\frac{U_\infty^2}{L}}_{\text{inertia}} \sim \underbrace{\frac{\mu\,U_\infty}{\delta^2}}_{\text{viscosity}}
    \qquad\Longrightarrow\qquad
    \delta \sim \frac{L}{\sqrt{Re_L}}, \qquad Re_L = \frac{\rho U_\infty L}{\mu}
    $$

    **Step 3: Define dimensionless variables**

    $$
    \begin{aligned}
    x^* &= \frac{x}{L}, &\quad
    y^* &= \frac{y}{\delta}, &\quad
    z^* &= \frac{z}{L}, \\[6pt]
    u^* &= \frac{u}{U_\infty}, &\quad
    v^* &= \frac{v}{\varepsilon U_\infty}, &\quad
    w^* &= \frac{w}{U_\infty}, \\[6pt]
    p^* &= \frac{p}{\rho_\infty U_\infty^2}, &\quad
    \rho^* &= \frac{\rho}{\rho_\infty}, &\quad
    \mu^* &= \frac{\mu}{\mu_\infty}, \\[6pt]
    T^* &= \frac{T}{T_\infty}, &\quad
    k^* &= \frac{k}{k_\infty}
    \end{aligned}
    $$

    where $\varepsilon = \delta/L = Re_L^{-1/2} \ll 1$.

    **Step 4: Substitute**

    Using $\partial/\partial x = (1/L)\,\partial/\partial x^*$, $\partial/\partial z = (1/L)\,\partial/\partial z^*$,
    and $\partial/\partial y = (1/\delta)\,\partial/\partial y^*$, substitute into each equation and divide
    out the common dimensional factor. Terms small by $\varepsilon$ or $\varepsilon^2$ appear with an
    explicit prefactor.

    **Continuity**

    $$
    \cancel{\frac{\rho_\infty U_\infty}{L}}\,\frac{\partial(\rho^* u^*)}{\partial x^*}
    + \cancel{\frac{\rho_\infty U_\infty}{L}}\,\frac{\partial(\rho^* v^*)}{\partial y^*}
    + \cancel{\frac{\rho_\infty U_\infty}{L}}\,\frac{\partial(\rho^* w^*)}{\partial z^*} = 0
    $$

    $$
    \frac{\partial(\rho^* u^*)}{\partial x^*}
    + \frac{\partial(\rho^* v^*)}{\partial y^*}
    + \frac{\partial(\rho^* w^*)}{\partial z^*} = 0
    $$

    All three terms are $\mathcal{O}(1)$ — nothing to drop.

    **x-momentum**

    $$
    \begin{aligned}
    &\cancel{\frac{\rho_\infty U_\infty^2}{L}}\rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*}
      + v^*\frac{\partial u^*}{\partial y^*} + w^*\frac{\partial u^*}{\partial z^*}\right)
    = -\cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial p^*}{\partial x^*}
    + \cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right) \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty U_\infty^2}{L}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(\frac{4\mu^*}{3}\frac{\partial u^*}{\partial x^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial v^*}{\partial y^*} + \frac{\partial w^*}{\partial z^*}\right)\right)
        + \frac{\partial}{\partial z^*}\!\left(\mu^*\frac{\partial u^*}{\partial z^*} + \mu^*\frac{\partial w^*}{\partial x^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*}\right)
    \right]
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*} + v^*\frac{\partial u^*}{\partial y^*} + w^*\frac{\partial u^*}{\partial z^*}\right)
    &= -\frac{\partial p^*}{\partial x^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right) \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial x^*}\!\left(\frac{4\mu^*}{3}\frac{\partial u^*}{\partial x^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial v^*}{\partial y^*} + \frac{\partial w^*}{\partial z^*}\right)\right)
        + \frac{\partial}{\partial z^*}\!\left(\mu^*\frac{\partial u^*}{\partial z^*} + \mu^*\frac{\partial w^*}{\partial x^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*}\right)
    \right]
    \end{aligned}
    $$

    **z-momentum**

    $$
    \begin{aligned}
    &\cancel{\frac{\rho_\infty U_\infty^2}{L}}\rho^*\!\left(u^*\frac{\partial w^*}{\partial x^*}
      + v^*\frac{\partial w^*}{\partial y^*} + w^*\frac{\partial w^*}{\partial z^*}\right)
    = -\cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial p^*}{\partial z^*}
    + \cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial w^*}{\partial y^*}\right) \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty U_\infty^2}{L}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial w^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial z^*}\right)
        + \frac{\partial}{\partial z^*}\!\left(\frac{4\mu^*}{3}\frac{\partial w^*}{\partial z^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*}\right)\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial z^*}\right)
    \right]
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial w^*}{\partial x^*} + v^*\frac{\partial w^*}{\partial y^*} + w^*\frac{\partial w^*}{\partial z^*}\right)
    &= -\frac{\partial p^*}{\partial z^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial w^*}{\partial y^*}\right) \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial w^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial z^*}\right)
        + \frac{\partial}{\partial z^*}\!\left(\frac{4\mu^*}{3}\frac{\partial w^*}{\partial z^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*}\right)\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial z^*}\right)
    \right]
    \end{aligned}
    $$

    **y-momentum**, common factor $\varepsilon\,\rho_\infty U_\infty^2/L$:

    $$
    \begin{gathered}
    \varepsilon\cancel{\frac{\rho_\infty U_\infty^2}{L}}\rho^*\!\left(u^*\frac{\partial v^*}{\partial x^*}
      + v^*\frac{\partial v^*}{\partial y^*} + w^*\frac{\partial v^*}{\partial z^*}\right)
    = -\cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{1}{\varepsilon}\frac{\partial p^*}{\partial y^*} \\[6pt]
    + \varepsilon\cancel{\frac{\rho_\infty U_\infty^2}{L}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\frac{4\mu^*}{3}\frac{\partial v^*}{\partial y^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial w^*}{\partial z^*}\right)\right)
        + \frac{\partial}{\partial z^*}\!\left(\mu^*\frac{\partial v^*}{\partial z^*} + \mu^*\frac{\partial w^*}{\partial y^*}\right)
    \right]
    \end{gathered}
    $$

    Multiplying through by $\varepsilon$:

    $$
    \begin{aligned}
    \varepsilon^2\,\rho^*\!\left(u^*\frac{\partial v^*}{\partial x^*} + v^*\frac{\partial v^*}{\partial y^*} + w^*\frac{\partial v^*}{\partial z^*}\right)
    &= -\frac{\partial p^*}{\partial y^*} \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\frac{4\mu^*}{3}\frac{\partial v^*}{\partial y^*}
          - \frac{2\mu^*}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial w^*}{\partial z^*}\right)\right)
        + \frac{\partial}{\partial z^*}\!\left(\mu^*\frac{\partial v^*}{\partial z^*} + \mu^*\frac{\partial w^*}{\partial y^*}\right)
    \right]
    \end{aligned}
    $$

    At leading order $\partial p^*/\partial y^* = 0$, so pressure is a function of $x^*$ and $z^*$ only.

    **Energy**

    $$
    \begin{aligned}
    &\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*}
      + v^*\frac{\partial T^*}{\partial y^*} + w^*\frac{\partial T^*}{\partial z^*}\right) \\[6pt]
    &\quad= \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\mathrm{Ec}\!\left(u^*\frac{\partial p^*}{\partial x^*}
      + w^*\frac{\partial p^*}{\partial z^*}\right)
    + \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right) \\[6pt]
    &\quad+ \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[
        \!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2}
        + \!\left(\frac{\partial w^*}{\partial y^*}\right)^{\!2}\right] \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{1}{\mathrm{Pr}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(k^*\frac{\partial T^*}{\partial x^*}\right)
        + \frac{\partial}{\partial z^*}\!\left(k^*\frac{\partial T^*}{\partial z^*}\right)\right] \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[
        2\!\left(\frac{\partial u^*}{\partial x^*}\right)^{\!2}
        + 2\!\left(\frac{\partial v^*}{\partial y^*}\right)^{\!2}
        + 2\!\left(\frac{\partial w^*}{\partial z^*}\right)^{\!2}
        + \!\left(\frac{\partial u^*}{\partial z^*} + \frac{\partial w^*}{\partial x^*}\right)^{\!2}
        - \frac{2}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*} + \frac{\partial w^*}{\partial z^*}\right)^{\!2}
    \right]
    \end{aligned}
    $$

    where $\mathrm{Pr} = \mu_\infty c_p / k_\infty$ and $\mathrm{Ec} = U_\infty^2/(c_p T_\infty)$.

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*} + v^*\frac{\partial T^*}{\partial y^*} + w^*\frac{\partial T^*}{\partial z^*}\right)
    &= \mathrm{Ec}\!\left(u^*\frac{\partial p^*}{\partial x^*} + w^*\frac{\partial p^*}{\partial z^*}\right)
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2}
      + \!\left(\frac{\partial w^*}{\partial y^*}\right)^{\!2}\right] \\[6pt]
    &+ \varepsilon^2\frac{1}{\mathrm{Pr}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(k^*\frac{\partial T^*}{\partial x^*}\right)
        + \frac{\partial}{\partial z^*}\!\left(k^*\frac{\partial T^*}{\partial z^*}\right)\right] \\[6pt]
    &+ \varepsilon^2\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[
        2\!\left(\frac{\partial u^*}{\partial x^*}\right)^{\!2}
        + 2\!\left(\frac{\partial v^*}{\partial y^*}\right)^{\!2}
        + 2\!\left(\frac{\partial w^*}{\partial z^*}\right)^{\!2}
        + \!\left(\frac{\partial u^*}{\partial z^*} + \frac{\partial w^*}{\partial x^*}\right)^{\!2}
        - \frac{2}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*} + \frac{\partial w^*}{\partial z^*}\right)^{\!2}
    \right]
    \end{aligned}
    $$

    **Step 5: drop $\mathcal{O}(\varepsilon^2)$ terms.**

    Setting $\varepsilon \to 0$ in each equation above:

    $$
    \frac{\partial(\rho^* u^*)}{\partial x^*}
    + \frac{\partial(\rho^* v^*)}{\partial y^*}
    + \frac{\partial(\rho^* w^*)}{\partial z^*} = 0
    $$

    $$
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*} + v^*\frac{\partial u^*}{\partial y^*} + w^*\frac{\partial u^*}{\partial z^*}\right)
    = -\frac{\partial p^*}{\partial x^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right)
    $$

    $$
    \rho^*\!\left(u^*\frac{\partial w^*}{\partial x^*} + v^*\frac{\partial w^*}{\partial y^*} + w^*\frac{\partial w^*}{\partial z^*}\right)
    = -\frac{\partial p^*}{\partial z^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial w^*}{\partial y^*}\right)
    $$

    $$\frac{\partial p^*}{\partial y^*} = 0$$

    $$
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*} + v^*\frac{\partial T^*}{\partial y^*} + w^*\frac{\partial T^*}{\partial z^*}\right)
    = \mathrm{Ec}\!\left(u^*\frac{\partial p^*}{\partial x^*} + w^*\frac{\partial p^*}{\partial z^*}\right)
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2}
      + \!\left(\frac{\partial w^*}{\partial y^*}\right)^{\!2}\right]
    $$

    **Step 6: re-dimensionalize.**

    Reversing the substitutions and applying $\partial p^*/\partial y^* = 0 \Rightarrow \partial p/\partial y = 0$
    recovers the dimensional BL equations shown at the top of this page.

[^liu_2021]: Liu, Z. (2021). Compressible Falkner–Skan–Cooke boundary layer on a flat plate. *Physics of Fluids*, 33(12). DOI: [10.1063/5.0075233](https://doi.org/10.1063/5.0075233)

