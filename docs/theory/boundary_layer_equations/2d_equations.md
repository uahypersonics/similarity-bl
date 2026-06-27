## Assumptions

- Steady, two-dimensional, laminar flow
- Calorically perfect gas (constant $\gamma$, $c_p$, $\mathrm{Pr}$)
- Temperature-dependent viscosity (Sutherland or power law)
- Thin boundary layer ($\delta \ll L$)

## Governing Equations

The compressible laminar boundary layer equations for a calorically perfect gas
are (see for example [^white_2006][^schlichting_2017]):

**Continuity**

$$
\frac{\partial (\rho u)}{\partial x} + \frac{\partial (\rho v)}{\partial y} = 0
$$

**x-momentum**

$$
\rho \left( u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} \right)
= -\frac{dp}{dx} + \frac{\partial}{\partial y}\!\left( \mu \frac{\partial u}{\partial y} \right)
$$

**y-momentum**

$$
\frac{\partial p}{\partial y} = 0
$$

**Energy**

$$
\rho c_p \left( u \frac{\partial T}{\partial x} + v \frac{\partial T}{\partial y} \right)
= u \frac{dp}{dx}
  + \frac{\partial}{\partial y}\!\left( k \frac{\partial T}{\partial y} \right)
  + \mu \left(\frac{\partial u}{\partial y}\right)^{\!2}
$$

**Perfect gas state equation**

$$
p = \rho R T
$$

From the y-momentum equation it follows that pressure is uniform across the layer and imposed entirely by the outer inviscid flow (see derivation below).

## Outer Inviscid Flow

The pressure gradient is set by the inviscid outer flow via the Euler momentum equation:

$$
-\frac{dp}{dx} = \rho_e u_e \frac{du_e}{dx}
$$

??? details "Derivation from the Navier-Stokes equations"

    ### Compressible Navier-Stokes Equations (2D)

    To obtain the boundary layer equations, start with the compressible Navier-Stokes (NS) equations for a Newtonian fluid with Stokes' hypothesis ($\lambda = -2\mu/3$).

    **Continuity**

    $$
    \frac{\partial \rho}{\partial t}
    + \frac{\partial (\rho u)}{\partial x}
    + \frac{\partial (\rho v)}{\partial y} = 0
    $$

    **x-momentum**

    $$
    \rho \left(\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}\right)
    = -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\frac{\partial v}{\partial y}\right)
    + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}
      + \mu\frac{\partial v}{\partial x}\right)
    $$

    **y-momentum**

    $$
    \rho \left(\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y}\right)
    = -\frac{\partial p}{\partial y}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial v}{\partial x}
      + \mu\frac{\partial u}{\partial y}\right)
    + \frac{\partial}{\partial y}\!\left(\frac{4\mu}{3}\frac{\partial v}{\partial y}
      - \frac{2\mu}{3}\frac{\partial u}{\partial x}\right)
    $$

    **Energy** (in terms of temperature)

    $$
    \rho c_p \left(\frac{\partial T}{\partial t} + u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
    = \left(\frac{\partial p}{\partial t} + u\frac{\partial p}{\partial x} + v\frac{\partial p}{\partial y}\right)
    + \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
    + \Phi
    $$

    where

    $$
    \Phi = \mu\left[2\left(\frac{\partial u}{\partial x}\right)^2 + 2\left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2 - \frac{2}{3}\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)^2\right]
    $$

    ### Steady Flow

    For steady flow, all $\partial/\partial t$ terms vanish.
    The continuity, momentum, and energy equations reduce to:

    $$
    \frac{\partial (\rho u)}{\partial x} + \frac{\partial (\rho v)}{\partial y} = 0
    $$

    $$
    \rho \left(u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}\right)
    = -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial x}\!\left(\frac{4\mu}{3}\frac{\partial u}{\partial x}
      - \frac{2\mu}{3}\frac{\partial v}{\partial y}\right)
    + \frac{\partial}{\partial y}\!\left(\mu\frac{\partial u}{\partial y}
      + \mu\frac{\partial v}{\partial x}\right)
    $$

    $$
    \rho \left(u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y}\right)
    = -\frac{\partial p}{\partial y}
    + \frac{\partial}{\partial x}\!\left(\mu\frac{\partial v}{\partial x}
      + \mu\frac{\partial u}{\partial y}\right)
    + \frac{\partial}{\partial y}\!\left(\frac{4\mu}{3}\frac{\partial v}{\partial y}
      - \frac{2\mu}{3}\frac{\partial u}{\partial x}\right)
    $$

    $$
    \rho c_p \left(u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}\right)
    = u\frac{\partial p}{\partial x} + v\frac{\partial p}{\partial y}
    + \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
    + \Phi
    $$

    ### Boundary Layer Scaling

    The steady NS equations above simplify when the boundary layer is thin relative to the streamwise length scale.
    Introduce a reference length $L$ and free-stream speed $U_\infty$, and let $\delta \ll L$ be the BL
    thickness.

    **Step 1: get $V$ scaling from continuity.**
    Both terms of the continuity equation must be the same order:

    $$
    \frac{U_\infty}{L} \sim \frac{V}{\delta}
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

    **Step 4: Substitute**

    Using $\partial/\partial x = (1/L)\,\partial/\partial x^*$ and $\partial/\partial y = (1/\delta)\,\partial/\partial y^*$,
    substitute into each equation and divide out the common dimensional factor. Terms small by
    $\varepsilon$ or $\varepsilon^2$ appear with an explicit prefactor.

    **Continuity**

    $$
    \cancel{\frac{\rho_\infty U_\infty}{L}}\,\frac{\partial(\rho^* u^*)}{\partial x^*}
    + \cancel{\frac{\rho_\infty U_\infty}{L}}\,\frac{\partial(\rho^* v^*)}{\partial y^*} = 0
    $$

    $$
    \frac{\partial(\rho^* u^*)}{\partial x^*} + \frac{\partial(\rho^* v^*)}{\partial y^*} = 0
    $$

    **x-momentum**

    $$
    \begin{aligned}
    &\cancel{\frac{\rho_\infty U_\infty^2}{L}}\rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*} + v^*\frac{\partial u^*}{\partial y^*}\right)
    = -\cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial p^*}{\partial x^*}
    + \cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right) \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty U_\infty^2}{L}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(\frac{4\mu^*}{3}\frac{\partial u^*}{\partial x^*}
        - \frac{2\mu^*}{3}\frac{\partial v^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*}\right)
    \right]
    \end{aligned}
    $$

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*} + v^*\frac{\partial u^*}{\partial y^*}\right)
    &= -\frac{\partial p^*}{\partial x^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right) \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial x^*}\!\left(\frac{4\mu^*}{3}\frac{\partial u^*}{\partial x^*}
        - \frac{2\mu^*}{3}\frac{\partial v^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*}\right)
    \right]
    \end{aligned}
    $$

    **y-momentum**

    $$
    \begin{gathered}
    \varepsilon\cancel{\frac{\rho_\infty U_\infty^2}{L}}\rho^*\!\left(u^*\frac{\partial v^*}{\partial x^*} + v^*\frac{\partial v^*}{\partial y^*}\right)
    = -\cancel{\frac{\rho_\infty U_\infty^2}{L}}\frac{1}{\varepsilon}\frac{\partial p^*}{\partial y^*} \\[6pt]
    + \varepsilon\cancel{\frac{\rho_\infty U_\infty^2}{L}}\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\frac{4\mu^*}{3}\frac{\partial v^*}{\partial y^*}
        - \frac{2\mu^*}{3}\frac{\partial u^*}{\partial x^*}\right)
    \right]
    \end{gathered}
    $$

    Multiplying through by $\varepsilon$:

    $$
    \begin{aligned}
    \varepsilon^2\,\rho^*\!\left(u^*\frac{\partial v^*}{\partial x^*} + v^*\frac{\partial v^*}{\partial y^*}\right)
    &= -\frac{\partial p^*}{\partial y^*} \\[6pt]
    &+ \varepsilon^2\!\left[
        \frac{\partial}{\partial x^*}\!\left(\mu^*\frac{\partial v^*}{\partial x^*} + \mu^*\frac{\partial u^*}{\partial y^*}\right)
        + \frac{\partial}{\partial y^*}\!\left(\frac{4\mu^*}{3}\frac{\partial v^*}{\partial y^*}
        - \frac{2\mu^*}{3}\frac{\partial u^*}{\partial x^*}\right)
    \right]
    \end{aligned}
    $$

    **Energy**

    $$
    \begin{aligned}
    &\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*} + v^*\frac{\partial T^*}{\partial y^*}\right) \\[6pt]
    &\quad= \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\mathrm{Ec}\,u^*\frac{\partial p^*}{\partial x^*}
    + \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2} \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{1}{\mathrm{Pr}}\frac{\partial}{\partial x^*}\!\left(k^*\frac{\partial T^*}{\partial x^*}\right) \\[6pt]
    &\quad+ \varepsilon^2\cancel{\frac{\rho_\infty c_p T_\infty U_\infty}{L}}\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[
            2\!\left(\frac{\partial u^*}{\partial x^*}\right)^{\!2}
            + 2\!\left(\frac{\partial v^*}{\partial y^*}\right)^{\!2}
            + \!\left(\frac{\partial v^*}{\partial x^*}\right)^{\!2}
            - \frac{2}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*}\right)^{\!2}
        \right]
    \end{aligned}
    $$

    where $\mathrm{Pr} = \mu_\infty c_p / k_\infty$ and $\mathrm{Ec} = U_\infty^2/(c_p T_\infty)$.

    $$
    \begin{aligned}
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*} + v^*\frac{\partial T^*}{\partial y^*}\right)
    &= \mathrm{Ec}\,u^*\frac{\partial p^*}{\partial x^*}
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2} \\[6pt]
    &+ \varepsilon^2\frac{1}{\mathrm{Pr}}\frac{\partial}{\partial x^*}\!\left(k^*\frac{\partial T^*}{\partial x^*}\right) \\[6pt]
    &+ \varepsilon^2\frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left[
        2\!\left(\frac{\partial u^*}{\partial x^*}\right)^{\!2}
        + 2\!\left(\frac{\partial v^*}{\partial y^*}\right)^{\!2}
        + \!\left(\frac{\partial v^*}{\partial x^*}\right)^{\!2}
        - \frac{2}{3}\!\left(\frac{\partial u^*}{\partial x^*} + \frac{\partial v^*}{\partial y^*}\right)^{\!2}
    \right]
    \end{aligned}
    $$

    **Step 5: drop $\mathcal{O}(\varepsilon^2)$ terms.**

    Setting $\varepsilon \to 0$ in each equation above:

    $$
    \frac{\partial(\rho^* u^*)}{\partial x^*} + \frac{\partial(\rho^* v^*)}{\partial y^*} = 0
    $$

    $$
    \rho^*\!\left(u^*\frac{\partial u^*}{\partial x^*} + v^*\frac{\partial u^*}{\partial y^*}\right)
    = -\frac{\partial p^*}{\partial x^*}
    + \frac{\partial}{\partial y^*}\!\left(\mu^*\frac{\partial u^*}{\partial y^*}\right)
    $$

    $$\frac{\partial p^*}{\partial y^*} = 0$$

    $$
    \rho^*\!\left(u^*\frac{\partial T^*}{\partial x^*} + v^*\frac{\partial T^*}{\partial y^*}\right)
    = \mathrm{Ec}\,u^*\frac{\partial p^*}{\partial x^*}
    + \frac{1}{\mathrm{Pr}}\frac{\partial}{\partial y^*}\!\left(k^*\frac{\partial T^*}{\partial y^*}\right)
    + \frac{\mathrm{Ec}}{\mathrm{Pr}}\,\mu^*\!\left(\frac{\partial u^*}{\partial y^*}\right)^{\!2}
    $$

    **Step 6: re-dimensionalize.**

    Reversing the substitutions ($u^* \to u/U_\infty$, $y^* \to y/\delta$, $p^* \to p/(\rho_\infty U_\infty^2)$, etc.)
    and applying $\partial p^*/\partial y^* = 0 \Rightarrow \partial p/\partial y = 0$ recovers the
    dimensional BL equations shown at the top of this page.

[^white_2006]: White, F. M. (2006). *Viscous Fluid Flow*, 3rd ed. McGraw-Hill, New York.
[^schlichting_2017]: Schlichting, H. & Gersten, K. (2017). *Boundary Layer Theory*, 9th ed. Springer. DOI: [10.1007/978-3-662-52919-5](https://doi.org/10.1007/978-3-662-52919-5)


