# Quasi-3D Boundary Layer

!!! note
    These equations are obtained from the
    [3D boundary layer equations](3d_equations.md) by assuming homogeneous flow in
    the spanwise direction ($\partial/\partial z = 0$). They are the starting point for the
    [Falkner-Skan-Cooke derivation](../falkner_skan_cooke/derivation.md).

## Assumptions

!!! abstract "Inherited from the 3D boundary layer equations"
    - Steady, three-dimensional, laminar flow
    - Cartesian coordinates: $x$ streamwise, $z$ crossflow (spanwise), $y$ wall-normal
    - Calorically perfect gas (constant $\gamma$, $c_p$, $\mathrm{Pr}$)
    - Temperature-dependent viscosity (Sutherland or power law)
    - Thin boundary layer ($\delta \ll L$)

**Additional assumption:**

$$\frac{\partial}{\partial z} = 0$$

The flow has three velocity components $(u, v, w)$ but no gradients in the spanwise
$z$-direction.

## Reduction

Starting from the [3D boundary layer equations](3d_equations.md), applying
$\partial/\partial z = 0$ term by term:

**Continuity.** The $\partial(\rho w)/\partial z$ term drops:

$$
\frac{\partial (\rho u)}{\partial x}
+ \frac{\partial (\rho v)}{\partial y}
+ \cancel{\frac{\partial (\rho w)}{\partial z}} = 0
$$

**x-momentum.** The $w\,\partial u/\partial z$ convective term drops:

$$
\rho \left( u \frac{\partial u}{\partial x}
           + v \frac{\partial u}{\partial y}
           + \cancel{w \frac{\partial u}{\partial z}} \right)
= -\frac{\partial p}{\partial x}
+ \frac{\partial}{\partial y}\!\left( \mu \frac{\partial u}{\partial y} \right)
$$

**z-momentum.** The $w\,\partial w/\partial z$ convective term drops. Assuming homogeneous flow in the spanwise direction
also implies no spanwise pressure gradient ($\partial p/\partial z = 0$), so the pressure
driving term drops too:

$$
\rho \left( u \frac{\partial w}{\partial x}
           + v \frac{\partial w}{\partial y}
           + \cancel{w \frac{\partial w}{\partial z}} \right)
= \cancel{-\frac{\partial p}{\partial z}}
+ \frac{\partial}{\partial y}\!\left( \mu \frac{\partial w}{\partial y} \right)
$$

**y-momentum.** Unchanged:

$$\frac{\partial p}{\partial y} = 0$$

**Energy.** The $w\,\partial T/\partial z$ convective term and the $w\,\partial p/\partial z$
pressure work term both drop:

$$
\begin{aligned}
\rho c_p \left( u \frac{\partial T}{\partial x}
               + v \frac{\partial T}{\partial y}
               + \cancel{w \frac{\partial T}{\partial z}} \right)
&= \left( u \frac{\partial p}{\partial x} + \cancel{w \frac{\partial p}{\partial z}} \right)
+ \frac{\partial}{\partial y}\!\left( k \frac{\partial T}{\partial y} \right) \\
&+ \mu \left[ \left(\frac{\partial u}{\partial y}\right)^{\!2}
           + \left(\frac{\partial w}{\partial y}\right)^{\!2} \right]
\end{aligned}
$$

## Governing Equations

!!! info ""

    **Continuity**

    $$
    \frac{\partial (\rho u)}{\partial x}
    + \frac{\partial (\rho v)}{\partial y} = 0
    $$

    **x-momentum**

    $$
    \rho \left( u \frac{\partial u}{\partial x}
               + v \frac{\partial u}{\partial y} \right)
    = -\frac{\partial p}{\partial x}
    + \frac{\partial}{\partial y}\!\left( \mu \frac{\partial u}{\partial y} \right)
    $$

    **z-momentum**

    $$
    \rho \left( u \frac{\partial w}{\partial x}
               + v \frac{\partial w}{\partial y} \right)
    = \frac{\partial}{\partial y}\!\left( \mu \frac{\partial w}{\partial y} \right)
    $$

    **y-momentum**

    $$\frac{\partial p}{\partial y} = 0$$

    **Energy**

    $$
    \rho c_p \left( u \frac{\partial T}{\partial x}
                   + v \frac{\partial T}{\partial y} \right)
    = u \frac{\partial p}{\partial x}
    + \frac{\partial}{\partial y}\!\left( k \frac{\partial T}{\partial y} \right)
    + \mu \left[ \left(\frac{\partial u}{\partial y}\right)^{\!2}
               + \left(\frac{\partial w}{\partial y}\right)^{\!2} \right]
    $$

## Outer Inviscid Flow

Applying $\partial/\partial z = 0$ to the 3D outer flow relations:

$$-\frac{\partial p}{\partial x} = \rho_e u_e \frac{\partial u_e}{\partial x}$$

The spanwise pressure gradient $-\partial p/\partial z$ vanishes, consistent with
the assumption of homogeneous flow in the spanwise direction. The edge crossflow velocity $w_e$ may still be nonzero (it
is prescribed as a boundary condition).

## Relation to 2D System

- Setting $w = 0$ (and $w_e = 0$) recovers the
  [2D boundary layer equations](2d_equations.md)

