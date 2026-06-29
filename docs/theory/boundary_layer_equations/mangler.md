# Mangler Transformation

The Mangler transformation [^mangler_1948] maps the
[axisymmetric BL equations](axisymmetric_equations.md) to an equivalent 2D
boundary layer, allowing the Falkner-Skan similarity solver to be applied
directly to bodies of revolution such as cones, ogives, and cylinders.

!!! warning
    The Mangler transformation applies to **2D (Falkner-Skan) flows only**.
    Extension to swept axisymmetric bodies (Falkner-Skan-Cooke) is not
    standard and is not covered here.

## Coordinate Transformation

Introduce transformed coordinates

$$
\tilde{x} = \frac{1}{L^2}\int_0^x r_0(\xi)^2\,d\xi,
\qquad
\tilde{y} = \frac{r_0(x)}{L}\,y
$$

where $r_0(x)$ is the local body radius, $y$ is the physical wall-normal
coordinate, and $L$ is an arbitrary reference length. Under this map the
axisymmetric continuity equation takes the standard 2D form in
$(\tilde{x}, \tilde{y})$, so the momentum and energy equations in tilde space
are identical to those of a 2D flat-plate or wedge flow.

## Inverse Transform

After solving the 2D similarity problem in tilde space (obtaining
$\tilde{y}(\eta)$ via the [Levy-Lees](../similarity_to_physical_coordinate_transform/levy_lees.md) inverse), the physical
wall-normal coordinate follows from inverting the second relation:

$$
y(x, \eta) = \frac{L}{r_0(x)}\,\tilde{y}(\eta)
$$

The streamwise direction requires inverting $\tilde{x}(x)$ for the body
geometry of interest; $y$ does not depend on this inversion.

## Composition with the Existing eta2y Transforms

The Mangler map wraps the standard Levy-Lees `eta2y` transform in three steps:

1. **Preprocess** — compute the Mangler-transformed station $\tilde{x}$:

   $$\tilde{x} = \frac{1}{L^2}\int_0^x r_0(\xi)^2\,d\xi$$

2. **Call** `eta2y` with $\tilde{x}$ in place of $x$ to obtain $\tilde{y}(\eta)$.

3. **Postprocess** — scale back to physical $y$:

   $$y = \frac{L}{r_0(x)}\,\tilde{y}$$

The `eta2y` transform itself is unchanged; Mangler enters only as a coordinate
substitution before and a radial rescaling after.

## Common Body Geometries

| Body | $r_0(x)$ | $\tilde{x}$ |
|---|---|---|
| Sharp cone (half-angle $\theta$) | $x\sin\theta$ | $\dfrac{\sin^2\theta}{3L^2}x^3$ |
| Cylinder (radius $R$) | $R$ | $\dfrac{R^2}{L^2}x$ |
| Ogive / arbitrary body | tabulated $r_0(x)$ | numerical quadrature |

For a **sharp cone** the Mangler $\tilde{x}$ grows as $x^3$, which compresses
the effective streamwise coordinate and produces a thinner equivalent 2D BL
relative to the physical arc length.

## Reference Length Choice

$L$ cancels in the final physical coordinate $y$ because $\tilde{y} \propto
1/L$ while the postprocessing scale is $L/r_0$. The choice of $L$ therefore
affects only the numerical magnitude of $\tilde{x}$ and $\tilde{y}$
individually, not the computed $y$. A convenient default is $L = 1\,\text{m}$
(or $L = r_\text{nose}$ for ogive bodies) to keep $\tilde{x}$ of order unity.

[^mangler_1948]: Mangler, W. (1948). Zusammenhang zwischen ebenen und rotationssymmetrischen Grenzschichten in kompressiblen Flüssigkeiten. *Zeitschrift für Angewandte Mathematik und Mechanik*, 28(4), 97–103. <https://doi.org/10.1002/zamm.19480280401>
