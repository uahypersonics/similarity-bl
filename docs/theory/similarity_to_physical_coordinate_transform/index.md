# Similarity to Physical Coordinate Transform

The Falkner-Skan and Falkner-Skan-Cooke solvers return profiles on $\eta$. This transform maps those profiles to physical $y$.

## Conventions

The similarity solutions use the nondimensional temperature profile

$$
\tau(\eta) = \frac{T(\eta)}{T_e}
$$

where $T_e$ is the boundary-layer edge temperature at the streamwise station of
interest. For a perfect gas boundary layer with $\partial p/\partial y = 0$,
the density relation is

$$
\rho = \frac{\rho_e}{\tau}
$$

(For the derivation see [Levy-Lees](../falkner_skan/derivation.md) and
[Illingworth-Stewartson](../falkner_skan/derivation_is.md).)

The edge velocity is represented by the Falkner-Skan power law

$$
u_e(x) = Cx^m
$$

and the Hartree pressure-gradient parameter is

$$
\beta_H = \frac{2m}{m+1}
$$

## General Inverse Form

Both implemented transforms are density-weighted wall-normal coordinates. After
differentiating the chosen definition of $\eta$ with respect to $y$, the density
factor is replaced using $\rho = \rho_e/\tau$. The remaining factors depend only
on the local station and the selected transform, so they are collected into a
single inverse length scale $\eta_s$.

With that notation, both implemented transforms can be written locally as

$$
\frac{\partial \eta}{\partial y} = \frac{\eta_s}{\tau}
$$

The scale $\eta_s$ is therefore not an additional assumption. It is just the
coefficient multiplying $1/\tau$ after the coordinate transform has been
differentiated. Its exact definition is different for the Levy-Lees and
Illingworth-Stewartson transforms.

??? note "Levy-Lees: η definition and η_s"

    $$
    \eta = \sqrt{\frac{(m+1)\rho_e u_e}{2\mu_e x}}\int_0^y \frac{\rho}{\rho_e}\,dy'
    $$

    Differentiating with respect to $y$ and substituting $\rho/\rho_e = 1/\tau$:

    $$
    \frac{\partial\eta}{\partial y} = \frac{1}{\tau}\underbrace{\sqrt{\frac{(m+1)\rho_e u_e}{2\mu_e x}}}_{\eta_{s,LL}}
    $$

??? note "Illingworth-Stewartson: η definition and η_s"

    $$
    \eta = \sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}\int_0^y \frac{\rho}{\rho_e}\,dy'
    $$

    where $\nu_{ref} = \mu_{ref}/\rho_{ref}$ is a prescribed reference kinematic viscosity
    (defaults to the local edge value $\nu_e = \mu_e/\rho_e$).

    Differentiating with respect to $y$ and substituting $\rho/\rho_e = 1/\tau$:

    $$
    \frac{\partial\eta}{\partial y} = \frac{1}{\tau}\underbrace{\sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}}_{\eta_{s,IS}}
    $$

    With the default reference state ($\nu_{ref} = \nu_e$), $\eta_{s,IS} = \eta_{s,LL}$.

Inverting gives

$$
\frac{\partial y}{\partial \eta} = \frac{\tau}{\eta_s}
$$

Integrating from the wall, where $\eta = 0$ and $y = 0$,

$$
y(\eta) = \frac{1}{\eta_s}\int_0^\eta \tau(\hat{\eta})\,d\hat{\eta}
$$

Thus the transform has two pieces:

1. Compute the scale $\eta_s$ for the chosen similarity transformation.
2. Integrate $\tau(\eta)$ over the similarity grid.

## Available Transforms

| Transform | Use |
|---|---|
| [Levy-Lees](levy_lees.md) | Default inverse map for the `falkner_skan` equation family |
| [Illingworth-Stewartson](illingworth_stewartson.md) | Default inverse map for the `falkner_skan_cooke` equation family |
