# Illingworth-Stewartson

The Illingworth-Stewartson inverse transform maps the similarity coordinate
$\eta$ back to the physical wall-normal coordinate $y$ using a reference
kinematic viscosity.

## Reference State

The reference kinematic viscosity is

$$
\nu_{ref} = \frac{\mu_{ref}}{\rho_{ref}}
$$

where $\rho_{ref}$ and $\mu_{ref}$ may be specified independently. In local
self-similar usage, `simbl` defaults these reference values to the local edge
state:

$$
\rho_{ref} = \rho_e, \qquad \mu_{ref} = \mu_e
$$

## IS Coordinate Definition

The Illingworth-Stewartson coordinate generalizes the Levy-Lees coordinate by
replacing the local edge kinematic viscosity $\nu_e$ with a prescribed
$\nu_{ref}$:

$$
\eta = \sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}\int_0^y\frac{\rho}{\rho_e}\,dy'
$$

This form follows from the IS stream-function variable and power-law edge
velocity $u_e = Cx^m$, consistent with the IS transformation used in the
[Falkner-Skan-Cooke derivation](../falkner_skan_cooke/derivation_is.md).

## Derivation of the Inverse Scale

Differentiate the IS coordinate with respect to $y$ at fixed $x$:

$$
\frac{\partial\eta}{\partial y}
= \sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}\cdot\frac{\rho}{\rho_e}
$$

Substitute $\rho/\rho_e = 1/\tau$ from the ideal gas law (see
[index](index.md)):

$$
\frac{\partial\eta}{\partial y}
= \frac{1}{\tau}\sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}
$$

Therefore the Illingworth-Stewartson inverse scale is

$$
\eta_{s,IS} = \sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}
$$

Inverting gives

$$
\frac{\partial y}{\partial\eta} = \frac{\tau}{\eta_{s,IS}}
$$

and integrating from the wall ($\eta=0$, $y=0$),

$$
y(\eta) = \frac{1}{\eta_{s,IS}}\int_0^\eta \tau(\hat{\eta})\,d\hat{\eta}
$$

## Reduction to Levy-Lees with Default Reference State

When the reference state defaults to the local edge state
($\rho_{ref} = \rho_e$, $\mu_{ref} = \mu_e$), the reference kinematic
viscosity is $\nu_{ref} = \nu_e = \mu_e/\rho_e$, and the IS scale becomes

$$
\eta_{s,IS}
= \sqrt{\frac{(m+1)u_e}{2\nu_e x}}
= \sqrt{\frac{(m+1)\rho_e u_e}{2\mu_e x}}
= \eta_{s,LL}
$$

which is exactly the [Levy-Lees scale](levy_lees.md). With the default
reference state, both transforms produce identical physical coordinates for any
$\tau$ profile. The IS variant is only distinct when an explicit non-local
reference state is supplied (for example, freestream values
$\rho_\infty$, $\mu_\infty$), which decouples the scale from local variations
in edge conditions along the body.