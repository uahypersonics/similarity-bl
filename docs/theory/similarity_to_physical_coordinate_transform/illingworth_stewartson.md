# Illingworth-Stewartson

The Illingworth-Stewartson inverse transform maps the similarity coordinate
$\eta$ back to the physical wall-normal coordinate $y$ using a reference
kinematic viscosity.

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

The implemented Illingworth-Stewartson inverse scale is

$$
\eta_{s,IS} = \sqrt{\frac{(m+1)u_e}{2\nu_{ref}x}}
$$

Using $\nu_{ref} = \mu_{ref}/\rho_{ref}$,

$$
\eta_{s,IS}
= \sqrt{\frac{(m+1)u_e\rho_{ref}}{2\mu_{ref}x}}
$$

The inverse coordinate map is

$$
y(\eta) = \frac{1}{\eta_{s,IS}}\int_0^\eta \tau(\hat{\eta})\,d\hat{\eta}
$$