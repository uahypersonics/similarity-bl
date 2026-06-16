# Levy-Lees

The Levy-Lees inverse transform maps the Falkner-Skan similarity coordinate
$\eta$ back to the physical wall-normal coordinate $y$.

The Levy-Lees coordinate used in the Falkner-Skan derivation is

$$
\eta = \frac{u_e}{\sqrt{2\xi}}\int_0^y \rho\,dy'
$$

with

$$
\xi = \int_0^x \rho_e \mu_e u_e\,dx'
$$

For a power-law edge velocity $u_e = Cx^m$, and with $\rho_e\mu_e$ constant over
the local similarity station,

$$
\xi = \rho_e\mu_e\int_0^x Cx'^m\,dx'
$$

$$
\xi = \rho_e\mu_e C\frac{x^{m+1}}{m+1}
$$

Using $u_e = Cx^m$,

$$
\xi = \frac{\rho_e\mu_e u_e x}{m+1}
$$

Differentiate the Levy-Lees coordinate with respect to $y$ at fixed $x$:

$$
\frac{\partial \eta}{\partial y}
= \frac{u_e}{\sqrt{2\xi}}\rho
$$

Substitute $\rho = \rho_e/\tau$:

$$
\frac{\partial \eta}{\partial y}
= \frac{\rho_e u_e}{\tau\sqrt{2\xi}}
$$

Now substitute the power-law expression for $\xi$:

$$
\frac{\partial \eta}{\partial y}
= \frac{1}{\tau}\frac{\rho_e u_e}{\sqrt{2\rho_e\mu_e u_e x/(m+1)}}
$$

Collecting terms gives

$$
\frac{\partial \eta}{\partial y}
= \frac{1}{\tau}\sqrt{\frac{(m+1)\rho_e u_e}{2\mu_e x}}
$$

Therefore the Levy-Lees inverse scale is

$$
\eta_{s,LL} = \sqrt{\frac{(m+1)\rho_e u_e}{2\mu_e x}}
$$

and the physical coordinate is

$$
y(\eta) = \frac{1}{\eta_{s,LL}}\int_0^\eta \tau(\hat{\eta})\,d\hat{\eta}
$$