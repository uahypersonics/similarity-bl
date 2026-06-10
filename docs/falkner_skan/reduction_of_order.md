# Reduction of Order

Starting from the [derived ODEs](derivation.md):

$$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

$$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau'
+ (\gamma-1)M_e^2\!\left[Cf''^2 - \beta_H\tau f'\right] = 0$$

where $C = \rho\mu/(\rho_e\mu_e)$ is the Chapman–Rubesin factor, $f' = u/u_e$ the
dimensionless velocity, $\tau = T/T_e$ the dimensionless temperature, and all primes
denote differentiation with respect to the wall-normal similarity variable $\eta$.

## Chapman–Rubesin factor

For a perfect gas $\rho/\rho_e = 1/\tau$, so the Chapman–Rubesin factor simplifies to:

$$C = \frac{\rho\mu}{\rho_e\mu_e} = \frac{\mu}{\mu_e\tau}$$

Differentiating with respect to $\eta$ via the quotient rule:

$$C' = \frac{d}{d\eta}\!\left(\frac{\mu}{\mu_e\tau}\right)
= \frac{1}{\mu_e}\frac{\mu'\tau - \mu\tau'}{\tau^2}$$

## Expanding the compound derivatives

Applying the product rule to $(Cf'')'$ and $\bigl(C\tau'/\mathrm{Pr}\bigr)'$:

$$(Cf'')' = C'f'' + Cf'''$$

$$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' = \frac{C'}{\mathrm{Pr}}\tau' + \frac{C}{\mathrm{Pr}}\tau''$$

Substituting into the two ODEs:

$$C'f'' + Cf''' + ff'' + \beta_H(\tau - f'^2) = 0$$

$$\frac{C'}{\mathrm{Pr}}\tau' + \frac{C}{\mathrm{Pr}}\tau'' + f\tau'
+ (\gamma-1)M_e^2\!\left[Cf''^2 - \beta_H\tau f'\right] = 0$$

## Substituting C and C'

Replacing $C = \mu/(\mu_e\tau)$ and $C' = (\mu'\tau - \mu\tau')/(\mu_e\tau^2)$:

**Momentum:**

$$\frac{\mu'\tau - \mu\tau'}{\mu_e\tau^2}\,f''
+ \frac{\mu}{\mu_e\tau}\,f'''
+ ff'' + \beta_H(\tau - f'^2) = 0$$

**Energy:**

$$\frac{(\mu'\tau - \mu\tau')}{\mathrm{Pr}\,\mu_e\tau^2}\,\tau'
+ \frac{\mu}{\mathrm{Pr}\,\mu_e\tau}\,\tau''
+ f\tau'
+ (\gamma-1)M_e^2\!\left[\frac{\mu}{\mu_e\tau}\,f''^2 - \beta_H\tau f'\right] = 0$$

## Isolating the highest derivatives

### Momentum

Multiply by $\mu_e$:

$$\frac{\mu'\tau - \mu\tau'}{\tau^2}\,f''
+ \frac{\mu}{\tau}\,f'''
+ \mu_e ff'' + \mu_e\beta_H(\tau - f'^2) = 0$$

Multiply by $\tau/\mu$:

$$\frac{\mu'\tau - \mu\tau'}{\mu\tau}\,f''
+ f'''
+ \frac{\mu_e\tau}{\mu}\,ff''
+ \frac{\mu_e\tau}{\mu}\,\beta_H(\tau - f'^2) = 0$$

Noting that $\mu_e\tau/\mu = 1/C$ and $(\mu'\tau - \mu\tau')/(\mu\tau) = \mu'/\mu - \tau'/\tau$, and rearranging:

$$\boxed{f''' = \frac{\tau'}{\tau}\,f'' - \frac{\mu'}{\mu}\,f'' - \frac{ff''}{C} - \frac{\beta_H(\tau - f'^2)}{C}}$$

### Energy

Multiply by $\mu_e$:

$$\frac{(\mu'\tau - \mu\tau')}{\mathrm{Pr}\,\tau^2}\,\tau'
+ \frac{\mu}{\mathrm{Pr}\,\tau}\,\tau''
+ \mu_e f\tau'
+ (\gamma-1)M_e^2\!\left[\frac{\mu}{\tau}\,f''^2 - \beta_H\mu_e\tau f'\right] = 0$$

Multiply by $\mathrm{Pr}\,\tau/\mu$, using $\mu_e\tau/\mu = 1/C$ and $\mu_e\tau^2/\mu = \tau/C$:

$$\frac{(\mu'\tau - \mu\tau')}{\mu\tau}\,\tau'
+ \tau''
+ \frac{\mathrm{Pr}\,f\tau'}{C}
+ \mathrm{Pr}(\gamma-1)M_e^2\!\left[f''^2 - \frac{\beta_H\tau f'}{C}\right] = 0$$

Applying $(\mu'\tau - \mu\tau')/(\mu\tau) = \mu'/\mu - \tau'/\tau$ and rearranging:

$$\boxed{\tau'' = \frac{\tau'^2}{\tau} - \frac{\mu'\tau'}{\mu}
- \frac{\mathrm{Pr}\,f\tau'}{C}
- \mathrm{Pr}(\gamma-1)M_e^2\!\left[f''^2 - \frac{\beta_H\tau f'}{C}\right]}$$

## Viscosity derivatives

The $\mu'/\mu$ term requires the derivative of viscosity with respect to $\eta$.
By the chain rule:

$$\mu' = \frac{d\mu}{d\eta} = \frac{d\mu}{dT}\frac{dT}{d\eta}$$

Since $T = T_e\tau(\eta)$:

$$\frac{dT}{d\eta} = T_e\frac{d\tau}{d\eta} = T_e\tau'$$

Therefore:

$$\mu' = \frac{d\mu}{dT}\,T_e\tau'$$

Substituting into the boxed expressions:

$$\boxed{f''' = \frac{\tau'}{\tau}\,f''
- \frac{T_e}{\mu}\frac{d\mu}{dT}\,\tau'f''
- \frac{ff''}{C}
- \frac{\beta_H(\tau - f'^2)}{C}}$$

$$\boxed{\tau'' = \frac{\tau'^2}{\tau}
- \frac{T_e}{\mu}\frac{d\mu}{dT}\,\tau'^2
- \frac{\mathrm{Pr}\,f\tau'}{C}
- \mathrm{Pr}(\gamma-1)M_e^2\!\left[f''^2 - \frac{\beta_H\tau f'}{C}\right]}$$

where $T = T_e\tau$ and $C = \mu(T)/(\mu_e\tau)$ are evaluated locally at each $\eta$.

## State variables

The third-order equation in $f$ and the second-order equation in $\tau$ yield a
fifth-order system. Define:

$$y_1 = f, \qquad y_2 = f', \qquad y_3 = f'', \qquad y_4 = \tau, \qquad y_5 = \tau'$$

## First-order ODE system

Using the state variables and the expressions above, with $T = T_e y_4$ and
$C = \mu(T)/(\mu_e y_4)$:

!!! info ""

    $$y_1' = y_2$$

    $$y_2' = y_3$$

    $$y_3' = \frac{y_5}{y_4}\,y_3
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5 y_3
    - \frac{y_1 y_3 + \beta_H(y_4 - y_2^2)}{C}$$

    $$y_4' = y_5$$

    $$y_5' = \frac{y_5^2}{y_4}
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5^2
    - \frac{\mathrm{Pr}\,y_1 y_5}{C}
    - \mathrm{Pr}(\gamma-1)M_e^2\!\left[y_3^2 - \frac{\beta_H y_4 y_2}{C}\right]$$

## Boundary conditions

$$\eta = 0:\quad f = 0,\quad f' = 0,\quad
\tau = \tau_w\ \text{(isothermal)}$$

$$\eta \to \infty:\quad f' \to 1,\quad \tau \to 1$$

For an adiabatic wall, replace $\tau = \tau_w$ at $\eta = 0$ with $\tau' = 0$ at $\eta = 0$.

