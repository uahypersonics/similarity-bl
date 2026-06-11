# Reduction of Order

Starting from the [FSC ODEs](index.md):

$$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

$$(Cg')' + fg' = 0$$

$$\begin{aligned}
&\left(\frac{C}{\mathrm{Pr}}\tau'\right)'
+ (S-1)\left(C(f'^2)'\right)'
+ (K-1)S\left(C(g^2)'\right)' \\
&\quad + f\!\left[\tau' + (S-1)(f'^2)' + (K-1)S(g^2)'\right] = 0
\end{aligned}$$

where $C = \rho\mu/(\rho_e\mu_e)$, $f' = u/u_e$, $g = w/w_e$, $\tau = T/T_e$, and all
primes denote differentiation with respect to $\eta$.

## Chapman–Rubesin factor

For a perfect gas $\rho/\rho_e = 1/\tau$, so the Chapman–Rubesin factor simplifies to:

$$C = \frac{\rho\mu}{\rho_e\mu_e} = \frac{\mu}{\mu_e\tau}$$

Differentiating with respect to $\eta$ via the quotient rule:

$$C' = \frac{d}{d\eta}\!\left(\frac{\mu}{\mu_e\tau}\right)
= \frac{1}{\mu_e}\frac{\mu'\tau - \mu\tau'}{\tau^2}$$

and therefore:

$$\frac{C'}{C} = \frac{\mu'}{\mu} - \frac{\tau'}{\tau}$$

## Expanding the compound derivatives

Applying the product rule to each term:

$$(Cf'')' = C'f'' + Cf'''$$

$$(Cg')' = C'g' + Cg''$$

$$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' = \frac{C'}{\mathrm{Pr}}\tau' + \frac{C}{\mathrm{Pr}}\tau''$$

$$\left(C(f'^2)'\right)' = \left(2Cf'f''\right)' = 2\!\left(C'f'f'' + Cf''^2 + Cf'f'''\right)$$

$$\left(C(g^2)'\right)' = \left(2Cgg'\right)' = 2\!\left(C'gg' + Cg'^2 + Cgg''\right)$$

Substituting into the three ODEs:

**x-momentum:**

$$C'f'' + Cf''' + ff'' + \beta_H(\tau - f'^2) = 0$$

**z-momentum:**

$$C'g' + Cg'' + fg' = 0$$

**Energy:**

$$\frac{C'}{\mathrm{Pr}}\tau' + \frac{C}{\mathrm{Pr}}\tau''
+ 2(S-1)\!\left(C'f'f'' + Cf''^2 + Cf'f'''\right)
+ 2(K-1)S\!\left(C'gg' + Cg'^2 + Cgg''\right)
+ f\!\left[\tau' + 2(S-1)f'f'' + 2(K-1)Sgg'\right] = 0$$

## Isolating the highest derivatives

### x-momentum → $f'''$

Divide by $C$ and rearrange, using $C'/C = \mu'/\mu - \tau'/\tau$:

$$\boxed{f''' = \frac{\tau'}{\tau}\,f'' - \frac{\mu'}{\mu}\,f''
- \frac{ff''}{C} - \frac{\beta_H(\tau - f'^2)}{C}}$$

### z-momentum → $g''$

Divide by $C$ and rearrange, using $C'/C = \mu'/\mu - \tau'/\tau$:

$$\boxed{g'' = \frac{\tau'}{\tau}\,g' - \frac{\mu'}{\mu}\,g' - \frac{fg'}{C}}$$

### Energy → $\tau''$

Solve for $\tau''$, multiply by $\mathrm{Pr}/C$, and substitute $C'/C = \mu'/\mu - \tau'/\tau$:

$$\boxed{\begin{aligned}
\tau'' &= \frac{\tau'^2}{\tau} - \frac{\mu'\tau'}{\mu} - \frac{\mathrm{Pr}\,f\tau'}{C} \\
&\quad - 2\mathrm{Pr}(S-1)\!\left[\left(\frac{\mu'}{\mu} - \frac{\tau'}{\tau}\right)f'f'' + f''^2 + f'f'''\right] \\
&\quad - 2\mathrm{Pr}(K-1)S\!\left[\left(\frac{\mu'}{\mu} - \frac{\tau'}{\tau}\right)gg' + g'^2 + gg''\right] \\
&\quad - \frac{2\mathrm{Pr}\,f}{C}\!\left[(S-1)f'f'' + (K-1)Sgg'\right]
\end{aligned}}$$

## Viscosity derivatives

The $\mu'/\mu$ term requires the derivative of viscosity with respect to $\eta$.
By the chain rule:

$$\mu' = \frac{d\mu}{d\eta} = \frac{d\mu}{dT}\frac{dT}{d\eta}$$

Since $T = T_e\tau(\eta)$:

$$\frac{dT}{d\eta} = T_e\frac{d\tau}{d\eta} = T_e\tau'$$

Therefore:

$$\mu' = \frac{d\mu}{dT}\,T_e\tau'$$

Substituting $\mu'/\mu = (T_e/\mu)(d\mu/dT)\,\tau'$ into the intermediate expressions gives the final forms:

!!! success ""

    $$f''' = \frac{\tau'}{\tau}\,f''
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,\tau'f''
    - \frac{ff''}{C}
    - \frac{\beta_H(\tau - f'^2)}{C}$$

    $$g'' = \frac{\tau'}{\tau}\,g'
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,\tau'g'
    - \frac{fg'}{C}$$

    $$\begin{aligned}
    \tau'' &= \frac{\tau'^2}{\tau} - \frac{T_e}{\mu}\frac{d\mu}{dT}\,\tau'^2 - \frac{\mathrm{Pr}\,f\tau'}{C} \\
    &\quad - 2\mathrm{Pr}(S-1)\!\left[\left(\frac{T_e}{\mu}\frac{d\mu}{dT} - \frac{1}{\tau}\right)\tau'f'f'' + f''^2 + f'f'''\right] \\
    &\quad - 2\mathrm{Pr}(K-1)S\!\left[\left(\frac{T_e}{\mu}\frac{d\mu}{dT} - \frac{1}{\tau}\right)\tau'gg' + g'^2 + gg''\right] \\
    &\quad - \frac{2\mathrm{Pr}\,f}{C}\!\left[(S-1)f'f'' + (K-1)Sgg'\right]
    \end{aligned}$$

where $\mu = \mu(T_e\tau)$ and $C = \mu/(\mu_e\tau)$ are evaluated locally at each $\eta$,
and $f'''$ and $g''$ are computed first.

## State variables

The third-order equation in $f$, the second-order equation in $g$, and the second-order equation in
$\tau$ yield a seventh-order system. Define:

$$y_1 = f, \quad y_2 = f', \quad y_3 = f'', \quad y_4 = \tau, \quad y_5 = \tau', \quad y_6 = g, \quad y_7 = g'$$

## First-order ODE system

!!! info ""

    $$y_1' = y_2$$

    $$y_2' = y_3$$

    $$y_3' = \frac{y_5}{y_4}\,y_3
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5 y_3
    - \frac{y_1 y_3 + \beta_H(y_4 - y_2^2)}{C}$$

    $$y_4' = y_5$$

    $$\begin{aligned}
    y_5' &= \frac{y_5^2}{y_4} - \frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5^2 - \frac{\mathrm{Pr}\,y_1 y_5}{C} \\
    &\quad - 2\mathrm{Pr}(S-1)\!\left[\left(\frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5 - \frac{y_5}{y_4}\right)y_2 y_3 + y_3^2 + y_2 y_3'\right] \\
    &\quad - 2\mathrm{Pr}(K-1)S\!\left[\left(\frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5 - \frac{y_5}{y_4}\right)y_6 y_7 + y_7^2 + y_6 y_7'\right] \\
    &\quad - \frac{2\mathrm{Pr}\,y_1}{C}\!\left[(S-1)y_2 y_3 + (K-1)S\,y_6 y_7\right]
    \end{aligned}$$

    $$y_6' = y_7$$

    $$y_7' = \frac{y_5}{y_4}\,y_7
    - \frac{T_e}{\mu}\frac{d\mu}{dT}\,y_5 y_7
    - \frac{y_1 y_7}{C}$$

where $C = \mu(T_e y_4)/(\mu_e y_4)$, $T = T_e y_4$, and $y_3'$, $y_7'$ are computed before $y_5'$.

## Boundary conditions

**Wall** ($\eta = 0$):

$$y_1 = 0, \qquad y_2 = 0, \qquad y_6 = 0$$

- **Isothermal**: $y_4 = T_w/T_e$
- **Adiabatic**: $y_5 = 0$

**Edge** ($\eta \to \infty$):

$$y_2 = 1, \qquad y_4 = 1, \qquad y_6 = 1$$

The unknown wall values $y_3(0) = f''(0)$, $y_7(0) = g'(0)$, and (for isothermal walls)
$y_5(0) = \tau'(0)$ are determined by shooting to match the edge conditions.

