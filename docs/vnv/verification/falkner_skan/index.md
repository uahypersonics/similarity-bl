# Falkner-Skan Verification

Analytic-limit tests that confirm the Falkner-Skan solver is solving the
governing equations correctly.

| Case | Description |
|---|---|
| [Blasius limit](blasius.md) | $M_e \to 0$, adiabatic, $\beta=0$ — recovers classical Blasius $f''(0)$ |
| [Crocco relation](crocco.md) | $\Pr=1$, adiabatic — pointwise $\tau(\eta)$ matches Crocco's theorem |
| [Hiemenz stagnation point](hiemenz/index.md) | $\beta=1$, $M_e \to 0$ — recovers stagnation-point wall shear |
| [Hartree table](hartree/index.md) | $\beta = 0, 0.2, 0.5, 1.0$, $M_e \to 0$ — matches Hartree (1937) tabulated $f''(0)$ |
| [Separation point](separation/index.md) | $\beta \approx -0.1988$ — $f''(0) \to 0$ at Hartree separation parameter |

