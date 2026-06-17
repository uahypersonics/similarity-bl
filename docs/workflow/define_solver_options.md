# Define Solver Options

Solver options describe how the similarity problem should be solved.

The most important option is the equation family:

| Equation family | Use case | Main profiles |
|---|---|---|
| `falkner_skan` | 2D flat plate or wedge-like flows | $u/u_e$, $T/T_e$, $\rho/\rho_e$ |
| `falkner_skan_cooke` | Swept flows with crossflow | $u/u_e$, $w/w_e$, $T/T_e$, $\rho/\rho_e$ |

Numerical options control the eta grid and shooting solve:

| Option | Meaning |
|---|---|
| `eta_max` | Maximum similarity coordinate |
| `n_points` | Number of eta grid points |
| `tolerance` | Shooting convergence tolerance |
| `max_iterations` | Maximum shooting iterations |
| `ode_method` | SciPy ODE integration method |

For most flat-plate and swept-plate studies, the defaults are a good starting
point. Increase `eta_max` if the profiles have not reached their edge values by
the end of the domain.
