# Run Solver

The solver combines the physical inputs and numerical options, then returns the
similarity solution and convergence information.

The native solution is expressed in similarity space using the coordinate
$\eta$.

Common output fields are:

| Field | Meaning |
|---|---|
| `eta` | Similarity coordinate |
| `fp` | Streamwise velocity ratio, $u/u_e$ |
| `tau` | Temperature ratio, $T/T_e$ |
| `1/tau` | Density ratio, $\rho/\rho_e$ |
| `g` | FSC crossflow ratio, $w/w_e$ |

The FSC crossflow field `g` is present for `falkner_skan_cooke` solutions. For
2D Falkner-Skan solutions, use `fp`, `tau`, and derived density.
