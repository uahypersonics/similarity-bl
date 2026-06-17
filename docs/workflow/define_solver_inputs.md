# Define Solver Inputs

Solver inputs describe the local physical similarity problem.

At minimum, a solve needs an edge Mach number and edge temperature. Most real
cases also specify a wall condition, gas constants, pressure-gradient parameter,
and viscosity model.

| Input | Meaning |
|---|---|
| `mach_edge` | Local edge Mach number |
| `temp_edge` | Local edge temperature in K |
| `wall_bc` | `adiabatic` or `isothermal` |
| `temp_wall` | Wall temperature in K for isothermal walls |
| `beta` | Hartree pressure-gradient parameter |
| `gamma` | Ratio of specific heats |
| `prandtl` | Prandtl number |
| `viscosity_model` | Dynamic viscosity model |
| `sweep_angle` | Sweep angle in degrees for FSC solves |

If edge conditions vary with streamwise location, define a separate local input
state at each station. `simbl` then solves a locally self-similar problem for
each station rather than one fully non-similar streamwise evolution.
