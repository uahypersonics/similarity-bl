"""Templates for generated SIMBL configuration files.

Two templates:
    CONFIG_TEMPLATE_FS  - 2D Falkner-Skan (flat plate, wedge)
    CONFIG_TEMPLATE_FSC - 3D Falkner-Skan-Cooke (swept wing, adds sweep_angle)

CONFIG_TEMPLATE is an alias for CONFIG_TEMPLATE_FS (default for `simbl init`).
"""

# --------------------------------------------------
# Falkner-Skan (2D) template
# --------------------------------------------------

CONFIG_TEMPLATE_FS = """\
# simbl configuration file for Falkner-Skan (2D)
# run: simbl solve simbl_config.toml

# equations: "falkner_skan" or "fs" for 2D, "falkner_skan_cooke" or "fsc" for 3D
equations = "falkner_skan"

# --------------------------------------------------
# conditions
# --------------------------------------------------

[conditions]
# edge mach number
mach_edge = 6.0
# edge temperature [K]
temp_edge = 300.0
# hartree parameter (beta = 2m/(m+1)))
beta = 0.0

# --------------------------------------------------
# gas properties
# --------------------------------------------------

[gas]
# ratio of specific heats
gamma = 1.4
# Prandtl number
prandtl = 0.72

# --------------------------------------------------
# wall boundary condition
# --------------------------------------------------

[wall]
# wall BC: "adiabatic" or "isothermal"
type = "adiabatic"
# temp_wall = 600.0  # wall temperature [K] (isothermal only)

# --------------------------------------------------
# viscosity model
# --------------------------------------------------

[viscosity]
# available models: sutherland, sutherland_low_temp, sutherland_blended, keyes, power_law
model = "sutherland"

# --------------------------------------------------
# numerical solver options
# --------------------------------------------------

[numerics]
# maximum similarity coordinate (eta = 8-20 usually sufficient)
eta_max = 15.0
# number of grid points
n_points = 500
# shooting convergence tolerance
tolerance = 1e-8
# maximum shooting iterations
max_iterations = 100
# ode integrator method passed to scipy.integrate.solve_ivp (LSODA, RK45, RK23, DOP853, Radau, BDF)
ode_method = "LSODA"
# maximum wall-clock time [s] for one solve call (increase for difficult cases)
max_solve_time = 20.0
"""

# --------------------------------------------------
# Falkner-Skan-Cooke (3D) template
# --------------------------------------------------

CONFIG_TEMPLATE_FSC = """\
# simbl configuration file for Falkner-Skan-Cooke (3D)
# run: simbl solve simbl_config.toml

# equations: "falkner_skan" or "fs" for 2D, "falkner_skan_cooke" or "fsc" for 3D
equations = "falkner_skan_cooke"

# --------------------------------------------------
# conditions
# --------------------------------------------------

[conditions]
# edge mach number
mach_edge = 6.0
# edge temperature [K]
temp_edge = 300.0
# hartree parameter (beta = 2m/(m+1)))
beta = 0.0
# sweep angle [deg] (0 = 2D, must be < 90)
sweep_angle = 45.0

# --------------------------------------------------
# gas properties
# --------------------------------------------------

[gas]
# ratio of specific heats
gamma = 1.4
# Prandtl number
prandtl = 0.72

# --------------------------------------------------
# wall boundary condition
# --------------------------------------------------

[wall]
# wall BC: "adiabatic" or "isothermal"
type = "adiabatic"
# temp_wall = 600.0  # wall temperature [K] (isothermal only)

# --------------------------------------------------
# viscosity model
# --------------------------------------------------

[viscosity]
# available models: sutherland, sutherland_low_temp, sutherland_blended, keyes, power_law
model = "sutherland"

# --------------------------------------------------
# numerical solver options
# --------------------------------------------------

[numerics]
# maximum similarity coordinate (eta = 8-20 usually sufficient)
eta_max = 15.0
# number of grid points
n_points = 500
# shooting convergence tolerance
tolerance = 1e-8
# maximum shooting iterations
max_iterations = 100
# ode integrator method passed to scipy.integrate.solve_ivp (LSODA, RK45, RK23, DOP853, Radau, BDF)
ode_method = "LSODA"
# maximum wall-clock time [s] for one solve call (increase for difficult cases)
max_solve_time = 20.0
"""

# --------------------------------------------------
# default alias
# --------------------------------------------------

# CONFIG_TEMPLATE defaults to FS (2D); used by simbl init when no --model is given
CONFIG_TEMPLATE = CONFIG_TEMPLATE_FS
