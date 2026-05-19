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
# SIMBL Configuration File - Falkner-Skan (2D)
# Compressible similarity boundary layer on a flat plate or wedge surface.
# Run: simbl solve simbl_config.toml

equations = "falkner_skan"

[conditions]
mach_edge = 6.0
temp_edge = 300.0  # edge temperature [K]
beta = 0.0         # Falkner-Skan pressure gradient parameter
               #   beta = 0.0  -> flat plate (Blasius)
               #   beta > 0.0  -> favorable gradient (wedge flow)
               #   beta < 0.0  -> adverse gradient (expanding flow)
gamma = 1.4
prandtl = 0.72

[wall]
# wall BC: "adiabatic" or "isothermal"
type = "adiabatic"
# temp_wall = 600.0  # wall temperature [K] (isothermal only)

[viscosity]
# available models: sutherland, sutherland_low_temp, sutherland_blended, keyes, power_law
model = "sutherland"
# optional: override built-in air() preset constants for sutherland
# mu_ref = 1.716e-5  # reference viscosity [Pa s]
# T_ref  = 273.15   # reference temperature [K]
# S      = 110.4    # Sutherland constant [K]

[numerics]
eta_max = 15.0       # maximum similarity coordinate (eta = 8-12 usually sufficient)
n_points = 500       # number of grid points
tolerance = 1e-8     # shooting convergence tolerance
max_iterations = 100 # maximum shooting iterations
ode_method = "LSODA" # ODE integrator (LSODA, RK45, BDF)
"""

# --------------------------------------------------
# Falkner-Skan-Cooke (3D swept wing) template
# --------------------------------------------------

CONFIG_TEMPLATE_FSC = """\
# SIMBL Configuration File - Falkner-Skan-Cooke (3D swept wing)
# Compressible similarity boundary layer on an infinite swept wing.
# Adds sweep_angle to the 2D Falkner-Skan problem, which introduces
# a crossflow velocity component w(eta) normal to the streamwise direction.
# Run: simbl solve simbl_config.toml

equations = "falkner_skan_cooke"

[conditions]
mach_edge = 6.0
temp_edge = 300.0  # edge temperature [K]
beta = 0.0         # Falkner-Skan pressure gradient parameter
               #   beta = 0.0  -> flat plate (Blasius)
               #   beta > 0.0  -> favorable gradient (wedge flow)
               #   beta < 0.0  -> adverse gradient (expanding flow)
sweep_angle = 45.0 # wing sweep angle [deg] (0 = 2D, must be < 90)
               # crossflow grows with sweep angle
gamma = 1.4
prandtl = 0.72

[wall]
# wall BC: "adiabatic" or "isothermal"
type = "adiabatic"
# temp_wall = 600.0  # wall temperature [K] (isothermal only)

[viscosity]
# available models: sutherland, sutherland_low_temp, sutherland_blended, keyes, power_law
model = "sutherland"
# optional: override built-in air() preset constants for sutherland
# mu_ref = 1.716e-5  # reference viscosity [Pa s]
# T_ref  = 273.15   # reference temperature [K]
# S      = 110.4    # Sutherland constant [K]

[numerics]
eta_max = 15.0       # maximum similarity coordinate (eta = 8-12 usually sufficient)
n_points = 500       # number of grid points
tolerance = 1e-8     # shooting convergence tolerance
max_iterations = 100 # maximum shooting iterations
ode_method = "LSODA" # ODE integrator (LSODA, RK45, BDF)
"""

# --------------------------------------------------
# default alias
# --------------------------------------------------

# CONFIG_TEMPLATE defaults to FS (2D); used by simbl init when no --model is given
CONFIG_TEMPLATE = CONFIG_TEMPLATE_FS
