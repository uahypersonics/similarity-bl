"""Debug script to check energy equation terms"""

from simbl import SimilarityInputs, SolverOptions
from simbl.solver.falkner_skan_cooke.ode import bl_ode
from simbl.solver.falkner_skan_cooke.builder import build_solver_problem
import numpy as np
from flow_state import get_transport_model

# Setup problem
inputs = SimilarityInputs(
    mach_edge=4.0,
    temp_edge=300.0,
    wall_bc='adiabatic',
    sweep_angle=30.0,
)

# Get transport model
visc_model = get_transport_model(inputs.viscosity_model, inputs.viscosity_model_kwargs)

# Test at a point in the boundary layer where crossflow is developing
eta = 2.0
# State: [f, fp, fpp, g_cf, gcf_p, tau, tau_p]
# Use approximate values from solution
y = np.array([0.5, 0.8, 0.3, 0.95, 0.1, 2.5, 0.0])

sweep_rad = np.radians(inputs.sweep_angle)
tan2_sweep = np.tan(sweep_rad)**2

# Call ODE function
dy = bl_ode(eta, y, inputs, visc_model, tan2_sweep)

print("Energy Equation Debug at eta = 2.0")
print("=" * 70)
print(f"State vector:")
print(f"  f = {y[0]:.6f}, fp = {y[1]:.6f}, fpp = {y[2]:.6f}")
print(f"  g_cf = {y[3]:.6f}, gcf_p = {y[4]:.6f}")
print(f"  tau = {y[5]:.6f}, tau_p = {y[6]:.6f}")
print()
print(f"Computed derivatives:")
print(f"  dy[6] (tau'') = {dy[6]:.10e}")
print()

# Now manually compute S and K and check the terms
M = inputs.mach_edge
gamma = inputs.gamma
cos_Lambda = np.cos(sweep_rad)
cos2_Lambda = cos_Lambda**2
v = 1.0

S = 1.0 + (gamma - 1.0) / 2.0 * v**2 * M**2 * cos2_Lambda
numerator = 1.0 + (gamma - 1.0) / 2.0 * M**2
denominator = 1.0 + (gamma - 1.0) / 2.0 * M**2 * cos2_Lambda
K = numerator / denominator

print(f"Liu parameters:")
print(f"  S = {S:.10f}")
print(f"  K = {K:.10f}")
print(f"  (S - 1) = {S - 1.0:.10f}")
print(f"  (K - 1) = {K - 1.0:.10f}")
print()

# Check crossflow contribution
gcf2_prime = 2.0 * y[3] * y[4]
print(f"Crossflow kinetic energy:")
print(f"  (g_cf^2)' = 2*g_cf*gcf_p = {gcf2_prime:.10e}")
print(f"  (K - 1) * S = {(K - 1.0) * S:.10f}")
print(f"  Expected crossflow contribution: (K-1)*S*N*(g_cf^2)' ~ {(K-1.0)*S*gcf2_prime:.10e}")
