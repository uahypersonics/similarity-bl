"""Test if S and K affect the solution for non-ZPG case"""

from simbl import SimilarityInputs, SolverOptions, solve_similarity
import numpy as np

print("Testing S and K effects with NON-ZERO pressure gradient")
print("=" * 70)

beta = 0.3  # Non-zero pressure gradient

for lam in [0.0, 30.0]:
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc='adiabatic',
        beta=beta,
        sweep_angle=lam,
    )

    options = SolverOptions(
        equations='falkner_skan_cooke',
        n_points=200,
        eta_max=10.0,
    )

    solution, info = solve_similarity(inputs, options)

    # Compute S and K
    cos2_Lambda = np.cos(np.radians(lam))**2
    S = 1.0 + 0.2 * 4.0 * cos2_Lambda
    K = 2.6 / (1.0 + 0.2 * 4.0 * cos2_Lambda)

    print(f"\nLambda = {lam}°, beta = {beta}:")
    print(f"  S = {S:.6f}, K = {K:.6f}")
    print(f"  hartree_beta * S = {beta * S:.6f}")
    print(f"  Converged: {info.converged}")
    print(f"  f''(0) = {solution.fpp[0]:.8f}")
    print(f"  gcf_p(0) = {solution.gcf_p[0]:.8f}")
    print(f"  tau(0) = {solution.tau[0]:.8f}")

# Compare
print("\n" + "=" * 70)
sol_0, _ = solve_similarity(
    SimilarityInputs(mach_edge=2.0, temp_edge=300.0, wall_bc='adiabatic', beta=beta, sweep_angle=0.0),
    SolverOptions(equations='falkner_skan_cooke', n_points=200, eta_max=10.0)
)
sol_30, _ = solve_similarity(
    SimilarityInputs(mach_edge=2.0, temp_edge=300.0, wall_bc='adiabatic', beta=beta, sweep_angle=30.0),
    SolverOptions(equations='falkner_skan_cooke', n_points=200, eta_max=10.0)
)

print(f"Difference in f''(0): {sol_30.fpp[0] - sol_0.fpp[0]:.10e}")
print(f"Difference in tau(0): {sol_30.tau[0] - sol_0.tau[0]:.10e}")

if abs(sol_30.tau[0] - sol_0.tau[0]) > 1e-6:
    print("\n✓ tau(0) DOES depend on sweep angle for non-ZPG!")
    print("  This confirms S and K are working in the momentum/energy equations.")
else:
    print("\n✗ tau(0) is still independent of sweep angle!")
    print("  This suggests a bug in the implementation.")
