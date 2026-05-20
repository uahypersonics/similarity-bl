"""Test that v = (xi/xi_ref)^m is computed correctly for flat plate (m=0)"""

from simbl import SimilarityInputs, SolverOptions, solve_similarity
import numpy as np

print("Testing v = (xi/xi_ref)^m for flat plate (beta=0, m=0)")
print("=" * 70)
print("\nFor flat plate: m = 0, so v = (xi/xi_ref)^0 = 1 for any xi/xi_ref")
print("This means f''(0) and tau(0) should be independent of xi/xi_ref\n")

results = []

for xi_ratio in [0.5, 1.0, 2.0, 5.0]:
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc='adiabatic',
        beta=0.0,  # flat plate
        xi_over_xi_ref=xi_ratio,
    )

    options = SolverOptions(
        equations='falkner_skan_cooke',
        n_points=150,
        eta_max=8.0,
    )

    solution, info = solve_similarity(inputs, options)

    # For m=0, v = xi_ratio^0 = 1 for any xi_ratio
    m = 0.0
    v = 1.0

    print(f"xi/xi_ref = {xi_ratio:.1f}:")
    print(f"  m = {m:.4f}, v = {v:.6f}")
    print(f"  f''(0) = {solution.fpp[0]:.8f}")
    print(f"  tau(0) = {solution.tau[0]:.8f}")
    print(f"  Converged: {info.converged}")
    print()

    results.append((xi_ratio, solution.fpp[0], solution.tau[0]))

# Check that results are identical (within numerical tolerance)
print("=" * 70)
fpp_values = [r[1] for r in results]
tau_values = [r[2] for r in results]

fpp_std = np.std(fpp_values)
tau_std = np.std(tau_values)

print(f"Standard deviation of f''(0): {fpp_std:.2e}")
print(f"Standard deviation of tau(0): {tau_std:.2e}")

if fpp_std < 1e-10 and tau_std < 1e-10:
    print("\n✓ PASS: f''(0) and tau(0) are independent of xi/xi_ref for flat plate")
    print("  This confirms v = 1 for m = 0 as expected from Liu Eq. 23")
else:
    print("\n✗ FAIL: Results vary with xi/xi_ref")
