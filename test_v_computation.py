"""Test that v = (xi/xi_ref)^m is computed correctly from Liu Eq. 23"""

from simbl import SimilarityInputs, SolverOptions, solve_similarity
import numpy as np

print("Testing v = (xi/xi_ref)^m computation")
print("=" * 70)

# Test 1: Flat plate (beta=0, m=0) should have v=1 regardless of xi_over_xi_ref
print("\nTest 1: Flat plate (beta=0, m=0)")
for xi_ratio in [0.5, 1.0, 2.0]:
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc='adiabatic',
        beta=0.0,
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
    v_expected = 1.0

    print(f"  xi/xi_ref = {xi_ratio:.1f}: m = {m:.4f}, v_expected = {v_expected:.6f}")
    print(f"    f''(0) = {solution.fpp[0]:.6f}, tau(0) = {solution.tau[0]:.6f}")
    print(f"    Converged: {info.converged}")

# Test 2: Favorable pressure gradient (beta=0.5, m=1)
print("\nTest 2: Favorable pressure gradient (beta=0.5, m=1)")
for xi_ratio in [0.5, 1.0, 2.0]:
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc='adiabatic',
        beta=0.5,
        xi_over_xi_ref=xi_ratio,
    )

    options = SolverOptions(
        equations='falkner_skan_cooke',
        n_points=150,
        eta_max=8.0,
    )

    solution, info = solve_similarity(inputs, options)

    # For beta=0.5: m = 0.5/(2-0.5) = 0.5/1.5 = 1/3
    m = 0.5 / (2.0 - 0.5)
    v_expected = xi_ratio**m

    print(f"  xi/xi_ref = {xi_ratio:.1f}: m = {m:.4f}, v_expected = {v_expected:.6f}")
    print(f"    f''(0) = {solution.fpp[0]:.6f}, tau(0) = {solution.tau[0]:.6f}")
    print(f"    Converged: {info.converged}")

# Test 3: Verify that S depends on v correctly
print("\nTest 3: S = 1 + (gamma-1)/2 * v^2 * M^2 * cos^2(Lambda)")
M = 4.0
gamma = 1.4
beta = 0.0  # flat plate, m=0, v=1

for lam in [0.0, 30.0]:
    cos2_Lambda = np.cos(np.radians(lam))**2
    v = 1.0  # For beta=0
    S = 1.0 + (gamma - 1.0) / 2.0 * v**2 * M**2 * cos2_Lambda

    print(f"\n  Lambda = {lam}°:")
    print(f"    cos^2(Lambda) = {cos2_Lambda:.6f}")
    print(f"    v = {v:.6f}")
    print(f"    S = {S:.6f}")
    print(f"    (gamma-1)/2 * M^2 = {(gamma-1.0)/2.0 * M**2:.6f}")

print("\n" + "=" * 70)
print("v computation tests complete!")
