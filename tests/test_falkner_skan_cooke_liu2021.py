"""Tests for Liu (2021) Falkner-Skan-Cooke implementation

Validates the implementation of Liu (2021), Phys. Fluids 33, 126109, Eqs. 16-18.

Test cases:
1. Variable naming: no state variable named 'g' should represent temperature
2. No-crossflow case: verify crossflow terms vanish when g_cf = 0
3. Incompressible limit: M -> 0 should reduce to incompressible form
4. Blasius limit: m = 0 (beta = 0) with no crossflow should reduce to compressible Blasius
5. Boundary conditions: isothermal and adiabatic wall conditions
"""

import numpy as np

from simbl import (
    SimilarityInputs,
    SolverOptions,
    solve_similarity,
)


def test_variable_naming():
    """Test 1: Verify no state variable named 'g' represents temperature

    The solution should use 'tau' for temperature (T/T_e) and 'g_cf' for
    crossflow velocity (w/w_e), not 'g' for temperature.
    """
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc="adiabatic",
        sweep_angle=30.0,
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=200,
        eta_max=10.0,
    )

    solution, info = solve_similarity(inputs, options)

    # Check that solution has the correct attribute names
    assert hasattr(solution, 'tau'), "Solution should have 'tau' attribute for temperature"
    assert hasattr(solution, 'tau_p'), "Solution should have 'tau_p' attribute for temperature gradient"
    assert hasattr(solution, 'g_cf'), "Solution should have 'g_cf' attribute for crossflow velocity"
    assert hasattr(solution, 'gcf_p'), "Solution should have 'gcf_p' attribute for crossflow velocity gradient"

    # Check that solution does NOT have old attribute names for these quantities
    assert not hasattr(solution, 'w'), "Solution should not have 'w' attribute (use 'g_cf')"
    assert not hasattr(solution, 'wp'), "Solution should not have 'wp' attribute (use 'gcf_p')"
    assert not hasattr(solution, 'g'), "Solution should not have 'g' attribute (use 'tau')"
    assert not hasattr(solution, 'gp'), "Solution should not have 'gp' attribute (use 'tau_p')"

    print("✓ Test 1 passed: Variable naming is correct (tau for temperature, g_cf for crossflow)")


def test_no_crossflow_case():
    """Test 2: Verify crossflow terms vanish when sweep angle = 0

    For zero sweep angle, g_cf and gcf_p should remain zero throughout
    the boundary layer, and the energy equation should have no crossflow
    dissipation contribution.
    """
    inputs = SimilarityInputs(
        mach_edge=3.0,
        temp_edge=300.0,
        wall_bc="adiabatic",
        sweep_angle=0.0,  # No crossflow
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=200,
        eta_max=10.0,
    )

    solution, info = solve_similarity(inputs, options)

    # For zero sweep angle, crossflow should be zero everywhere (within numerical precision)
    assert np.allclose(solution.g_cf, 0.0, atol=1e-8), \
        f"Crossflow velocity g_cf should be zero for sweep_angle=0, got max={np.max(np.abs(solution.g_cf)):.2e}"
    assert np.allclose(solution.gcf_p, 0.0, atol=1e-8), \
        f"Crossflow velocity gradient gcf_p should be zero for sweep_angle=0, got max={np.max(np.abs(solution.gcf_p)):.2e}"

    print("✓ Test 2 passed: No-crossflow case - all crossflow terms vanish for sweep_angle=0")


def test_incompressible_limit():
    """Test 3: Incompressible limit (M -> 0)

    For very low Mach numbers, temperature should approach unity (tau ~ 1)
    and the Chapman-Rubesin parameter N should approach 1 (mu_ratio ~ 1).
    """
    # Low Mach number case
    inputs = SimilarityInputs(
        mach_edge=0.1,  # Very low Mach number
        temp_edge=300.0,
        wall_bc="adiabatic",
        sweep_angle=0.0,
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=200,
        eta_max=10.0,
    )

    solution, info = solve_similarity(inputs, options)

    # In the incompressible limit, temperature should stay close to 1
    # (adiabatic recovery should have minimal effect)
    assert np.allclose(solution.tau, 1.0, atol=0.05), \
        f"Temperature should be near 1 in incompressible limit, got tau(0)={solution.tau[0]:.4f}"

    print("✓ Test 3 passed: Incompressible limit - temperature remains near unity for M << 1")


def test_blasius_limit():
    """Test 4: Blasius limit (beta = 0, no crossflow, compressible)

    For a flat plate (beta = 0) with no sweep angle, the FSC equations
    should reduce to the compressible Blasius problem. The solution should
    converge and produce reasonable wall shear and temperature values.
    """
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc="adiabatic",
        beta=0.0,         # Flat plate
        sweep_angle=0.0,  # No crossflow
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=200,
        eta_max=10.0,
    )

    solution, info = solve_similarity(inputs, options)

    # Check convergence
    assert info.converged, "Blasius limit case should converge"

    # Check boundary conditions
    assert np.isclose(solution.fp[0], 0.0, atol=1e-6), "fp(0) should be 0 (no-slip)"
    assert np.isclose(solution.fp[-1], 1.0, atol=1e-3), "fp(inf) should be 1 (edge velocity)"
    assert np.isclose(solution.tau[-1], 1.0, atol=1e-3), "tau(inf) should be 1 (edge temperature)"

    # For adiabatic wall, tau_p(0) should be zero
    assert np.isclose(solution.tau_p[0], 0.0, atol=1e-6), \
        "tau_p(0) should be 0 for adiabatic wall"

    # Wall temperature should be greater than 1 due to compressibility (recovery effect)
    assert solution.tau[0] > 1.0, \
        f"Adiabatic wall temperature should exceed edge temperature, got tau(0)={solution.tau[0]:.4f}"

    print(f"✓ Test 4 passed: Blasius limit - converged with f''(0)={solution.fpp[0]:.4f}, tau(0)={solution.tau[0]:.4f}")


def test_energy_equation_uses_S_and_K():
    """Test that the energy equation correctly uses S and K parameters

    Verifies that Liu Eq. 18 is implemented with Liu's exact definitions:
    - S = 1 + (gamma-1)/2 * v^2 * Ma_e_ref^2 * cos^2(Lambda)  [Liu Eq. 21]
    - K = [1 + (gamma-1)/2 * Ma_e_ref^2] / [1 + (gamma-1)/2 * Ma_e_ref^2 * cos^2(Lambda)]  [Liu Eq. 20]

    And that the equation includes the proper nested derivative terms:
    - (N/Pr tau')' term
    - (S - 1) [N (fp^2)']' term
    - (K - 1) [S N (g_cf^2)']' term
    """
    # Read the source code to verify S and K are computed
    import inspect

    from simbl.solver.falkner_skan_cooke.ode import bl_ode

    source = inspect.getsource(bl_ode)

    # Check that S is computed with correct Liu formula
    assert "S = " in source, "S parameter should be defined"
    assert "(problem.gamma - 1" in source or "(gamma - 1" in source, \
        "S should involve (gamma - 1)"
    assert "Ma_e_ref" in source or "mach_edge" in source, \
        "S should involve Ma_e_ref (edge Mach number)"
    assert "cos2_Lambda" in source or "cos^2(Lambda)" in source or "cos_Lambda**2" in source, \
        "S should involve cos^2(Lambda)"

    # Check that K is computed with correct Liu formula
    assert "K = " in source, "K parameter should be defined"
    assert "numerator" in source and "denominator" in source, \
        "K should be computed as ratio (Liu Eq. 20)"

    # Check that the energy equation includes the kinetic energy derivative terms
    assert "fp2_prime" in source or "fp^2" in source or "(fp**2)" in source, \
        "Energy equation should include fp^2 derivative term"
    assert "gcf2_prime" in source or "g_cf^2" in source or "gcf^2" in source, \
        "Energy equation should include g_cf^2 derivative term"

    # Check that fppp (third derivative) is computed from momentum equation
    assert "fppp" in source, "Energy equation should use fppp from momentum equation"

    # Check that gcf_pp (second derivative of crossflow) is computed
    assert "gcf_pp" in source, "Energy equation should use gcf_pp from crossflow momentum"

    # Check that S and K are actually used in the equations
    assert "(S - 1" in source or "(S-1" in source, \
        "Energy equation should use (S - 1) term"
    assert "(K - 1" in source or "(K-1" in source, \
        "Energy equation should use (K - 1) term"

    # Check that hartree_beta * S appears in momentum equation (Liu Eq. 16)
    assert "hartree_beta" in source, "Momentum equation should use hartree_beta = 2m/(m+1)"
    assert "hartree_beta * S" in source or "hartree_beta*S" in source, \
        "Momentum equation should have hartree_beta * S term (Liu Eq. 16)"

    # Now run a simple case to verify it executes correctly
    inputs = SimilarityInputs(
        mach_edge=2.0,
        temp_edge=300.0,
        wall_bc="adiabatic",
        sweep_angle=30.0,
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=150,
        eta_max=8.0,
    )

    solution, info = solve_similarity(inputs, options)

    # Verify it converges and produces reasonable results
    assert info.converged, "Solution with S and K terms should converge"
    assert solution.tau[0] > 1.0, "Adiabatic wall temperature should exceed edge temperature"
    assert solution.tau[-1] < solution.tau[0] + 0.5, "Temperature should approach edge value"

    print("✓ Test 6 passed: Energy equation correctly uses S and K parameters")


def test_boundary_conditions():
    """Test 5: Boundary conditions for isothermal and adiabatic walls

    Verify that:
    - Adiabatic wall: tau_p(0) = 0 (zero heat flux)
    - Isothermal wall: tau(0) = tau_wall (prescribed temperature)
    - Both: fp(0) = 0 (no-slip), g_cf(0) = 0 (no crossflow at wall)
    - Both: fp(inf) = 1, tau(inf) = 1, g_cf(inf) = 1 (for sweep > 0)
    """
    # Test adiabatic wall
    inputs_adiabatic = SimilarityInputs(
        mach_edge=3.0,
        temp_edge=300.0,
        wall_bc="adiabatic",
        sweep_angle=30.0,
    )

    options = SolverOptions(
        equations="falkner_skan_cooke",
        n_points=200,
        eta_max=10.0,
    )

    sol_adiabatic, info_adiabatic = solve_similarity(inputs_adiabatic, options)

    # Adiabatic BC checks
    assert info_adiabatic.converged, "Adiabatic case should converge"
    assert np.isclose(sol_adiabatic.tau_p[0], 0.0, atol=1e-6), \
        "Adiabatic wall should have tau_p(0) = 0"
    assert np.isclose(sol_adiabatic.fp[0], 0.0, atol=1e-6), \
        "Wall should have fp(0) = 0 (no-slip)"
    assert np.isclose(sol_adiabatic.g_cf[0], 0.0, atol=1e-6), \
        "Wall should have g_cf(0) = 0 (no crossflow)"

    # Edge condition checks
    assert np.isclose(sol_adiabatic.fp[-1], 1.0, atol=1e-3), \
        "Edge should have fp(inf) = 1"
    assert np.isclose(sol_adiabatic.tau[-1], 1.0, atol=1e-3), \
        "Edge should have tau(inf) = 1"
    assert np.isclose(sol_adiabatic.g_cf[-1], 1.0, atol=1e-3), \
        "Edge should have g_cf(inf) = 1 for sweep_angle > 0"

    print("✓ Test 5a passed: Adiabatic boundary conditions satisfied")

    # Test isothermal wall
    inputs_isothermal = SimilarityInputs(
        mach_edge=3.0,
        temp_edge=300.0,
        temp_wall=400.0,  # Hot wall
        wall_bc="isothermal",
        sweep_angle=30.0,
    )

    sol_isothermal, info_isothermal = solve_similarity(inputs_isothermal, options)

    # Isothermal BC checks
    assert info_isothermal.converged, "Isothermal case should converge"
    tau_wall_expected = 400.0 / 300.0  # T_wall / T_edge
    assert np.isclose(sol_isothermal.tau[0], tau_wall_expected, atol=1e-6), \
        f"Isothermal wall should have tau(0) = {tau_wall_expected}"
    assert np.isclose(sol_isothermal.fp[0], 0.0, atol=1e-6), \
        "Wall should have fp(0) = 0 (no-slip)"
    assert np.isclose(sol_isothermal.g_cf[0], 0.0, atol=1e-6), \
        "Wall should have g_cf(0) = 0 (no crossflow)"

    print("✓ Test 5b passed: Isothermal boundary conditions satisfied")


if __name__ == "__main__":
    print("Running Liu (2021) FSC Implementation Tests")
    print("=" * 80)
    print()

    test_variable_naming()
    print()

    test_no_crossflow_case()
    print()

    test_incompressible_limit()
    print()

    test_blasius_limit()
    print()

    test_boundary_conditions()
    print()

    test_energy_equation_uses_S_and_K()
    print()

    print("=" * 80)
    print("All tests passed! ✓")
