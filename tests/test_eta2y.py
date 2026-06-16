"""Tests for eta-to-y similarity coordinate transforms."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
import numpy as np
import pytest
from typer.testing import CliRunner

from simbl.cli import cli
from simbl.transform import eta2y


# --------------------------------------------------
# tests
# --------------------------------------------------
def test_eta2y_levy_lees_constant_tau() -> None:
    """Check Levy-Lees eta2y for constant tau against analytic scaling."""

    # build constant-temperature profile
    eta = np.array([0.0, 1.0, 2.0, 3.0])
    tau = np.ones_like(eta)

    # build dimensional edge state
    x = 2.0
    dens_edge = 4.0
    uvel_edge = 9.0
    visc_edge = 0.5

    # compute expected analytic inverse transform
    eta_scale = np.sqrt(dens_edge * uvel_edge / (2.0 * visc_edge * x))
    expected_y = eta / eta_scale

    # compute eta2y transform
    y = eta2y(
        eta=eta,
        tau=tau,
        x=x,
        dens_edge=dens_edge,
        uvel_edge=uvel_edge,
        visc_edge=visc_edge,
        transform="levy_lees",
    )

    # validate result
    np.testing.assert_allclose(y, expected_y)


def test_eta2y_illingworth_stewartson_integrates_tau() -> None:
    """Check Illingworth-Stewartson eta2y integrates tau with trapezoid rule."""

    # build profile with eta_scale equal to one
    eta = np.array([0.0, 1.0, 2.0])
    tau = np.array([1.0, 2.0, 3.0])
    expected_y = np.array([0.0, 1.5, 4.0])

    # compute eta2y transform
    y = eta2y(
        eta=eta,
        tau=tau,
        x=1.0,
        dens_edge=2.0,
        uvel_edge=1.0,
        visc_edge=1.0,
        transform="illingworth_stewartson",
    )

    # validate result
    np.testing.assert_allclose(y, expected_y)


def test_eta2y_levy_lees_and_illingworth_stewartson_use_distinct_scales() -> None:
    """Check LL and IS transforms can use distinct reference scales."""

    # build constant-temperature profile
    eta = np.array([0.0, 1.0, 2.0])
    tau = np.ones_like(eta)

    # compute both transforms with the same local edge state
    y_levy_lees = eta2y(
        eta=eta,
        tau=tau,
        x=2.0,
        dens_edge=4.0,
        uvel_edge=9.0,
        visc_edge=0.5,
        transform="levy_lees",
    )
    y_illingworth_stewartson = eta2y(
        eta=eta,
        tau=tau,
        x=2.0,
        dens_edge=4.0,
        uvel_edge=9.0,
        visc_edge=0.5,
        dens_ref=1.0,
        visc_ref=1.0,
        transform="illingworth_stewartson",
    )

    # validate the transforms are not accidentally sharing one scale formula
    assert not np.allclose(y_levy_lees, y_illingworth_stewartson)


def test_eta2y_omitted_transform_uses_fsc_default() -> None:
    """Check omitted transform maps FSC equations to Illingworth-Stewartson."""

    # build simple input profile
    eta = np.array([0.0, 1.0, 2.0])
    tau = np.array([1.0, 1.5, 2.0])

    # compute using omitted and explicit transform names
    y_default = eta2y(
        eta=eta,
        tau=tau,
        x=1.0,
        dens_edge=2.0,
        uvel_edge=1.0,
        visc_edge=1.0,
        equations="falkner_skan_cooke",
    )
    y_explicit = eta2y(
        eta=eta,
        tau=tau,
        x=1.0,
        dens_edge=2.0,
        uvel_edge=1.0,
        visc_edge=1.0,
        transform="illingworth_stewartson",
    )

    # validate dispatch is equivalent
    np.testing.assert_allclose(y_default, y_explicit)


def test_eta2y_rejects_invalid_transform() -> None:
    """Check invalid transform names fail clearly."""

    # build simple input profile
    eta = np.array([0.0, 1.0])
    tau = np.array([1.0, 1.0])

    # validate bad transform raises
    with pytest.raises(ValueError, match="Unknown eta2y transform"):
        eta2y(
            eta=eta,
            tau=tau,
            x=1.0,
            dens_edge=2.0,
            uvel_edge=1.0,
            visc_edge=1.0,
            transform="bad_transform",
        )


def test_eta2y_rejects_mismatched_shapes() -> None:
    """Check eta and tau must describe the same one-dimensional grid."""

    # build eta and tau arrays with incompatible shapes
    eta = np.array([0.0, 1.0, 2.0])
    tau = np.array([1.0, 1.0])

    # validate shape mismatch raises
    with pytest.raises(ValueError, match="eta and tau must have the same shape"):
        eta2y(
            eta=eta,
            tau=tau,
            x=1.0,
            dens_edge=2.0,
            uvel_edge=1.0,
            visc_edge=1.0,
            transform="levy_lees",
        )


def test_eta2y_omitted_transform_requires_equations() -> None:
    """Check omitted transform requires equation family."""

    # build simple input profile
    eta = np.array([0.0, 1.0])
    tau = np.array([1.0, 1.0])

    # validate missing equations raises
    with pytest.raises(ValueError, match="equations must be provided"):
        eta2y(
            eta=eta,
            tau=tau,
            x=1.0,
            dens_edge=2.0,
            uvel_edge=1.0,
            visc_edge=1.0,
        )


def test_cli_eta2y_reads_profile_table(tmp_path) -> None:
    """Check CLI reads a profile table and writes eta/y output."""

    # create profile table with Tecplot-style header and eta, u, tau columns
    profile_path = tmp_path / "profile.dat"
    profile_path.write_text(
        "TITLE = \"profile\"\n"
        "VARIABLES = \"eta\", \"u_ue\", \"T_Te\"\n"
        "ZONE T=\"solution\"\n"
        "0.0 0.0 1.0\n"
        "1.0 0.5 2.0\n"
        "2.0 1.0 3.0\n"
    )
    output_path = tmp_path / "y.dat"

    # run cli command using tau column from the Tecplot-style profile
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "eta2y",
            str(profile_path),
            "--x",
            "1.0",
            "--dens-edge",
            "2.0",
            "--uvel-edge",
            "1.0",
            "--visc-edge",
            "1.0",
            "--transform",
            "levy_lees",
            "--tau-column",
            "2",
            "--output",
            str(output_path),
        ],
    )

    # validate command completed and wrote expected output
    assert result.exit_code == 0, result.output
    output_table = np.loadtxt(output_path)
    np.testing.assert_allclose(output_table[:, 0], np.array([0.0, 1.0, 2.0]))
    np.testing.assert_allclose(output_table[:, 1], np.array([0.0, 1.5, 4.0]))
