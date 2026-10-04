"""Tests for initializing SIMBL configuration from a FlowState."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import pytest
from flow_state import from_mach_pres_temp, solve
from flow_state.gas import PerfectGas
from flow_state.io import write_json
from flow_state.transport import PowerLaw
from typer.testing import CliRunner

from simbl.cli import cli
from simbl.config import config_init, config_load


# --------------------------------------------------
# FlowState config initialization tests
# --------------------------------------------------
def test_config_init_imports_complete_flow_state(tmp_path) -> None:
    """A complete FlowState should populate all matching SIMBL properties."""

    # build and write a complete edge state
    flow_path = tmp_path / "edge_conditions.json"
    state = solve(
        mach=5.4,
        pres=3200.0,
        temp=187.5,
        pr=0.71,
    )
    write_json(state, flow_path)

    # initialize and reload the SIMBL config
    config_path = tmp_path / "simbl_config.toml"
    warnings = config_init(config_path, flow_state=flow_path)
    config = config_load(config_path)

    # validate imported physics and source provenance
    config_text = config_path.read_text(encoding="utf-8")
    assert warnings == []
    assert config.conditions.mach_edge == pytest.approx(state.mach)
    assert config.conditions.temp_edge == pytest.approx(state.temp)
    assert config.gas.gamma == pytest.approx(state.gamma)
    assert config.gas.prandtl == pytest.approx(state.pr)
    assert config.viscosity.model == state.transport_model.model_type
    assert config.viscosity.mu_ref == pytest.approx(state.transport_model.parameters["mu_ref"])
    assert config.viscosity.T_ref == pytest.approx(state.transport_model.parameters["T_ref"])
    assert config.viscosity.S == pytest.approx(state.transport_model.parameters["S"])
    assert f"# initialized from FlowState: {flow_path}" in config_text


def test_config_init_warns_for_missing_optional_properties(tmp_path) -> None:
    """Missing Prandtl and transport data should retain defaults with warnings."""

    # build a valid state without optional transport properties
    flow_path = tmp_path / "edge_conditions.json"
    state = from_mach_pres_temp(
        mach=4.0,
        pres=2500.0,
        temp=210.0,
        gas=PerfectGas.air(),
        transport=None,
        pr=None,
    )
    write_json(state, flow_path)

    # initialize and reload the SIMBL config
    config_path = tmp_path / "simbl_config.toml"
    warnings = config_init(config_path, flow_state=flow_path)
    config = config_load(config_path)

    # validate explicit warnings and retained defaults
    assert len(warnings) == 2
    assert "Prandtl number" in warnings[0]
    assert "transport model" in warnings[1]
    assert config.gas.prandtl == pytest.approx(0.72)
    assert config.viscosity.model == "sutherland"


def test_config_init_warns_for_unrepresentable_transport_parameters(tmp_path) -> None:
    """Unsupported transport parameters should use the registered model default."""

    # build a state whose power-law exponent is not representable in the schema
    flow_path = tmp_path / "edge_conditions.json"
    state = from_mach_pres_temp(
        mach=4.0,
        pres=2500.0,
        temp=210.0,
        gas=PerfectGas.air(),
        transport=PowerLaw.air(m=0.81),
        pr=0.71,
    )
    write_json(state, flow_path)

    # initialize and reload the SIMBL config
    config_path = tmp_path / "simbl_config.toml"
    warnings = config_init(config_path, flow_state=flow_path)
    config = config_load(config_path)

    # validate model selection and explicit loss-of-parameters warning
    assert len(warnings) == 1
    assert "m" in warnings[0]
    assert "registered default" in warnings[0]
    assert config.viscosity.model == "power_law"
    assert config.viscosity.mu_ref is None
    assert config.viscosity.T_ref is None


def test_cli_init_displays_flow_state_warnings(tmp_path) -> None:
    """CLI should display consistency warnings returned by config initialization."""

    # build a state without optional properties
    flow_path = tmp_path / "edge_conditions.json"
    state = from_mach_pres_temp(
        mach=4.0,
        pres=2500.0,
        temp=210.0,
        gas=PerfectGas.air(),
        transport=None,
        pr=None,
    )
    write_json(state, flow_path)

    # initialize through the public CLI alias
    config_path = tmp_path / "simbl_config.toml"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["init", "--flow", str(flow_path), "--output", str(config_path)],
    )

    # validate warning visibility and successful output
    assert result.exit_code == 0
    assert "Warning: FlowState does not contain a Prandtl number" in result.output
    assert "Warning: FlowState does not contain a transport model" in result.output
    assert config_path.is_file()
