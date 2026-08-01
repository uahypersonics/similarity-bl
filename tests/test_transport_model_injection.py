"""Tests for viscosity-model injection at the similarity solver boundary."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from types import SimpleNamespace

from flow_state.transport import get_transport_model

from simbl import SimilarityInputs, solve_similarity
from simbl.solver import main as solver_main
from simbl.solver.falkner_skan import builder, solution


# --------------------------------------------------
# viscosity-model injection tests
# --------------------------------------------------
def test_solve_similarity_uses_injected_viscosity_model(monkeypatch) -> None:
    """An injected viscosity model should bypass string-based construction."""

    # build an existing flow-state model and a minimal similarity problem
    visc_model = get_transport_model("keyes")
    inputs = SimilarityInputs(mach_edge=4.0, temp_edge=300.0)
    captured_models = []

    # replace numerical work so only transport-model routing is exercised
    def _known_problem(problem, received_model, **kwargs):
        captured_models.append(received_model)
        return object()

    expected_result = SimpleNamespace(converged=True)
    expected_solution = SimpleNamespace(fpp=[0.5])

    def _unexpected_factory(*args, **kwargs):
        raise AssertionError("string factory should not be called")

    monkeypatch.setattr(builder, "build_solver_problem", _known_problem)
    monkeypatch.setattr(
        solver_main,
        "shooting_method",
        lambda solver_problem, options: expected_result,
    )
    monkeypatch.setattr(solution, "build_solution", lambda result: expected_solution)
    monkeypatch.setattr(solver_main, "get_transport_model", _unexpected_factory)

    # solve with the existing model object
    actual_solution, actual_result = solve_similarity(inputs, visc_model=visc_model)

    # check
    assert captured_models == [visc_model]
    assert actual_solution is expected_solution
    assert actual_result is expected_result


def test_solve_similarity_constructs_configured_viscosity_model(monkeypatch) -> None:
    """String-based callers should retain the existing model factory path."""

    # build the configured model returned by the flow-state registry
    expected_model = get_transport_model("keyes")
    inputs = SimilarityInputs(
        mach_edge=4.0,
        temp_edge=300.0,
        viscosity_model="keyes",
    )
    factory_calls = []
    captured_models = []

    # replace model construction and numerical work to isolate routing
    def _known_factory(name, **kwargs):
        factory_calls.append((name, kwargs))
        return expected_model

    def _known_problem(problem, received_model, **kwargs):
        captured_models.append(received_model)
        return object()

    expected_result = SimpleNamespace(converged=True)
    expected_solution = SimpleNamespace(fpp=[0.5])

    monkeypatch.setattr(solver_main, "get_transport_model", _known_factory)
    monkeypatch.setattr(builder, "build_solver_problem", _known_problem)
    monkeypatch.setattr(
        solver_main,
        "shooting_method",
        lambda solver_problem, options: expected_result,
    )
    monkeypatch.setattr(solution, "build_solution", lambda result: expected_solution)

    # solve through the existing string configuration path
    actual_solution, actual_result = solve_similarity(inputs)

    # check
    assert factory_calls == [("keyes", {})]
    assert captured_models == [expected_model]
    assert actual_solution is expected_solution
    assert actual_result is expected_result
