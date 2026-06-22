"""Solve command for the simbl CLI."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from pathlib import Path
from typing import Annotated

try:
    import typer
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The simbl CLI requires typer. Install it with: pip install similarity-bl[cli]"
    ) from exc

# --------------------------------------------------
# default file names
# --------------------------------------------------
DEFAULT_CONFIG = "simbl_config.toml"


# --------------------------------------------------
# solve command: load config, apply overrides, run solver
#
# simbl solve takes an optional config file argument, and optional cli overrides for config values
# --------------------------------------------------
def cmd_solve(
    config: Annotated[
        Path | None,
        typer.Argument(help="TOML configuration file."),
    ] = None,
    # --------------------------------------------------
    # optional cli overrides (applied on top of config file values)
    # --------------------------------------------------
    mach: Annotated[
        float | None,
        typer.Option("--mach", "-m", help="Edge Mach number."),
    ] = None,
    temp_edge: Annotated[
        float | None,
        typer.Option("--temp-edge", help="Edge temperature [K]."),
    ] = None,
    beta: Annotated[
        float | None,
        typer.Option("--beta", "-b", help="Falkner-Skan pressure gradient parameter."),
    ] = None,
    wall: Annotated[
        str | None,
        typer.Option("--wall", "-w", help="Wall BC: adiabatic or isothermal."),
    ] = None,
    temp_wall: Annotated[
        float | None,
        typer.Option("--temp-wall", help="Wall temperature [K] (isothermal only)."),
    ] = None,
    equations: Annotated[
        str | None,
        typer.Option("--equations", help="Governing equations: falkner_skan or falkner_skan_cooke."),
    ] = None,
    eta_max: Annotated[
        float | None,
        typer.Option("--eta-max", help="Maximum similarity coordinate."),
    ] = None,
    n_points: Annotated[
        int | None,
        typer.Option("--n-points", help="Number of grid points."),
    ] = None,
    # --------------------------------------------------
    # output options
    # --------------------------------------------------
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output file (.dat or .json)."),
    ] = Path("simbl.dat"),
    quiet: Annotated[
        bool,
        typer.Option("--quiet", "-q", help="Suppress console output."),
    ] = False,
) -> None:
    """Solve the similarity boundary layer equations.

    Examples:

        simbl solve config.toml

        simbl solve config.toml --mach 6.0 --output result.dat

        simbl solve --mach 4.0 --wall adiabatic --output result.dat
    """

    # load necessary functions and classes for config handling and solving
    from simbl.config.config_ops import config_load, config_to_inputs
    from simbl.config.schema import SolverConfig
    from simbl.io import write
    from simbl.solver import solve_similarity

    # load config file
    # resolve which config file to use: explicit argument, then auto-discovered default
    if config is not None:
        config_path = Path(config)
    else:
        config_path = Path(DEFAULT_CONFIG) if Path(DEFAULT_CONFIG).exists() else None

    if config_path is None:
        typer.echo("No configuration file found.", err=True)
        typer.echo("  Run `simbl init` to generate a template, or pass a config file as an argument.", err=True)
        raise typer.Exit(1)

    if not config_path.exists():
        typer.echo(f"Config file not found: {config_path}", err=True)
        raise typer.Exit(1)

    try:
        cfg = config_load(config_path)
    except Exception as error:
        typer.echo(f"Error parsing {config_path}: {error}", err=True)
        raise typer.Exit(1) from None

    # apply cli overrides via mutable dict
    cfg_dict = cfg.model_dump()
    if mach is not None:
        cfg_dict["conditions"]["mach_edge"] = mach
    if temp_edge is not None:
        cfg_dict["conditions"]["temp_edge"] = temp_edge
    if beta is not None:
        cfg_dict["conditions"]["beta"] = beta
    if wall is not None:
        cfg_dict["wall"]["type"] = wall
    if temp_wall is not None:
        cfg_dict["wall"]["temp_wall"] = temp_wall
    if equations is not None:
        cfg_dict["equations"] = equations
    if eta_max is not None:
        cfg_dict["numerics"]["eta_max"] = eta_max
    if n_points is not None:
        cfg_dict["numerics"]["n_points"] = n_points

    # re-validate the config with the cli overrides applied, to catch any errors in the overrides
    try:
        cfg = SolverConfig(**cfg_dict)
    except Exception as error:
        typer.echo(f"Invalid configuration: {error}", err=True)
        raise typer.Exit(1) from None

    # convert validated config to solver inputs and run
    problem, options = config_to_inputs(cfg)

    # run the solver, catch any exceptions and print an error message
    try:
        sol, info = solve_similarity(problem, options)
    except Exception as error:
        typer.echo(f"Solver error: {error}", err=True)
        raise typer.Exit(1) from None

    # print summary to console
    if not quiet:
        typer.echo(f"Mach = {problem.mach_edge}, beta = {problem.beta}, wall = {problem.wall_bc}")
        typer.echo(f"  f''(0) = {sol.fpp[0]:.6f}")
        typer.echo(f"  tau(0)   = {sol.tau[0]:.6f}")
        typer.echo(f"  tau'(0)  = {sol.taup[0]:.6f}")
        # g'(0) is only present for Falkner-Skan-Cooke solutions
        if hasattr(sol, "gp"):
            typer.echo(f"  g'(0)  = {sol.gp[0]:.6f}")
        typer.echo(f"  Converged: {info.converged} ({info.iterations} iterations)")

    # write output file
    try:
        write(sol, output, problem=problem, shooting_result=info)
    except Exception as error:
        typer.echo(f"Error writing output: {error}", err=True)
        raise typer.Exit(1) from None
    if not quiet:
        typer.echo(f"  Output: {output}")



