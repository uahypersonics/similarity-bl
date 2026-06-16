"""Command-line interface for simbl.

Commands:
    simbl init   - generate a template config file
    simbl solve  - solve similarity BL equations from a config file
    simbl eta2y  - map eta profiles to physical y coordinates
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Annotated

import numpy as np

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
DEFAULT_OUTPUT = "solution.dat"

# --------------------------------------------------
# set up cli using typer
#
# typer builds the CLI from function signatures using Annotated:
#
#   param: Annotated[type, typer.Option("--flag", help="...")] = default
#
#   - type        : expected Python type; typer coerces input and rejects invalid values
#   - typer.Option: named flag (--flag value); typer.Argument is positional (no --)
#   - default     : value used if the flag is not passed on the command line
#
# the function docstring becomes the --help text
# the @cli.command("name") decorator registers the function as a subcommand
# --------------------------------------------------
cli = typer.Typer(
    name="simbl",
    help="SIMBL - Compressible Similarity Boundary Layer Solver",
    no_args_is_help=True,
    add_completion=False,
)


# --------------------------------------------------
# version callback
# user can run `simbl --version` or `simbl -V` to see the version number
# nothing is done and the program exits immediately after printing the version.
# --------------------------------------------------
def _version_callback(value: bool) -> None:
    """Print version and exit"""
    if value:
        from simbl import __version__

        typer.echo(f"simbl {__version__}")
        raise typer.Exit()


# --------------------------------------------------
# main callback: runs when no subcommand is given, just prints help
# --------------------------------------------------
@cli.callback()
def main(
    ctx: typer.Context,
    version: Annotated[
        bool | None,
        typer.Option(
            "--version",
            "-V",
            help="Show version and exit.",
            callback=_version_callback,
            is_eager=True,
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Enable informational output."),
    ] = False,
) -> None:
    """simbl: Compressible Similarity Boundary Layer Solver"""

    # store verbose flag in context for access by subcommands
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose

    # configure logging if verbose is requested
    if verbose:
        simbl_logger = logging.getLogger("simbl")
        simbl_logger.setLevel(logging.INFO)

        # attach a handler that prints to stderr
        # guard prevents adding duplicate handlers if callback runs more than once
        has_stream_handler = any(
            isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
            for h in simbl_logger.handlers
        )
        if not has_stream_handler:
            handler = logging.StreamHandler(sys.stderr)
            handler.setLevel(logging.INFO)
            handler.setFormatter(logging.Formatter("[%(levelname)-7s] %(name)s: %(message)s"))
            simbl_logger.addHandler(handler)


# --------------------------------------------------
# init command: write a blank template config file
#
# simbl init works without any additional options
# it defaults to writing simbl_config.toml in the current directory
#
# the user can specify:
# - a different output path with --output
# - a different template type with --equations (fs or fsc)
# - force overwrite an existing file with --force
# --------------------------------------------------

# decorator to register the command with typer, and define its options and arguments
@cli.command("init")
def cmd_init(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output config file path."),
    ] = Path(DEFAULT_CONFIG),
    equations: Annotated[
        str,
        typer.Option("--equations", help="Template type: fs (Falkner-Skan) or fsc (Falkner-Skan-Cooke)."),
    ] = "fs",
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing file."),
    ] = False,
) -> None:
    """Generate a template configuration file"""

    # load the config_ops.config_init function to access the config templates and write the file
    from simbl.config import config_init

    try:
        config_init(output, equations=equations, force=force)
    except FileExistsError as e:
        typer.echo(str(e), err=True)
        typer.echo("Use --force to overwrite.", err=True)
        raise typer.Exit(1) from None
    except ValueError as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(1) from None

    typer.echo(f"Created: {output}")
    typer.echo(f"Edit the file, then run: simbl solve {output}")


# --------------------------------------------------
# solve command: load config, apply overrides, run solver
#
# simbl solve takes an optional config file argument, and optional CLI overrides for config values
# --------------------------------------------------

# decorator to register the command with typer, and define its options and arguments
@cli.command("solve")
def cmd_solve(
    config: Annotated[
        Path | None,
        typer.Argument(help="TOML configuration file."),
    ] = None,
    # --------------------------------------------------
    # optional CLI overrides (applied on top of config file values)
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
        Path | None,
        typer.Option("--output", "-o", help="Output file (.dat or .json)."),
    ] = None,
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
    from simbl.config.schema import ConditionsConfig, SolverConfig
    from simbl.io import write
    from simbl.solver import solve_similarity

    # --------------------------------------------------
    # load config file or build a minimal default
    # --------------------------------------------------
    if config is not None:

        # config file given: attempt to load it
        if not config.exists():
            # config file path given but not found -> print error and exit
            typer.echo(f"Config file not found: {config}", err=True)
            raise typer.Exit(1)
        try:
            # load the config file, parse it into a SolverConfig dataclass, and validate it
            cfg = config_load(config)
        except Exception as e:
            # config file found but error parsing -> print error and exit
            typer.echo(f"Error parsing {config}: {e}", err=True)
            raise typer.Exit(1) from None
    else:

        # no config file given: check for DEFAULT_CONFIG in the current directory
        default_path = Path(DEFAULT_CONFIG)

        if default_path.exists():
            # found simbl_config.toml in cwd -> load it silently (natural workflow after simbl init)
            try:
                # load the config file, parse it into a SolverConfig dataclass, and validate it
                cfg = config_load(default_path)
            except Exception as e:
                # config file found but error parsing -> print error and exit
                typer.echo(f"Error parsing {default_path}: {e}", err=True)
                raise typer.Exit(1) from None
        elif mach is not None and temp_edge is not None:
            # no config file given, no default config found, but user supplied minimum required conditions via CLI flags -> build config from scratch
            cfg = SolverConfig(conditions=ConditionsConfig(mach_edge=mach, temp_edge=temp_edge))
        else:
            # no config file (provided or default) and insufficient CLI flags -> cannot proceed
            typer.echo("No configuration found.", err=True)
            typer.echo(
                "  Run `simbl init` to generate a template, or pass --mach and --temp-edge.",
                err=True,
            )
            raise typer.Exit(1)

    # --------------------------------------------------
    # apply CLI overrides via mutable dict
    #
    # user overrides take precedence over config file values, which take precedence over defaults
    # --------------------------------------------------
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

    # --------------------------------------------------
    # re-validate the config with the CLI overrides applied, to catch any errors in the overrides
    # --------------------------------------------------
    try:
        cfg = SolverConfig(**cfg_dict)
    except Exception as e:
        typer.echo(f"Invalid configuration: {e}", err=True)
        raise typer.Exit(1) from None

    # --------------------------------------------------
    # convert validated config to solver inputs and run
    #
    # wired up in config layer (config_ops.py)
    # --------------------------------------------------
    problem, options = config_to_inputs(cfg)

    # --------------------------------------------------
    # run the solver, catch any exceptions and print an error message
    # --------------------------------------------------
    try:
        sol, info = solve_similarity(problem, options)
    except Exception as e:
        typer.echo(f"Solver error: {e}", err=True)
        raise typer.Exit(1) from None

    # --------------------------------------------------
    # print summary to console
    # --------------------------------------------------
    if not quiet:
        typer.echo(f"Mach = {problem.mach_edge}, beta = {problem.beta}, wall = {problem.wall_bc}")
        typer.echo(f"  f''(0) = {sol.fpp[0]:.6f}")
        typer.echo(f"  g(0)   = {sol.g[0]:.6f}")
        typer.echo(f"  g'(0)  = {sol.gp[0]:.6f}")
        # w'(0) is only present for Falkner-Skan-Cooke (swept-wing) solutions
        if hasattr(sol, "wp"):
            typer.echo(f"  w'(0)  = {sol.wp[0]:.6f}")
        typer.echo(f"  Converged: {info.converged} ({info.iterations} iterations)")

    # --------------------------------------------------
    # write output file if requested
    # --------------------------------------------------
    if output is not None:
        try:
            write(sol, output, problem=problem, shooting_result=info)
        except Exception as e:
            typer.echo(f"Error writing output: {e}", err=True)
            raise typer.Exit(1) from None
        if not quiet:
            typer.echo(f"  Output: {output}")


# --------------------------------------------------
# eta2y command: map similarity coordinate to physical wall-normal coordinate
# --------------------------------------------------
@cli.command("eta2y")
def cmd_eta2y(
    profile: Annotated[
        Path,
        typer.Argument(help="Profile table containing eta and tau columns."),
    ],
    x: Annotated[
        float,
        typer.Option("--x", help="Streamwise station location in meters."),
    ],
    dens_edge: Annotated[
        float,
        typer.Option("--dens-edge", help="Dimensional edge density."),
    ],
    uvel_edge: Annotated[
        float,
        typer.Option("--uvel-edge", help="Dimensional edge streamwise velocity."),
    ],
    visc_edge: Annotated[
        float,
        typer.Option("--visc-edge", help="Dimensional edge dynamic viscosity."),
    ],
    beta: Annotated[
        float,
        typer.Option("--beta", help="Hartree pressure-gradient parameter."),
    ] = 0.0,
    dens_ref: Annotated[
        float | None,
        typer.Option("--dens-ref", help="Reference density for Illingworth-Stewartson."),
    ] = None,
    visc_ref: Annotated[
        float | None,
        typer.Option("--visc-ref", help="Reference dynamic viscosity for Illingworth-Stewartson."),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output table path. Prints to stdout if omitted."),
    ] = None,
    equations: Annotated[
        str | None,
        typer.Option("--equations", help="Equation family used when --transform is omitted."),
    ] = None,
    transform: Annotated[
        str | None,
        typer.Option("--transform", help="Transform name: levy_lees or illingworth_stewartson."),
    ] = None,
    eta_column: Annotated[
        int,
        typer.Option("--eta-column", help="Zero-based eta column index."),
    ] = 0,
    tau_column: Annotated[
        int,
        typer.Option("--tau-column", help="Zero-based tau column index."),
    ] = 1,
) -> None:
    """Map eta to physical y using a similarity coordinate transform."""

    # load transform helper only when the command runs
    from simbl.transform import eta2y

    # read eta and tau from the profile table
    try:
        eta, tau = _read_eta_tau_columns(
            profile=profile,
            eta_column=eta_column,
            tau_column=tau_column,
        )
    except Exception as e:
        typer.echo(f"Error reading profile table: {e}", err=True)
        raise typer.Exit(1) from None

    # compute physical wall-normal coordinate
    try:
        y = eta2y(
            eta=eta,
            tau=tau,
            x=x,
            dens_edge=dens_edge,
            uvel_edge=uvel_edge,
            visc_edge=visc_edge,
            beta=beta,
            dens_ref=dens_ref,
            visc_ref=visc_ref,
            transform=transform,
            equations=equations,
        )
    except ValueError as e:
        typer.echo(f"Error mapping eta to y: {e}", err=True)
        raise typer.Exit(1) from None

    # write eta/y table to file or stdout
    try:
        _write_eta2y_output(output=output, eta=eta, y=y)
    except Exception as e:
        typer.echo(f"Error writing eta2y output: {e}", err=True)
        raise typer.Exit(1) from None

    if output is not None:
        typer.echo(f"Wrote: {output}")


# --------------------------------------------------
# eta2y command helpers
# --------------------------------------------------
def _read_eta_tau_columns(
    profile: Path,
    eta_column: int,
    tau_column: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Read eta and tau columns from a whitespace- or comma-delimited table."""

    # validate inputs
    if not profile.exists():
        raise FileNotFoundError(profile)
    if eta_column < 0:
        raise ValueError("eta_column must be nonnegative")
    if tau_column < 0:
        raise ValueError("tau_column must be nonnegative")

    # initialize output lists
    eta_values: list[float] = []
    tau_values: list[float] = []

    # read table line by line so Tecplot-style headers can be ignored
    for line in profile.read_text().splitlines():
        stripped_line = line.strip()

        # skip blank and comment lines
        if not stripped_line or stripped_line.startswith("#"):
            continue

        # skip common Tecplot header lines
        lower_line = stripped_line.lower()
        if lower_line.startswith(("title", "variables", "zone", "auxdata")):
            continue

        # split on whitespace after normalizing comma-delimited rows
        normalized_line = stripped_line.replace(",", " ")
        columns = normalized_line.split()

        # skip nonnumeric header rows that are not caught above
        try:
            eta_value = float(columns[eta_column])
            tau_value = float(columns[tau_column])
        except (IndexError, ValueError):
            continue

        # store numeric row values
        eta_values.append(eta_value)
        tau_values.append(tau_value)

    # validate enough data were read to form a profile
    if len(eta_values) < 2:
        raise ValueError("profile table must contain at least two numeric rows")

    eta = np.asarray(eta_values, dtype=np.float64)
    tau = np.asarray(tau_values, dtype=np.float64)

    return eta, tau


def _write_eta2y_output(output: Path | None, eta: np.ndarray, y: np.ndarray) -> None:
    """Write eta/y table to a file or stdout."""

    # build output table
    table = np.column_stack((eta, y))

    # write to stdout when no output path is requested
    if output is None:
        np.savetxt(sys.stdout, table, header="eta y", comments="# ")
        return

    # write to requested output path
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, table, header="eta y", comments="# ")


# --------------------------------------------------
# entry point
# --------------------------------------------------
if __name__ == "__main__":
    cli()
