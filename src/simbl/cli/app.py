"""Typer application setup for the simbl command-line interface."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import logging
import sys
from typing import Annotated

try:
    import typer
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The simbl CLI requires typer. Install it with: pip install similarity-bl[cli]"
    ) from exc

# --------------------------------------------------
# package imports
# --------------------------------------------------
from simbl.cli.eta2y import cmd_eta2y
from simbl.cli.init import cmd_init
from simbl.cli.solve import cmd_solve

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
    """Print version and exit."""

    # print version and exit when requested
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
    """simbl: Compressible Similarity Boundary Layer Solver."""

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
            isinstance(handler, logging.StreamHandler)
            and not isinstance(handler, logging.FileHandler)
            for handler in simbl_logger.handlers
        )
        if not has_stream_handler:
            handler = logging.StreamHandler(sys.stderr)
            handler.setLevel(logging.INFO)
            handler.setFormatter(logging.Formatter("[%(levelname)-7s] %(name)s: %(message)s"))
            simbl_logger.addHandler(handler)


# --------------------------------------------------
# command registration
# --------------------------------------------------
cli.command("init")(cmd_init)
cli.command("solve")(cmd_solve)
cli.command("eta2y")(cmd_eta2y)


# --------------------------------------------------
# entry point
# --------------------------------------------------
if __name__ == "__main__":
    cli()
