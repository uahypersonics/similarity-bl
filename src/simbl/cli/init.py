"""Configuration-template command for the simbl CLI."""

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
def cmd_init(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output config file path."),
    ] = Path(DEFAULT_CONFIG),
    equations: Annotated[
        str,
        typer.Option(
            "--equations",
            help="Template type: fs (Falkner-Skan) or fsc (Falkner-Skan-Cooke).",
        ),
    ] = "fs",
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing file."),
    ] = False,
) -> None:
    """Generate a template configuration file."""

    # load the config_ops.config_init function to access the config templates and write the file
    from simbl.config import config_init

    # write template config and report validation errors
    try:
        config_init(output, equations=equations, force=force)
    except FileExistsError as error:
        typer.echo(str(error), err=True)
        typer.echo("Use --force to overwrite.", err=True)
        raise typer.Exit(1) from None
    except ValueError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from None

    # report created file and next command
    typer.echo(f"Created: {output}")
    typer.echo(f"Edit the file, then run: simbl solve {output}")
