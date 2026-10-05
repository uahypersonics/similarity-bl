"""Physical-grid baseflow workflow for the simbl CLI."""

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

from simbl.baseflow import (
    BASEFLOW_CONFIG_TEMPLATE,
    load_baseflow_config,
    map_baseflow_to_grid,
)

# --------------------------------------------------
# baseflow command group
# --------------------------------------------------
baseflow_app = typer.Typer(
    name="baseflow",
    help="Map converged similarity profiles onto structured physical grids.",
    no_args_is_help=True,
)

DEFAULT_BASEFLOW_CONFIG = Path("simbl_baseflow.toml")


# --------------------------------------------------
# initialize baseflow config
# --------------------------------------------------
@baseflow_app.command("init")
def cmd_baseflow_init(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output config file path."),
    ] = DEFAULT_BASEFLOW_CONFIG,
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite an existing config."),
    ] = False,
) -> None:
    """Write a focused physical-grid baseflow config."""

    # prevent accidental config replacement
    if output.exists() and not force:
        typer.echo(f"File already exists: {output}", err=True)
        typer.echo("Use --force to overwrite.", err=True)
        raise typer.Exit(1)

    # write and validate the generated config
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(BASEFLOW_CONFIG_TEMPLATE, encoding="utf-8")
    load_baseflow_config(output)

    typer.echo(f"Created: {output}")
    typer.echo(f"Edit the file, then run: simbl baseflow run {output}")


# --------------------------------------------------
# run baseflow mapping
# --------------------------------------------------
@baseflow_app.command("run")
def cmd_baseflow_run(
    config: Annotated[
        Path,
        typer.Argument(help="Physical-grid baseflow TOML configuration."),
    ] = DEFAULT_BASEFLOW_CONFIG,
) -> None:
    """Map a SIMBL JSON solution onto a structured HDF5 grid."""

    # load, validate, and execute the focused workflow
    try:
        baseflow_config = load_baseflow_config(config)
        output_path = map_baseflow_to_grid(baseflow_config)
    except (ImportError, NotImplementedError, OSError, TypeError, ValueError) as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(1) from None

    typer.echo(f"Wrote: {output_path}")
