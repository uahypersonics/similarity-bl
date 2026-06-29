"""CLI subcommand: simbl examples"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from pathlib import Path

import typer

from simbl.io.examples import EXAMPLES_REGISTRY, get_example_text


# --------------------------------------------------
# examples command: list or copy example config files
# --------------------------------------------------
def cmd_examples(
    name: str | None = typer.Argument(
        None,
        help="Example name to copy to current directory. Omit to list all examples.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Overwrite existing file.",
    ),
) -> None:
    """List available examples or copy one to the current directory.

    \b
    Examples:
        simbl examples                          # list all available examples
        simbl examples flat_plate_mach5_adiabatic  # copy to ./flat_plate_mach5_adiabatic.toml
    """

    # list all examples if no name given
    if name is None:
        typer.echo("Available examples:")
        for ex_name, (_, description) in EXAMPLES_REGISTRY.items():
            typer.echo(f"  {ex_name:<35}  {description}")
        return

    # validate name through the EXAMPLES_REGISTRY (defined in simbl.io.examples)
    if name not in EXAMPLES_REGISTRY:
        typer.echo(
            f"Error: unknown example {name!r}. Run 'simbl examples' to list available examples.",
            err=True,
        )
        raise typer.Exit(code=1)

    # copy TOML to current directory
    dest = Path(f"{name}.toml")
    if dest.exists() and not force:
        typer.echo(f"Error: {dest} already exists. Remove it first or use --force to overwrite.", err=True)
        raise typer.Exit(code=1)

    dest.write_text(get_example_text(name), encoding="utf-8")
    typer.echo(f"Written {dest}")
    typer.echo(f"Run with:  simbl solve {dest}")
