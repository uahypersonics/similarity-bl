"""Eta-to-physical-y command for the simbl CLI."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

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
# eta2y command: map similarity coordinate to physical wall-normal coordinate
# --------------------------------------------------
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
    except Exception as error:
        typer.echo(f"Error reading profile table: {error}", err=True)
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
    except ValueError as error:
        typer.echo(f"Error mapping eta to y: {error}", err=True)
        raise typer.Exit(1) from None

    # write eta/y table to file or stdout
    try:
        _write_eta2y_output(output=output, eta=eta, y=y)
    except Exception as error:
        typer.echo(f"Error writing eta2y output: {error}", err=True)
        raise typer.Exit(1) from None

    # report output path when a file was written
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
    """Read eta and tau columns from a whitespace- or comma-delimited table.

    Args:
        profile: Input table path.
        eta_column: Zero-based eta column index.
        tau_column: Zero-based tau column index.

    Returns:
        Eta and tau arrays.

    Raises:
        FileNotFoundError: If the input profile does not exist.
        ValueError: If column indices are invalid or no numeric profile is found.
    """

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

    # convert lists to numpy arrays
    eta = np.asarray(eta_values, dtype=np.float64)
    tau = np.asarray(tau_values, dtype=np.float64)

    return eta, tau


def _write_eta2y_output(output: Path | None, eta: np.ndarray, y: np.ndarray) -> None:
    """Write eta/y table to a file or stdout.

    Args:
        output: Output path, or None to print to stdout.
        eta: Similarity coordinate array.
        y: Physical wall-normal coordinate array.
    """

    # build output table
    table = np.column_stack((eta, y))

    # write to stdout when no output path is requested
    if output is None:
        np.savetxt(sys.stdout, table, header="eta y", comments="# ")
        return

    # write to requested output path
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, table, header="eta y", comments="# ")
