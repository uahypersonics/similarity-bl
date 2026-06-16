"""Map similarity coordinate eta to physical wall-normal coordinate y."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

from simbl.solver.equations import VALID_EQUATIONS

# --------------------------------------------------
# transform registries
# --------------------------------------------------
DEFAULT_ETA2Y_TRANSFORMS: dict[str, str] = {
    "falkner_skan": "levy_lees",
    "falkner_skan_cooke": "illingworth_stewartson",
}

Eta2YFunction = Callable[..., NDArray[np.float64]]


# --------------------------------------------------
#  eta2y dispatcher
# --------------------------------------------------
def eta2y(
    eta: NDArray[np.float64],
    tau: NDArray[np.float64],
    x: float,
    dens_edge: float,
    uvel_edge: float,
    visc_edge: float,
    *,
    beta: float = 0.0,
    dens_ref: float | None = None,
    visc_ref: float | None = None,
    transform: str | None = None,
    equations: str | None = None,
) -> NDArray[np.float64]:
    """Map similarity coordinate eta to physical wall-normal coordinate y.

    Args:
        eta: Similarity coordinate array.
        tau: Nondimensional temperature profile, tau = T / T_e.
        x: Streamwise station location in meters.
        dens_edge: Dimensional edge density.
        uvel_edge: Dimensional edge streamwise velocity.
        visc_edge: Dimensional edge dynamic viscosity.
        beta: Hartree pressure-gradient parameter, beta = 2m / (m + 1).
        dens_ref: Reference density for Illingworth-Stewartson. Defaults to ``dens_edge``.
        visc_ref: Reference dynamic viscosity for Illingworth-Stewartson. Defaults to ``visc_edge``.
        transform: Transform name. If omitted, the default is selected from ``equations``.
        equations: Governing equation family used when ``transform`` is omitted.

    Returns:
        Physical wall-normal coordinate array in meters.

    Raises:
        ValueError: If inputs are inconsistent or the transform cannot be resolved.
    """

    # validate the eta/tau coordinate pair before dispatching to the selected transform implementation
    if eta.ndim != 1:
        raise ValueError("eta must be a one-dimensional array")
    if tau.ndim != 1:
        raise ValueError("tau must be a one-dimensional array")
    if eta.shape != tau.shape:
        raise ValueError("eta and tau must have the same shape")
    if eta.size < 2:
        raise ValueError("eta and tau must contain at least two points")

    # resolve transform name from explicit name or equation-family default
    transform_name = _resolve_transform_name(transform=transform, equations=equations)

    # dispatch to concrete transform implementation
    transform_function = ETA2Y_TRANSFORMS[transform_name]

    # call transform function
    y = transform_function(
        eta=eta,
        tau=tau,
        x=x,
        dens_edge=dens_edge,
        uvel_edge=uvel_edge,
        visc_edge=visc_edge,
        beta=beta,
        dens_ref=dens_ref,
        visc_ref=visc_ref,
    )

    return y


# --------------------------------------------------
# eta 2 y transform for levy lees
# --------------------------------------------------
def eta2y_levy_lees(
    eta: NDArray[np.float64],
    tau: NDArray[np.float64],
    x: float,
    dens_edge: float,
    uvel_edge: float,
    visc_edge: float,
    beta: float = 0.0,
    dens_ref: float | None = None,
    visc_ref: float | None = None,
) -> NDArray[np.float64]:
    """Map eta to y using the flat-plate Levy-Lees inverse transform.

    Args:
        eta: Similarity coordinate array.
        tau: Nondimensional temperature profile, tau = T / T_e.
        x: Streamwise station location in meters.
        dens_edge: Dimensional edge density.
        uvel_edge: Dimensional edge streamwise velocity.
        visc_edge: Dimensional edge dynamic viscosity.
        beta: Hartree pressure-gradient parameter, beta = 2m / (m + 1).
        dens_ref: Unused. Accepted for dispatcher compatibility.
        visc_ref: Unused. Accepted for dispatcher compatibility.

    Returns:
        Physical wall-normal coordinate array in meters.
    """

    # compute Levy-Lees inverse scale from xi = rho_e * mu_e * u_e * x / (m + 1)
    m = beta / (2.0 - beta)
    eta_scale = np.sqrt((m + 1.0) * dens_edge * uvel_edge / (2.0 * visc_edge * x))

    # integrate inverse coordinate map
    y = _integrate_tau_over_eta(
        eta=eta,
        tau=tau,
        eta_scale=eta_scale,
    )

    return y


# --------------------------------------------------
# eta 2 y transform for illingworth stewartson
# --------------------------------------------------
def eta2y_illingworth_stewartson(
    eta: NDArray[np.float64],
    tau: NDArray[np.float64],
    x: float,
    dens_edge: float,
    uvel_edge: float,
    visc_edge: float,
    beta: float = 0.0,
    dens_ref: float | None = None,
    visc_ref: float | None = None,
) -> NDArray[np.float64]:
    """Map eta to y using the Illingworth-Stewartson inverse transform.

    Args:
        eta: Similarity coordinate array.
        tau: Nondimensional temperature profile, tau = T / T_e.
        x: Streamwise station location in meters.
        dens_edge: Dimensional edge density.
        uvel_edge: Dimensional edge streamwise velocity.
        visc_edge: Dimensional edge dynamic viscosity.
        beta: Hartree pressure-gradient parameter, beta = 2m / (m + 1).
        dens_ref: Reference density. Defaults to ``dens_edge`` for local scaling.
        visc_ref: Reference dynamic viscosity. Defaults to ``visc_edge`` for local scaling.

    Returns:
        Physical wall-normal coordinate array in meters.
    """

    # default reference state to local edge state for locally self-similar usage
    if dens_ref is None:
        dens_ref = dens_edge
    if visc_ref is None:
        visc_ref = visc_edge

    # compute Illingworth-Stewartson inverse scale from eta = bar_eta * sqrt((m + 1) U_e / (2 nu_ref x))
    m = beta / (2.0 - beta)
    nu_ref = visc_ref / dens_ref
    eta_scale = np.sqrt((m + 1.0) * uvel_edge / (2.0 * nu_ref * x))

    # integrate inverse coordinate map
    y = _integrate_tau_over_eta(
        eta=eta,
        tau=tau,
        eta_scale=eta_scale,
    )

    return y


ETA2Y_TRANSFORMS: dict[str, Eta2YFunction] = {
    "levy_lees": eta2y_levy_lees,
    "illingworth_stewartson": eta2y_illingworth_stewartson,
}


# --------------------------------------------------
# function to check which transform should be used
# --------------------------------------------------
def _resolve_transform_name(transform: str | None, equations: str | None) -> str:
    """Resolve explicit transform names."""

    # resolve omitted transform from equation family
    if transform is None:
        if equations is None:
            raise ValueError("equations must be provided when transform is omitted")

        # normalize the equation name, strip whitespace, force to all lowercase and replace hyphens with underscores for flexible matching
        equations_key = equations.strip().lower().replace("-", "_")

        # validate the equation name and look up the default transform for that equation family
        if equations_key not in VALID_EQUATIONS:
            valid_equations = ", ".join(sorted(VALID_EQUATIONS))
            raise ValueError(f"Unknown equations '{equations}'. Choose: {valid_equations}")

        # look up the default transform for the validated equation family
        transform_key = DEFAULT_ETA2Y_TRANSFORMS[equations_key]

    else:
        # use the explicitly provided transform name
        transform_key = transform.strip().lower().replace("-", "_")

        # sanity check: validate the explicit transform name is known
        if transform_key not in ETA2Y_TRANSFORMS:
            valid_names = ", ".join(sorted(ETA2Y_TRANSFORMS))
            raise ValueError(f"Unknown eta2y transform '{transform}'. Choose: {valid_names}")

    return transform_key


def _integrate_tau_over_eta(
    eta: NDArray[np.float64],
    tau: NDArray[np.float64],
    eta_scale: float,
) -> NDArray[np.float64]:
    """Integrate dy/deta = tau / eta_scale."""

    # integrate dy/deta = tau / eta_scale using a cumulative trapezoid rule
    delta_eta = np.diff(eta)
    average_tau = 0.5 * (tau[1:] + tau[:-1])
    segment_integrals = average_tau * delta_eta
    cumulative_integral = np.concatenate(([0.0], np.cumsum(segment_integrals)))
    y = cumulative_integral / eta_scale

    return y
