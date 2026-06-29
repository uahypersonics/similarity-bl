"""Registry of built-in example configurations for simbl."""

from __future__ import annotations

from importlib.resources import files

# -- registry: name -> (filename, one-line description) --
EXAMPLES_REGISTRY: dict[str, tuple[str, str]] = {
    "blasius": (
        "blasius.toml",
        "Blasius flat-plate boundary layer (near-incompressible, M=0.1)",
    ),
    "flat_plate_mach5_adiabatic": (
        "flat_plate_mach5_adiabatic.toml",
        "Mach 5 adiabatic flat plate (T_edge=68 K)",
    ),
    "flat_plate_mach5_isothermal": (
        "flat_plate_mach5_isothermal.toml",
        "Mach 5 cold-wall flat plate (T_wall=300 K, T_edge=68 K)",
    ),
    "wedge_mach3": (
        "wedge_mach3.toml",
        "Mach 3 wedge - favourable pressure gradient (beta=0.33)",
    ),
    "swept_plate_mach6": (
        "swept_plate_mach6.toml",
        "Mach 6 swept flat plate, 45 deg sweep (Falkner-Skan-Cooke)",
    ),
}


def available_example_names() -> list[str]:
    """Return list of available example names."""
    return list(EXAMPLES_REGISTRY.keys())


def get_example_text(name: str) -> str:
    """Return TOML text for the named example.

    Parameters
    ----------
    name : str
        Example name (key in EXAMPLES_REGISTRY).

    Returns
    -------
    str
        TOML file contents as a string.

    Raises
    ------
    KeyError
        If name is not in the registry.
    """
    if name not in EXAMPLES_REGISTRY:
        raise KeyError(f"Unknown example {name!r}. Available: {available_example_names()}")
    fname, _ = EXAMPLES_REGISTRY[name]
    return (files("simbl") / "examples" / fname).read_text(encoding="utf-8")
