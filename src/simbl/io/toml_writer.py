"""Raw TOML file writer.

Thin wrapper around tomli_w. Accepts a plain dict — no schema knowledge here.
Callers in config/ are responsible for serializing objects to dict first.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from pathlib import Path

import tomli_w


# --------------------------------------------------
# _write_toml: write a plain dict to a TOML file
# --------------------------------------------------
def _write_toml(data: dict, fname: Path) -> None:
    """Write a plain dict to a TOML file

    Parameters
    ----------
    data : dict
        Data to serialize. Must be TOML-compatible (no None values, etc.).
    fname : Path
        Output file path. Parent directory must exist.
    """
    with open(fname, "wb") as f:
        tomli_w.dump(data, f)
