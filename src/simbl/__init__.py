"""
SIMBL: Compressible Similarity Boundary Layer Solver.

A modern Python tool for solving compressible boundary layer similarity
equations using the Lees-Dorodnitsyn transformation.
"""

# --------------------------------------------------
# stdlib imports
# --------------------------------------------------
import logging
from importlib.metadata import PackageNotFoundError, version

# --------------------------------------------------
# package imports
# --------------------------------------------------
from simbl.config import SolverConfig, config_init, config_load, config_save, config_to_inputs
from simbl.io import write
from simbl.solver import (
    FalknerSkanCookeSolution,
    FalknerSkanSolution,
    ShootingResult,
    solve_similarity,
)
from simbl.solver.bcs import WallBCType
from simbl.solver.inputs import SimilarityInputs
from simbl.solver.options import SolverOptions, default_options
from simbl.transform import eta2y

# --------------------------------------------------
# configure package-level logger
# NullHandler suppresses "No handler found" warnings when the user's
# application has not configured logging
# --------------------------------------------------
logging.getLogger("simbl").addHandler(logging.NullHandler())

# --------------------------------------------------
# load version
# --------------------------------------------------
try:
    __version__ = version("similarity-bl")
except PackageNotFoundError:
    __version__ = "unknown"

# --------------------------------------------------
# public api
# --------------------------------------------------

__all__ = [
    "FalknerSkanCookeSolution",
    "FalknerSkanSolution",
    "ShootingResult",
    "SimilarityInputs",
    "SolverConfig",
    "SolverOptions",
    "default_options",
    "WallBCType",
    "__version__",
    "config_init",
    "config_load",
    "config_save",
    "config_to_inputs",
    "eta2y",
    "solve_similarity",
    "write",
]

