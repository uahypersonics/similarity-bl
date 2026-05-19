"""Similarity boundary layer solver package

Public API re-exported here for convenience:
- solve_similarity: main entry point
- SolverOptions: numerical settings
- FalknerSkanSolution / FalknerSkanCookeSolution: solution dataclasses
- ShootingResult: convergence info
"""

from simbl.solver.falkner_skan.solution import FalknerSkanSolution
from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
from simbl.solver.main import solve_similarity
from simbl.solver.options import SolverOptions
from simbl.solver.shooting import ShootingResult

__all__ = [
    "FalknerSkanCookeSolution",
    "FalknerSkanSolution",
    "ShootingResult",
    "SolverOptions",
    "solve_similarity",
]
