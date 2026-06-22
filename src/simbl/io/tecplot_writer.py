"""Tecplot ASCII writer for similarity solutions

Writes dimensionless physical profiles to Tecplot ASCII format:
  eta, u/u_e, T/T_e, rho/rho_e, (v/u_e)*sqrt(Re_x)
Crossflow velocity w/w_e is included when present (FSC model).
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

# imports only used in type annotations, not at runtime (improved performance)
# avoids circular imports between solver modules
if TYPE_CHECKING:
    from simbl.solver.falkner_skan.solution import FalknerSkanSolution
    from simbl.solver.falkner_skan_cooke.solution import FalknerSkanCookeSolution
    from simbl.solver.inputs import SimilarityInputs
    from simbl.solver.shooting import ShootingResult


# --------------------------------------------------
# write_tecplot: write similarity profiles to Tecplot ASCII format
# --------------------------------------------------
def _write_tecplot(
    solution: FalknerSkanSolution | FalknerSkanCookeSolution,
    fname: Path,
    problem: SimilarityInputs | None = None,
    title: str | None = None,
    zone_name: str | None = None,
    shooting_result: ShootingResult | None = None,
) -> None:
    """Write similarity solution to Tecplot ASCII format

    Parameters
    ----------
    solution : FalknerSkanSolution | FalknerSkanCookeSolution
        Similarity solution to write.
    fname : Path
        Output file path.
    problem : SimilarityInputs, optional
        Problem specification for metadata.
    title : str, optional
        Title for the Tecplot file header.
    zone_name : str, optional
        Zone name in Tecplot file.
    shooting_result : ShootingResult, optional
        Shooting method convergence info.
    """

    # --------------------------------------------------
    # compute scaled wall-normal velocity and density ratio for output
    #
    # rho/rho_e = 1/tau(eta)
    # - ideal gas at constant pressure: p = rho*R*T = const
    # - so rho/rho_e = T_e/T = 1/tau(eta)
    #
    # (v/u_e)*sqrt(Re_x) = 0.5*(eta*f' - f)
    # - from integrating the continuity equation in similarity space
    # - v = 0.5 * sqrt(nu_e * u_e / x) * (eta*f' - f)
    # - dividing by u_e and multiplying by sqrt(Re_x) = sqrt(u_e*x/nu_e) cancels the physical scales
    # - result is purely a function of the similarity profiles f and f'
    # --------------------------------------------------
    v_ue_sqrtRex = 0.5 * (solution.eta * solution.fp - solution.f)

    # Handle both old (g) and new (tau) naming for temperature
    if hasattr(solution, 'tau'):
        rho_rhoe = 1.0 / solution.tau
    else:
        rho_rhoe = 1.0 / solution.g

    # --------------------------------------------------
    # header information for Tecplot file: title, variable names, auxiliary data
    # --------------------------------------------------

    # title
    if title is None:

        if problem is not None:
            # set title with mach and beta if problem inputs are available
            title = f"Similarity Solution M={problem.mach_edge:.2f} beta={problem.beta:.3f}"
        else:
            # generic title if no problem inputs
            title = "Similarity Solution"

    # zone name
    if zone_name is None:

        if problem is not None:
            # set zone name with mach and beta if problem inputs are available
            zone_name = f"mach_{problem.mach_edge:.1f}_beta{problem.beta:.2f}_{problem.wall_bc}"
        else:
            # generic zone name if no problem inputs
            zone_name = "solution"

    # variable list: w_we appended for FSC solutions
    # Support both old (w) and new (g_cf) naming for crossflow
    variables = ["eta", "u_ue", "T_Te", "rho_rhoe", "v_ue_sqrtRex"]
    if hasattr(solution, "g_cf") or hasattr(solution, "w"):
        variables += ["w_we"]

    # --------------------------------------------------
    # write file
    # --------------------------------------------------
    with open(fname, "w") as f:

        # tecplot header block
        f.write(f'TITLE = "{title}"\n')
        var_str = ", ".join(f'"{v}"' for v in variables)
        f.write(f"VARIABLES = {var_str}\n")
        f.write(f'ZONE T="{zone_name}", I={len(solution.eta)}, F=POINT\n')

        # auxiliary data: problem inputs (if provided)
        if problem is not None:
            f.write(f'AUXDATA mach_edge = "{problem.mach_edge}"\n')
            f.write(f'AUXDATA temp_edge = "{problem.temp_edge}"\n')
            f.write(f'AUXDATA prandtl = "{problem.prandtl}"\n')
            f.write(f'AUXDATA gamma = "{problem.gamma}"\n')
            f.write(f'AUXDATA beta = "{problem.beta}"\n')
            f.write(f'AUXDATA wall_bc_type = "{problem.wall_bc}"\n')
            if problem.sweep_angle != 0.0:
                f.write(f'AUXDATA sweep_angle = "{problem.sweep_angle}"\n')

        # auxiliary data: wall values from solution
        f.write(f'AUXDATA fpp_wall = "{solution.fpp[0]}"\n')

        # Handle both old and new naming
        if hasattr(solution, 'taup'):
            f.write(f'AUXDATA taup_wall = "{solution.taup[0]}"\n')
            f.write(f'AUXDATA tau_wall = "{solution.tau[0]}"\n')
        elif hasattr(solution, 'gp'):
            f.write(f'AUXDATA gp_wall = "{solution.gp[0]}"\n')
            f.write(f'AUXDATA g_wall = "{solution.g[0]}"\n')

        if hasattr(solution, "gcf_p"):
            f.write(f'AUXDATA gcfp_wall = "{solution.gcf_p[0]}"\n')
        elif hasattr(solution, "wp"):
            f.write(f'AUXDATA wp_wall = "{solution.wp[0]}"\n')

        # auxiliary data: convergence info (if provided)
        if shooting_result is not None:
            f.write(f'AUXDATA converged = "{shooting_result.converged}"\n')
            f.write(f'AUXDATA iterations = "{shooting_result.iterations}"\n')

        f.write(f'AUXDATA n_points = "{len(solution.eta)}"\n')
        f.write(f'AUXDATA generated = "{datetime.now().isoformat()}"\n')

        # --------------------------------------------------
        # write data
        #
        # u/u_e = f'(eta)
        # - by definition of the stream function: u = u_e * f'(eta)
        #
        # T/T_e = tau(eta) or g(eta)
        # - by definition of the energy variable: tau = T/T_e (new) or g = T/T_e (old)
        #
        # rho/rho_e = 1/tau(eta) or 1/g(eta)
        # - precomputed above, stored in rho_rhoe
        #
        # (v/u_e)*sqrt(Re_x) = 0.5*(eta*f' - f)
        # - precomputed above from continuity equation, stored in v_ue_sqrtRex
        #
        # w/w_e = g_cf(eta) or w(eta)
        # - Falkner-Skan-Cooke only
        # --------------------------------------------------

        # Get temperature array (support both old and new naming)
        temp_ratio = solution.tau if hasattr(solution, 'tau') else solution.g

        # Get crossflow array (support both old and new naming)
        crossflow = None
        if hasattr(solution, 'g_cf'):
            crossflow = solution.g_cf
        elif hasattr(solution, 'w'):
            crossflow = solution.w

        for i in range(len(solution.eta)):
            row = (
                f"{solution.eta[i]:18.10E} "
                f"{solution.fp[i]:18.10E} "
                f"{temp_ratio[i]:18.10E} "
                f"{rho_rhoe[i]:18.10E} "
                f"{v_ue_sqrtRex[i]:18.10E} "
            )
            if crossflow is not None:
                row += f"{crossflow[i]:18.10E} "
            row = row.rstrip() + "\n"
            f.write(row)


