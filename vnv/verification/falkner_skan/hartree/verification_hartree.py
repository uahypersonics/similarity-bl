"""Hartree (1937) verification for the Falkner-Skan solver.

Runs the FS solver at Mach=0.01 (effectively incompressible), adiabatic wall,
for the beta values tabulated by Hartree (1937) Table I and compares f''(0)
against the published values.

Reference
---------
Hartree, D. R. (1937). On an equation occurring in Falkner and Skan's
approximate treatment of the equations of the boundary layer.
Proc. Cambridge Phil. Soc., 33(2), 223-239.

"""


# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from simbl import SimilarityInputs, default_options, solve_similarity

# --------------------------------------------------
# reference data from Hartree (1937) Table I
# --------------------------------------------------
data_ref = {
    0.3:  0.7745,
    0.4:  0.854,
    0.6:  0.995,
    0.8:  1.121,
    1.2:  1.335,
    1.6:  1.522,
}

# --------------------------------------------------
# set conditions
# --------------------------------------------------

# mach number: 0.01 is effectively incompressible
mach   = 0.01
# edge temperature [K]
temp_edge = 300.0

# --------------------------------------------------
# main: run verification
# --------------------------------------------------

# output for user
print(f"\nHartree (1937) verification (FS solver: adiabatic wall, Mach={mach}, T_e={temp_edge} K)\n")

print(f"{'beta':>6}  {'ref f´´(0)':>10}  {'solver':>10}  {'f´´(0) guess':>12}  {'g guess':>9}  {'rel err %':>10}  status")
print("-" * 75)

# set tolerances for pass/warn/fail

# pass tolerance: < 0.5% relative error
PASS_TOL = 0.5
# warn tolerance: < 2.0% relative error
WARN_TOL = 2.0

# initialize counters for pass/warn/fail/no convergence
n_pass = n_warn = n_fail = n_nc = 0


profiles: list[tuple[float, object, object, object, float, float]] = []   # (beta, eta, fp, fpp, initial_fpp_seed, initial_gvar_seed)

# seed for the first beta
initial_fpp  = 0.7
initial_gvar = 1.0

# walk the Hartree betas in ascending order.
# between consecutive betas, sub-step in 0.01
# The seed carries forward
prev_beta = 0.0

# loop over Hartree beta values
for beta, ref in sorted(data_ref.items()):

    # continuation: walk from prev_beta to beta in 0.01 steps
    n_betas = max(1, round((beta - prev_beta) / 0.02))
    beta_arr = np.linspace(prev_beta, beta, n_betas + 1)[1:]

    # capture the seed used for this target beta (before sub-stepping updates it)
    initial_fpp_seed  = initial_fpp
    initial_gvar_seed = initial_gvar

    # initialize is_converged flag for current beta value
    # if any sub-step fails to converge, the flag is set to False and the beta is reported as NO_CONV
    is_converged = True

    # loop over
    for beta_loc in beta_arr:

        # set up inputs and options for the FS solver
        input = SimilarityInputs(mach_edge=mach, temp_edge=temp_edge, wall_bc="adiabatic", beta=float(beta_loc))

        # set solver options
        options = default_options(input)

        # run solver
        solution, result = solve_similarity(input, options, initial_fpp=initial_fpp, initial_gvar=initial_gvar)

        # check convergence, if not converged, break out of the sub-step loop
        if not result.converged:
            is_converged = False
            break

        # update seed with converged shooting vars for the next sub-step
        initial_fpp  = float(result.shooting_vars[0])
        initial_gvar = float(result.shooting_vars[1])

    # report results after substep loop
    if not is_converged:
        print(f"{beta:>6.2f}  {ref:>10.5f}  {'---':>10}  {'---':>12}  {'---':>9}  {'---':>10}  NO CONV")
        n_nc += 1
        # do NOT update prev_beta — retry from the last good position
        continue

    # update prev_beta to the current beta for the next iteration
    prev_beta = beta

    fpp_wall = float(result.solution[2, 0])
    profiles.append((beta, result.eta, result.solution[1, :], result.solution[2, :], initial_fpp_seed, initial_gvar_seed))

    rel_err = abs(fpp_wall - ref) / ref * 100.0
    if rel_err < PASS_TOL:
        status = "PASS"
        n_pass += 1
    elif rel_err < WARN_TOL:
        status = "WARN"
        n_warn += 1
    else:
        status = "FAIL"
        n_fail += 1

    print(f"{beta:>6.2f}  {ref:>10.5f}  {fpp_wall:>10.5f}  {initial_fpp_seed:>12.4f}  {initial_gvar_seed:>9.4f}  {rel_err:>10.4f}  {status}")

print("-" * 75)
print(f"Result: {n_pass} PASS  {n_warn} WARN  {n_fail} FAIL  {n_nc} NO_CONV"
      f"  (of {len(data_ref)} cases)")

# exit with nonzero if any failures or convergence issues (avoid plotting in that case)
if n_fail or n_nc:
    sys.exit(1)

# --------------------------------------------------
# plots: one figure per beta, two panels (f', f'')
# --------------------------------------------------

mpl.rcParams.update({
    "font.family":        "serif",
    "mathtext.fontset":   "cm",
    "axes.grid":          True,
    "grid.color":         "#b0b0b0",
    "grid.linestyle":     "-",
    "grid.linewidth":     0.7,
    "grid.alpha":         1.0,
    "axes.spines.top":    True,
    "axes.spines.right":  True,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "figure.dpi":         150,
})

# plots saved directly next to this script
plots_dir = Path(__file__).parent

saved: list[tuple[float, float, Path]] = []

# unpack profiles and plot each beta
for beta, eta, fp, fpp, ig_fpp, ig_gvar in profiles:

    ref = data_ref[beta]

    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(3.25, 2.5), sharey=True,
        gridspec_kw={"wspace": 0.25},
    )

    # left panel: f'(eta)
    ax_l.plot(fp, eta, color='blue', linewidth=1.5)

    # set axes labels, limits, ticks
    ax_l.set_xlabel(r"$f'$",  fontsize=11)
    ax_l.set_ylabel(r"$\eta$", fontsize=11)
    ax_l.set_xlim(0.0, 1.2)
    ax_l.set_xticks([0.0, 0.6, 1.2])
    ax_l.set_ylim(0, 4)
    ax_l.set_yticks([0, 1, 2, 3, 4])
    ax_l.tick_params(labelsize=11)

    # right panel: f''(eta)
    ax_r.plot(fpp, eta, color="blue", linewidth=1.5)

    # mark f''(0) at the wall with a red filled circle
    ax_r.scatter([float(fpp[0])], [float(eta[0])], color="red", s=25, zorder=5)

    # label in the upper-right corner: red circle marker + value (no arrow)
    ax_r.text(
        0.97, 0.97,
        rf"$\bullet\ f''(0) = {float(fpp[0]):.4f}$",
        transform=ax_r.transAxes,
        fontsize=9, color="red",
        ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", linewidth=0.8, alpha=0.9),
    )

    # set axes labels, limits, ticks
    ax_r.set_xlabel(r"$f''$", fontsize=11)
    ax_r.tick_params(labelleft=False, labelsize=11)
    ax_r.set_xlim(0.0, 2.0)
    ax_r.set_xticks([0.0, 1.0, 2.0])

    fig.suptitle(
        rf"$M_e = {mach}$,  $T_e = {temp_edge:.0f}\ \mathrm{{K}}$,  $\beta = {beta}$",
        fontsize=11, y=1.01,
    )

    # set output path
    beta_str = f"{beta:.1f}".replace(".", "pt")
    fpath = plots_dir / f"profiles_beta_{beta_str}.png"

    # save figure
    fig.savefig(fpath, dpi=150, bbox_inches="tight")

    # close figure to free memory
    plt.close(fig)

    # store for markdown table
    saved.append((beta, float(fpp[0]), fpath, ig_fpp, ig_gvar))

# -- write results.md -------------------------------------------------------
results_path = plots_dir / "results.md"
with results_path.open("w") as f:
    f.write("# Hartree (1937) Verification Results\n\n")
    f.write(f"Solver: Falkner-Skan, adiabatic wall, Mach={mach}, T_e={temp_edge} K\n\n")
    f.write("| β | f''(0) ref | f''(0) solver | f''(0) guess | g guess | Profile |\n")
    f.write("|---|-----------|---------------|-------------|---------|---------|\n")
    for beta, fpp_wall, path, ig_fpp, ig_gvar in saved:
        ref = data_ref[beta]
        f.write(
            f"| {beta} | {ref} | {fpp_wall:.5f} "
            f"| {ig_fpp:.4f} | {ig_gvar:.4f} "
            f"| [view]({path.name}) |\n"
        )
print(f"\nResults written to {results_path}")
