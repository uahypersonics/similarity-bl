"""beta=0.6 initial-guess comparison: cold start vs continuation seed."""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
import warnings
from dataclasses import replace
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from simbl import SimilarityInputs, default_options, solve_similarity

# --------------------------------------------------
# conditions
# --------------------------------------------------
mach = 0.01
temp_edge = 300.0
beta_target = 0.6
beta_seed = 0.4
eta_max_vals = [6, 8, 10, 12]

input = SimilarityInputs(mach_edge=mach, temp_edge=temp_edge, wall_bc="adiabatic", beta=beta_target)

# --------------------------------------------------
# compute the initial seed from a smaller beta value
# --------------------------------------------------

# step 1: beta = 0.4 cold start (f''(0) = 0.7, g = 1.0)
input_seed  = SimilarityInputs(mach_edge=mach, temp_edge=temp_edge, wall_bc="adiabatic", beta=beta_seed)
options_seed = default_options(input_seed)
solution_seed, result_seed = solve_similarity(input_seed, options_seed, initial_fpp=0.7, initial_gvar=1.0)

# step 2: walk to beta = 0.6 in steps of 0.02

# initial seed from converged beta = 0.4 solution
fpp_seed  = float(result_seed.shooting_vars[0])
gvar_seed = float(result_seed.shooting_vars[1])

for beta_loc in np.linspace(beta_seed, beta_target, round((beta_target - beta_seed) / 0.02) + 1)[1:]:

    input_loc  = SimilarityInputs(mach_edge=mach, temp_edge=temp_edge, wall_bc="adiabatic", beta=float(beta_loc))

    options_loc = default_options(input_loc)

    solution_loc, result_loc = solve_similarity(input_loc, options_loc, initial_fpp=fpp_seed, initial_gvar=gvar_seed)

    # update the seed for the next iteration if converged otherwise raise error and exit
    if result_loc.converged:
        fpp_seed  = float(result_loc.shooting_vars[0])
        gvar_seed = float(result_loc.shooting_vars[1])
    else:
        raise RuntimeError(f"beta continuation failed at beta={beta_loc:.2f}")

print(f"continuation seed: f''(0)={fpp_seed:.4f},  g={gvar_seed:.4f}")

# --------------------------------------------------
# solve cold start and continuation for each eta_max
# --------------------------------------------------

# {eta_max: (res_cold, res_cont)}
results = {}

# loop over eta_max values, solve with cold start and continuation seed
for eta_max in eta_max_vals:

    # overwrite the default options with the current eta_max value (eta max will otherwise be determined based on beta)
    options = replace(default_options(input), eta_max=float(eta_max))

    # compute with cold start: suppress the spurious-solution warning (it is expected to converge to a spurious branch for large eta_max values)
    with warnings.catch_warnings(record=True):

        warnings.simplefilter("always")

        solution_cold, result_cold = solve_similarity(input, options, initial_fpp=0.7, initial_gvar=1.0)

    # compute with continuation seed
    solution_cont, result_cont = solve_similarity(input, options, initial_fpp=fpp_seed, initial_gvar=gvar_seed)

    # store results for plotting
    results[eta_max] = (result_cold, result_cont)

    # print summary to console
    print(
        f"eta_max={eta_max:2d}  "
        f"cold: converged={result_cold.converged}, f''(eta_max)={result_cold.solution[2,-1]:+.3f}  |  "
        f"cont: converged={result_cont.converged}, f''(eta_max)={result_cont.solution[2,-1]:.2e}"
    )

# --------------------------------------------------
# plot results
# --------------------------------------------------

# set global matplotlib rcParams for consistent figure style (latex type font, grid line thickness, etc.)
mpl.rcParams.update({
    "font.family":      "serif",
    "mathtext.fontset": "cm",
    "axes.grid":        True,
    "grid.color":       "#b0b0b0",
    "grid.linewidth":   0.7,
    "grid.alpha":       1.0,
    "figure.dpi":       150,
})

# generate 1 plot for each eta_max value, showing f' and f'' vs eta for cold start and continuation seed
for eta_max in eta_max_vals:

    # unpack results for this eta_max
    result_cold, result_cont = results[eta_max]

    # set up figure and axes (2 panels, shared y-axis)
    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(3.25, 2.5), sharey=True,
        gridspec_kw={"wspace": 0.25},
    )

    # left panel: f'(eta)

    # continuation plotted first as solid blue (reference)
    line_cont, = ax_l.plot(result_cont.solution[1, :], result_cont.eta, color="blue", linewidth=1.5,label=r"$\beta$-continuation seed")

    # cold start plotted on top as red dashed line
    line_cold, = ax_l.plot(result_cold.solution[1, :], result_cold.eta, color="red", linewidth=1.5,linestyle="--", label=r"cold start ($f''(0) = 0.7$)")

    # add horizontal line at eta_max
    ax_l.axhline(eta_max, color="k", linewidth=0.8, linestyle="--", alpha=0.4)
    # add vertical line at f'(eta)=0
    ax_l.axvline(0, color="k", linewidth=0.6, linestyle=":")

    # set axes labels, limits, ticks, and font sizes
    ax_l.set_xlabel(r"$f'$", fontsize=11)
    ax_l.set_ylabel(r"$\eta$", fontsize=11)
    ax_l.set_xlim(-1.5, 1.5)
    ax_l.set_ylim(0, 12)
    ax_l.set_xticks([-1.5, 0.0, 1.5])
    ax_l.tick_params(labelsize=11)

    # right panel: f''(eta)

    ax_r.plot(result_cont.solution[2, :], result_cont.eta, color="blue", linewidth=1.5)
    ax_r.plot(result_cold.solution[2, :], result_cold.eta, color="red", linewidth=1.5,linestyle="--")

    # add horizontal line at eta_max
    ax_r.axhline(eta_max, color="k", linewidth=0.8, linestyle="--", alpha=0.4)
    # add vertical line at f''(eta)=0
    ax_r.axvline(0, color="k", linewidth=0.6, linestyle=":")

    # set axes labels, limits, ticks, and font sizes
    ax_r.set_xlabel(r"$f''$", fontsize=11)
    ax_r.tick_params(labelleft=False, labelsize=11)
    ax_r.set_xlim(-2.0, 2.0)
    ax_r.set_xticks([-2.0, 0.0, 2.0])

    # legend above the figure, outside the axes
    fig.legend(
        handles=[line_cont, line_cold],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.06),
        ncol=2,
        fontsize=8,
        framealpha=0.9,
    )

    # add a title above the figure
    fig.suptitle(
        rf"$\beta = {beta_target}$,  $M_e = {mach}$,  $T_e = {temp_edge:.0f}\ \mathrm{{K}}$,"
        rf"  $\eta_{{max}} = {eta_max}$",
        fontsize=11, y=1.17,
    )

    # save figure to file
    fpath = Path(__file__).parent / f"initial_guess_sensitivity_beta_0pt6_eta_max_{eta_max:02d}.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {fpath.name}")
