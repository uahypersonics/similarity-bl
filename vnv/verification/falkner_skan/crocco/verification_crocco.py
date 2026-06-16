# BSD-3-Clause License - see LICENSE file
"""Crocco relation verification for the Falkner-Skan solver.

For Pr = 1, adiabatic wall, the total enthalpy is constant across the boundary
layer (Crocco's theorem). This gives an exact analytic relation between the
temperature and velocity profiles:

    tau(eta) = T/T_e = 1 + (gamma - 1)/2 * Me^2 * (1 - fp^2)

where fp = u/u_e. This holds for any pressure gradient (any beta) when Pr = 1
and the wall is adiabatic.

"""

# --------------------------------------------------
# import necessary modules
# --------------------------------------------------
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from simbl import SimilarityInputs, SolverOptions, solve_similarity

# --------------------------------------------------
# use LaTeX font rendering
# --------------------------------------------------
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,
})

# --------------------------------------------------
# set output directory (figures alongside this script)
# --------------------------------------------------
out_dir = Path(__file__).resolve().parent

# --------------------------------------------------
# cases: (label, mach_edge, beta)
# --------------------------------------------------
GAMMA = 1.4

cases = [
    ("A", 1.5, 0.0),
    ("B", 3.0, 0.0),
]

# --------------------------------------------------
# run all cases
# --------------------------------------------------
print(' ')
print("=" * 50)
print("Crocco relation verification (FS solver, Pr = 1)")
print("=" * 50)

for label, me, beta in cases:

    print(f"\n  Case {label}: Me = {me}, beta = {beta:.4f}")

    # set solver inputs
    inputs = SimilarityInputs(
        mach_edge=me,
        temp_edge=300.0,
        wall_bc="adiabatic",
        beta=beta,
        prandtl=1.0,
        gamma=GAMMA,
    )

    # set solver options
    options = SolverOptions(equations="falkner_skan")

    # run solver
    sol, info = solve_similarity(inputs, options)

    # compute Crocco exact profile
    tau_exact = 1.0 + (GAMMA - 1) / 2.0 * me**2 * (1.0 - sol.fp**2)

    # compute pointwise difference
    diff = sol.tau - tau_exact
    diff_scaled = diff * 1.0e8
    max_diff = float(np.max(np.abs(diff)))

    # print results
    print(f"    Converged   : {info.converged}")
    print(f"    max|diff| : {max_diff:.2e}")

    # --------------------------------------------------
    # plot profiles and difference
    # --------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.5))

    # panel 1: tau profiles
    ax1.plot(sol.tau, sol.eta, "k-", lw=1.5, label="Similarity")
    ax1.plot(tau_exact, sol.eta, "r--", lw=1.2, label="Crocco")
    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"$\eta$")
    ax1.set_xlim(0.5, 3.0)
    ax1.set_xticks(np.arange(0.5, 3.1, 0.5))
    ax1.set_ylim(0.0, 8.0)
    ax1.set_yticks(np.arange(0.0, 8.1, 2.0))
    ax1.legend(frameon=False)
    ax1.grid(True, ls=":", lw=0.5)

    # panel 2: residual
    ax2.plot(diff_scaled, sol.eta, "k-", lw=1.2)
    ax2.axvline(0, color="gray", lw=0.8, ls="--")
    ax2.set_xlabel(r"$(\tau_\mathrm{Similarity} - \tau_\mathrm{Crocco}) \times 10^8$")
    ax2.set_ylabel(r"$\eta$")
    ax2.set_xlim(-4.0, 4.0)
    ax2.set_xticks(np.arange(-4.0, 4.1, 1.0))
    ax2.set_ylim(0.0, 8.0)
    ax2.set_yticks(np.arange(0.0, 8.1, 2.0))
    ax2.grid(True, ls=":", lw=0.5)

    fig.tight_layout()

    me_str = f"{me:.2f}".replace(".", "pt")

    fname = out_dir / f"crocco_me_{me_str}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved {fname.name}")
