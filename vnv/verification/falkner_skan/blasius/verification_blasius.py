"""Blasius verification for the Falkner-Skan solver.

At low Mach number (Me -> 0), adiabatic wall, beta = 0 (flat plate), the
compressible FS system reduces to the classical Blasius equation:

    f''' + 0.5 * f * f'' = 0,  f''(0) = 0.469600 (Blasius constant)

Note: f''(0) in the compressible Levy-Lees system is NOT the same as the
incompressible Blasius f''(0) = 0.332057. In the Levy-Lees transform the
skin friction is scaled differently. At Me = 0.01 the value converges to
the compressible-limit wall shear, which equals 0.4696 in the Levy-Lees
non-dimensionalization. See documentation for details ().

"""

# --------------------------------------------------
# import necessary modules
# --------------------------------------------------
from pathlib import Path

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
# set reference value
# --------------------------------------------------

# Blasius wall shear in Levy-Lees compressible scaling (Me -> 0 limit)
FPP0_REFERENCE = 0.469600

# --------------------------------------------------
# set Mach numbers to sweep
# --------------------------------------------------

MACH_EDGES = [0.1, 0.01, 0.001, 0.0001]

# --------------------------------------------------
# helper: solve and plot profiles for one Me
# --------------------------------------------------

def run_case(mach_edge, out_dir):
    """Solve the FS system at the given Me, plot profiles, return metrics.

    Args:
        mach_edge: edge Mach number
        out_dir: Path to save the figure

    Returns:
        tuple: (fpp0, tau_max_dev, converged)
    """

    # set solver inputs
    inputs = SimilarityInputs(
        mach_edge=mach_edge,
        temp_edge=300.0,
        wall_bc="adiabatic",
        beta=0.0,
        prandtl=0.72,
        gamma=1.4,
    )

    # set solver options
    options = SolverOptions(equations="falkner_skan")

    # run solver
    sol, info = solve_similarity(inputs, options)

    # extract wall shear f''(0)
    fpp0 = float(sol.fpp[0])

    # compute max deviation of tau from 1.0 (should be ~1 everywhere at Me -> 0)
    tau_max_dev = float(np.max(np.abs(sol.tau - 1.0)))

    # print results on screen
    print(f"\n  Me = {mach_edge}")
    print(f"  Converged      : {info.converged}")
    print(f"  Iterations     : {info.iterations}")
    print(f"  f''(0)         : {fpp0:.6f}  (reference: {FPP0_REFERENCE:.6f})")
    print(f"  |f''(0) - ref| : {abs(fpp0 - FPP0_REFERENCE):.2e}")
    print(f"  max|tau - 1|   : {tau_max_dev:.2e}")

    # plot profiles
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3), sharey=True)

    # set eta limit
    eta_max = 8

    # f' (velocity)
    axes[0].plot(sol.fp, sol.eta, color="blue", linewidth=1.5)
    axes[0].set_xlabel(r"$f'$")
    axes[0].set_ylabel(r"$\eta$")
    axes[0].set_xlim(0, 1.2)
    axes[0].set_xticks([0.0, 0.4, 0.8, 1.2])
    axes[0].set_ylim(0, eta_max)
    axes[0].set_yticks(range(0, eta_max + 1, 2))
    axes[0].grid(True)

    # f'' (shear)
    axes[1].plot(sol.fpp, sol.eta, color="blue", linewidth=1.5)
    axes[1].axvline(FPP0_REFERENCE, color="red", linestyle="--", linewidth=1.2,
                    label=r"$f''(0)_\mathrm{ref}$")
    axes[1].set_xlabel(r"$f''$")
    axes[1].set_xlim(0, 0.5)
    axes[1].set_xticks([0.00, 0.25, 0.50])
    axes[1].set_ylim(0, eta_max)
    axes[1].set_yticks(range(0, eta_max + 1, 2))
    axes[1].legend(fontsize=10, loc="upper center")
    axes[1].grid(True)

    # (tau - 1) scaled to O(1) based on actual max deviation
    tau_exp = int(np.floor(-np.log10(tau_max_dev)))
    tau_scale = 10 ** tau_exp
    axes[2].plot((sol.tau - 1.0) * tau_scale, sol.eta, color="blue", linewidth=1.5)
    axes[2].axvline(0.0, color="red", linestyle="--", linewidth=1.2,
                    label=r"$\tau - 1 = 0$ (exact)")
    axes[2].set_xlabel(rf"$(\tau - 1) \times 10^{{{tau_exp}}}$")
    axes[2].set_xlim(-1, 1)
    axes[2].set_xticks([-1, 0, 1])
    axes[2].set_ylim(0, eta_max)
    axes[2].set_yticks(range(0, eta_max + 1, 2))
    axes[2].legend(fontsize=10, loc="upper center")
    axes[2].grid(True)

    fig.tight_layout()

    # write figure next to the script
    me_str = f"{mach_edge:.4f}".replace(".", "pt")
    fname = out_dir / f"verification_blasius_me_{me_str}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname.name}")

    return fpp0, tau_max_dev, info.converged


# --------------------------------------------------
# run all cases
# --------------------------------------------------

print(" ")
print("=" * 50)
print("Blasius verification (FS solver, Me -> 0 limit)")
print("=" * 50)

# output directory for figures is the same directory as this script
out_dir = Path(__file__).resolve().parent

# collect results for convergence plot
fpp0_devs = []
tau_devs = []

for me in MACH_EDGES:
    fpp0, tau_max_dev, converged = run_case(me, out_dir)
    fpp0_devs.append(abs(fpp0 - FPP0_REFERENCE))
    tau_devs.append(tau_max_dev)

print("\n" + "=" * 50)

# --------------------------------------------------
# convergence plot: errors vs Me
# --------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(6.5, 3))

# f''(0) deviation vs Me
axes[0].loglog(MACH_EDGES, fpp0_devs, "o-", color="blue", linewidth=1.5)
axes[0].set_xlabel(r"$M_e$")
axes[0].set_ylabel(r"$|f''(0) - f''(0)_\mathrm{ref}|$")
axes[0].set_xlim(1e-5, 10)
axes[0].set_ylim(1e-10, 1e-2)
axes[0].grid(True, which="major")

# max|tau - 1| vs Me
axes[1].loglog(MACH_EDGES, tau_devs, "o-", color="blue", linewidth=1.5)
axes[1].set_xlabel(r"$M_e$")
axes[1].set_ylabel(r"$\max|\tau - 1|$")
axes[1].set_xlim(1e-5, 10)
axes[1].set_ylim(1e-10, 1e-2)
axes[1].grid(True, which="major")

fig.tight_layout()

fname = out_dir / "verification_blasius_convergence.png"
fig.savefig(fname, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\n  Saved {fname.name}")

