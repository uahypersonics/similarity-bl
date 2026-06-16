# BSD-3-Clause License - see LICENSE file
"""Falkner-Skan similarity solution vs CFD baseflow comparison.
"""

# --------------------------------------------------
# import necessary modules
# --------------------------------------------------
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from cfd_io.readers.hdf5 import read_hdf5
from flow_state.transport import get_transport_model

from simbl import SimilarityInputs, SolverOptions, eta2y, solve_similarity

# --------------------------------------------------
# use LaTeX font rendering
# --------------------------------------------------
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,
})

# --------------------------------------------------
# set directories/paths
# --------------------------------------------------

# resolve script directory
script_dir = Path(__file__).resolve().parent

# case directory (parent of script dir)
case_dir = script_dir.parent

# set data directory
data_dir = case_dir / "data"

# set figure directory
figures_dir = case_dir / "figures"

# generate figure directory if it does not already exist
figures_dir.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# read hdf5 file
# --------------------------------------------------
ds = read_hdf5(data_dir / "base_flow_cfdpp.hdf5")

# get attributes from dataset ds
MACH_INF = float(ds.attrs["mach"])
RE1_INF = float(ds.attrs["re1"])
TEMP_INF = float(ds.attrs["temp_fs"])
TEMP_WALL = float(ds.attrs["temp_wall"])
GAMMA = float(ds.attrs["gamma"])
PR = float(ds.attrs["pr"])
R_GAS = float(ds.attrs["rgas"])
visc_model = get_transport_model(str(ds.attrs["visc_law"]))

# derived freestream quantities (ideal gas)
a_inf = np.sqrt(GAMMA * R_GAS * TEMP_INF)
u_inf = MACH_INF * a_inf
mu_inf = visc_model.mu(TEMP_INF)
rho_inf = RE1_INF * mu_inf / u_inf

# get grid and flow field data from dataset ds (squeeze to remove singleton dimensions)
x    = np.squeeze(ds.grid.x)
y    = np.squeeze(ds.grid.y)
uvel = np.squeeze(ds.flow["uvel"].data)
temp = np.squeeze(ds.flow["temp"].data)
dens = np.squeeze(ds.flow["dens"].data)

# --------------------------------------------------
# falkner skan solve
# --------------------------------------------------
sol, info = solve_similarity(
    SimilarityInputs(
        mach_edge=MACH_INF,
        temp_edge=TEMP_INF,
        wall_bc="isothermal",
        temp_wall=TEMP_WALL,
        prandtl=PR,
        gamma=GAMMA,
    ),
    SolverOptions(equations="falkner_skan"),
)

# output for user
print(f"FS converged: {info.converged}, iterations: {info.iterations}")

# compute rho/rho_e from similarity temperature profile
rho_fs = 1.0 / sol.tau

# --------------------------------------------------
# plot 1: FS profiles in eta-space
# --------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3), sharey=True)

# u/u_e vs eta
axes[0].plot(sol.fp,  sol.eta, "k-", lw=1.5)
axes[0].set_xlabel(r"$u/u_e$")
axes[0].set_xlim(0.0, 1.2)
axes[0].set_xticks(np.arange(0.0, 1.21, 0.4))

# T/T_e vs eta
axes[1].plot(sol.tau,   sol.eta, "k-", lw=1.5)
axes[1].set_xlabel(r"$T/T_e$")
axes[1].set_xlim(0.5, 4.5)
axes[1].set_xticks(np.arange(0, 5.01, 1.0))

# rho/rho_e vs eta
axes[2].plot(rho_fs,  sol.eta, "k-", lw=1.5)
axes[2].set_xlabel(r"$\rho/\rho_e$")
axes[2].set_xlim(0.0, 1.2)
axes[2].set_xticks(np.arange(0.0, 1.21, 0.4))

# shared y axis
axes[0].set_ylabel(r"$\eta$")
axes[0].set_ylim(0.0, 6.0)
axes[0].set_yticks(np.arange(0, 7, 1))

for ax in axes:
    ax.grid(True, ls=":", lw=0.5)
fig.tight_layout()

# save figure
fig.savefig(figures_dir / "fs_profiles_eta.png", dpi=150)
print("Saved fs_profiles_eta.png")

# --------------------------------------------------
# plot 2: Similarity vs CFD++
# --------------------------------------------------
x_loc = [0.1, 0.2, 0.3, 0.4, 0.5]

# pre-compute profiles at each x station
stations = []

for x_target in x_loc:
    i = int(np.argmin(np.abs(x[:, 0] - x_target)))
    x_i = float(x[i, 0])

    # edge = top-of-domain (outer boundary, non-dim)
    uvel_e = float(uvel[i, -1])
    temp_e = float(temp[i, -1])
    dens_e = float(dens[i, -1])

    # dimensional edge values for Levy-Lees eta scale
    temp_e_dim = temp_e * TEMP_INF
    uvel_e_dim = uvel_e * u_inf
    dens_e_dim = dens_e * rho_inf
    visc_e = visc_model.mu(temp_e_dim)

    # CFD++ profiles normalized by edge quantities
    y_i = y[i, :]
    uvel_cfd = uvel[i, :] / uvel_e
    temp_cfd = temp[i, :] / temp_e
    dens_cfd = dens[i, :] / dens_e

    # map density-weighted similarity eta back to physical y
    y_fs = eta2y(
        eta=sol.eta,
        tau=sol.tau,
        x=x_i,
        dens_edge=dens_e_dim,
        uvel_edge=uvel_e_dim,
        visc_edge=visc_e,
        equations="falkner_skan",
    )

    # convert to mm for cleaner tick labels
    y_i  = y_i  * 1e3
    y_fs = y_fs * 1e3

    # fixed y extent (mm)
    y_max = 4.0

    stations.append(dict(x_i=x_i, y_i=y_i, y_fs=y_fs, y_max=y_max,
                         uvel_cfd=uvel_cfd, temp_cfd=temp_cfd, dens_cfd=dens_cfd))

# one figure per quantity
quantities = [
    ("uvel_cfd", "fp",  r"$u/u_e$",       "similarity_vs_cfdpp_uvel.png"),
    ("temp_cfd", "tau", r"$T/T_e$",        "similarity_vs_cfdpp_temp.png"),
    ("dens_cfd",  None, r"$\rho/\rho_e$",  "similarity_vs_cfdpp_dens.png"),
]

for cfd_key, fs_attr, xlabel, fname in quantities:
    fig, axes = plt.subplots(1, len(stations), figsize=(6.5, 2.5), sharey=True)

    lines = None
    for col, st in enumerate(stations):
        ax = axes[col]
        q_cfd = st[cfd_key]
        # fs_attr is None for rho (derived from sol.g)
        q_fs = rho_fs if fs_attr is None else getattr(sol, fs_attr)

        l1, = ax.plot(q_cfd, st["y_i"],  "k-",  lw=1.0, label="CFD++")
        l2, = ax.plot(q_fs,  st["y_fs"], "r--", lw=1.2, label="Similarity")
        ax.set_ylim(0.0, st["y_max"])
        ax.set_yticks(np.arange(0, 5, 1))
        ax.set_title(f"x = {st['x_i']:.2f} m")
        ax.grid(True, ls=":", lw=0.5)
        if col == 0:
            ax.set_ylabel("y (mm)")
            lines = [l1, l2]
        if col == len(stations) // 2:
            ax.set_xlabel(xlabel)

    # single legend centered above all panels
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.legend(lines, ["CFD++", "Similarity"], loc="upper center",
               ncol=2, bbox_to_anchor=(0.5, 0.99), frameon=False)
    fig.savefig(figures_dir / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")
