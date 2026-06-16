# BSD-3-Clause License - see LICENSE file
"""Falkner-Skan-Cooke similarity solution vs swept flat plate CFD reference."""

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
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 12,
    }
)

# --------------------------------------------------
# set directories/paths
# --------------------------------------------------
script_dir = Path(__file__).resolve().parent
case_dir = script_dir.parent
data_dir = case_dir / "data"
figures_dir = case_dir / "figures"

# generate figure directory if it does not already exist
figures_dir.mkdir(parents=True, exist_ok=True)

# set reference data path
reference_path = data_dir / "base_flow_cfdpp.hdf5"


# --------------------------------------------------
# read reference hdf5 data (cfd++ base flow)
# --------------------------------------------------
ds = read_hdf5(reference_path)

# --------------------------------------------------
# get attributes from dataset ds
# --------------------------------------------------
attrs = ds.attrs

# extract specific attributes needed for similarity solve and CFD profile normalization
mach_inf = float(attrs["mach"])
temp_inf = float(attrs["temp_inf"])
dens_inf = float(attrs["dens_inf"])
temp_wall = float(attrs["t_wall"])
gamma = float(attrs["gamma"])
prandtl = float(attrs["pr"])
rgas = float(attrs["rgas"])
beta = 0.0

# set viscosity model
visc_model = get_transport_model(str(attrs["visc_law"]))

# compute freestream velocity used by the nondimensional CFD fields
a_inf = np.sqrt(gamma * rgas * temp_inf)
uvel_inf = mach_inf * a_inf

# --------------------------------------------------
# get grid and flow fields from dataset ds
# --------------------------------------------------
x = ds.grid.x
y = ds.grid.y

uvel = ds.flow["uvel"].data
vvel = ds.flow["vvel"].data
wvel = ds.flow["wvel"].data
temp = ds.flow["temp"].data
dens = ds.flow["dens"].data

# --------------------------------------------------
# determine boundary layer thickness (99% of freestream w-velocity) at each streamwise station
# --------------------------------------------------
n_x = wvel.shape[0]
j99 = np.zeros(n_x, dtype=int)

delta99 = np.zeros(n_x)

# extract boundary layer edge at 99% of the freestream w-velocity
for i in range(n_x):
    w99 = 0.99 * wvel[i, -1]
    j99[i] = np.searchsorted(wvel[i, :], w99)
    j99[i] = min(j99[i], wvel.shape[1] - 1)
    delta99[i] = y[i, j99[i]]

# --------------------------------------------------
# extract edge conditions at each streamwise station
# --------------------------------------------------
uvel_e = uvel[np.arange(n_x), j99]
vvel_e = vvel[np.arange(n_x), j99]
wvel_e = wvel[np.arange(n_x), j99]
temp_e = temp[np.arange(n_x), j99]
dens_e = dens[np.arange(n_x), j99]

# --------------------------------------------------
# compute dimensional edge conditions
# --------------------------------------------------
temp_e_dim = temp_e * temp_inf
dens_e_dim = dens_e * dens_inf
uvel_e_dim = uvel_e * uvel_inf
vvel_e_dim = vvel_e * uvel_inf
wvel_e_dim = wvel_e * uvel_inf

# compute velocity magnitude at edge
vel_mag_e = np.sqrt(uvel_e_dim**2 + vvel_e_dim**2 + wvel_e_dim**2)

# compute speed of sound at edge
a_e = np.sqrt(gamma * rgas * temp_e_dim)

# compute mach number at edge
mach_e = vel_mag_e / a_e

# compute sweep angle at edge in degrees
sweep_angle_edge_deg = np.rad2deg(np.arctan2(wvel_e, uvel_e))

# compute viscosity at edge using the specified viscosity model
visc_e = np.array([visc_model.mu(float(temp_value)) for temp_value in temp_e_dim])

x_wall = x[:, 0]

# --------------------------------------------------
# nondimensionalize CFD profiles using local edge conditions
# --------------------------------------------------
uvel_cfd = uvel / uvel_e[:, None]
wvel_cfd = wvel / wvel_e[:, None]
temp_cfd = temp / temp_e[:, None]
dens_cfd = dens / dens_e[:, None]
y_cfd_mm = y * 1.0e3
y_max_mm = 3.0 * delta99 * 1.0e3

# --------------------------------------------------
# solve FSC at each reference station
# --------------------------------------------------
uvel_sim = None
wvel_sim = None
temp_sim = None
dens_sim = None
eta_sim = None
y_sim_mm = None

for i, x_i in enumerate(x_wall):
    # set up similarity solver inputs
    inputs = SimilarityInputs(
        mach_edge=float(mach_e[i]),
        temp_edge=float(temp_e_dim[i]),
        wall_bc="isothermal",
        temp_wall=temp_wall,
        prandtl=prandtl,
        gamma=gamma,
        beta=beta,
        sweep_angle=float(sweep_angle_edge_deg[i]),
    )

    # set solver options
    options = SolverOptions(equations="falkner_skan_cooke")

    # run solver
    sol, info = solve_similarity(inputs, options)

    # output for user
    print(
        f"x = {x_i:.3f} m, "
        f"Me = {mach_e[i]:.4f}, "
        f"Lambda = {sweep_angle_edge_deg[i]:.2f} deg, "
        f"converged = {info.converged}"
    )

    # map density-weighted similarity eta back to physical y
    y_sim_m = eta2y(
        eta=sol.eta,
        tau=sol.tau,
        x=float(x_wall[i]),
        dens_edge=float(dens_e_dim[i]),
        uvel_edge=float(uvel_e_dim[i]),
        visc_edge=float(visc_e[i]),
        equations="falkner_skan_cooke",
    )

    # compute density from constant-pressure similarity relation
    dens_sim = 1.0 / sol.tau

    # allocate similarity arrays after the first solve sets the eta-grid size
    if i == 0:
        n_eta = sol.eta.size
        uvel_sim = np.zeros((n_x, n_eta))
        wvel_sim = np.zeros((n_x, n_eta))
        temp_sim = np.zeros((n_x, n_eta))
        dens_sim_all = np.zeros((n_x, n_eta))
        eta_sim = np.zeros((n_x, n_eta))
        y_sim_mm = np.zeros((n_x, n_eta))

    # store similarity solution for plotting
    uvel_sim[i, :] = sol.fp
    wvel_sim[i, :] = sol.g
    temp_sim[i, :] = sol.tau
    dens_sim_all[i, :] = dens_sim
    eta_sim[i, :] = sol.eta
    y_sim_mm[i, :] = y_sim_m * 1.0e3

# --------------------------------------------------
# plot FSC profiles in eta-space at each station
# --------------------------------------------------
for i, x_i in enumerate(x_wall):
    fig, axes = plt.subplots(1, 4, figsize=(6.2, 2.15), sharey=True)

    axes[0].plot(uvel_sim[i, :], eta_sim[i, :], "k-", lw=1.5)
    axes[0].set_xlabel(r"$u/u_e$")
    axes[0].set_xlim(0.0, 1.5)
    axes[0].set_xticks(np.arange(0.0, 1.51, 0.5))

    axes[1].plot(wvel_sim[i, :], eta_sim[i, :], "k-", lw=1.5)
    axes[1].set_xlabel(r"$w/w_e$")
    axes[1].set_xlim(0.0, 1.5)
    axes[1].set_xticks(np.arange(0.0, 1.51, 0.5))

    axes[2].plot(temp_sim[i, :], eta_sim[i, :], "k-", lw=1.5)
    axes[2].set_xlabel(r"$T/T_e$")
    axes[2].set_xlim(0.0, 5.0)
    axes[2].set_xticks(np.arange(0.0, 6.1, 2.0))

    axes[3].plot(dens_sim_all[i, :], eta_sim[i, :], "k-", lw=1.5)
    axes[3].set_xlabel(r"$\rho/\rho_e$")
    axes[3].set_xlim(0.0, 1.2)
    axes[3].set_xticks(np.arange(0.0, 1.21, 0.6))

    axes[0].set_ylabel(r"$\eta$")
    axes[0].set_ylim(0.0, 6.0)
    axes[0].set_yticks(np.arange(0.0, 6.1, 2.0))

    fig.suptitle(
        rf"$x = {x_i:.3f}\ \mathrm{{m}}$, "
        rf"$M_e = {mach_e[i]:.3f}$, "
        rf"$T_e = {temp_e_dim[i]:.1f}\ \mathrm{{K}}$, "
        rf"$\Lambda_e = {sweep_angle_edge_deg[i]:.2f}^\circ$",
        fontsize=10,
        y=0.95,
    )

    for ax in axes:
        ax.grid(True, ls=":", lw=0.5)

    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.235, top=0.83, wspace=0.34)

    x_str = f"{x_i:.3f}".replace(".", "pt")
    filename = f"fsc_profiles_eta_x_{x_str}.png"
    fig.savefig(figures_dir / filename, dpi=150, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Saved {filename}")

# --------------------------------------------------
# plot similarity vs CFD profiles in physical space
# --------------------------------------------------

# plot streamwise velocity profiles
fig, axes = plt.subplots(1, n_x, figsize=(6.5, 2.5), sharey=True)

for i, x_i in enumerate(x_wall):
    ax = axes[i]
    ax.plot(uvel_cfd[i, :], y_cfd_mm[i, :], "b-", lw=1.0, label="CFD++")
    ax.plot(uvel_sim[i, :], y_sim_mm[i, :], "r--", lw=1.2, label="Similarity")
    ax.set_ylim(0.0, y_max_mm[i])
    ax.set_title(f"x = {x_i:.2f} m")
    ax.grid(True, ls=":", lw=0.5)

axes[0].set_ylabel("y (mm)")
axes[n_x // 2].set_xlabel(r"$u/u_e$", labelpad=1)
fig.tight_layout(rect=[0, 0, 1, 0.9])
lines, labels = axes[0].get_legend_handles_labels()
fig.legend(lines, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.99), frameon=False)
filename = "similarity_vs_cfd_uvel.png"
fig.savefig(figures_dir / filename, dpi=150, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)

# plot spanwise velocity profiles
fig, axes = plt.subplots(1, n_x, figsize=(6.5, 2.5), sharey=True)

for i, x_i in enumerate(x_wall):
    ax = axes[i]
    ax.plot(wvel_cfd[i, :], y_cfd_mm[i, :], "b-", lw=1.0, label="CFD++")
    ax.plot(wvel_sim[i, :], y_sim_mm[i, :], "r--", lw=1.2, label="Similarity")
    ax.set_ylim(0.0, y_max_mm[i])
    ax.set_title(f"x = {x_i:.2f} m")
    ax.grid(True, ls=":", lw=0.5)

axes[0].set_ylabel("y (mm)")
axes[n_x // 2].set_xlabel(r"$w/w_e$", labelpad=1)
fig.tight_layout(rect=[0, 0, 1, 0.9])
lines, labels = axes[0].get_legend_handles_labels()
fig.legend(lines, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.99), frameon=False)
filename = "similarity_vs_cfd_wvel.png"
fig.savefig(figures_dir / filename, dpi=150, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)

# plot temperature profiles
fig, axes = plt.subplots(1, n_x, figsize=(6.5, 2.5), sharey=True)

for i, x_i in enumerate(x_wall):
    ax = axes[i]
    ax.plot(temp_cfd[i, :], y_cfd_mm[i, :], "b-", lw=1.0, label="CFD++")
    ax.plot(temp_sim[i, :], y_sim_mm[i, :], "r--", lw=1.2, label="Similarity")
    ax.set_ylim(0.0, y_max_mm[i])
    ax.set_title(f"x = {x_i:.2f} m")
    ax.grid(True, ls=":", lw=0.5)

axes[0].set_ylabel("y (mm)")
axes[n_x // 2].set_xlabel(r"$T/T_e$", labelpad=1)
fig.tight_layout(rect=[0, 0, 1, 0.9])
lines, labels = axes[0].get_legend_handles_labels()
fig.legend(lines, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.99), frameon=False)
filename = "similarity_vs_cfd_temp.png"
fig.savefig(figures_dir / filename, dpi=150, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)

# plot density profiles
fig, axes = plt.subplots(1, n_x, figsize=(6.5, 2.5), sharey=True)

for i, x_i in enumerate(x_wall):
    ax = axes[i]
    ax.plot(dens_cfd[i, :], y_cfd_mm[i, :], "b-", lw=1.0, label="CFD++")
    ax.plot(dens_sim_all[i, :], y_sim_mm[i, :], "r--", lw=1.2, label="Similarity")
    ax.set_ylim(0.0, y_max_mm[i])
    ax.set_title(f"x = {x_i:.2f} m")
    ax.grid(True, ls=":", lw=0.5)

axes[0].set_ylabel("y (mm)")
axes[n_x // 2].set_xlabel(r"$\rho/\rho_e$", labelpad=1)
fig.tight_layout(rect=[0, 0, 1, 0.9])
lines, labels = axes[0].get_legend_handles_labels()
fig.legend(lines, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.99), frameon=False)
filename = "similarity_vs_cfd_dens.png"
fig.savefig(figures_dir / filename, dpi=150, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)
