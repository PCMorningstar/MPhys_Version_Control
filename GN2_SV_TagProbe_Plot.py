import numpy as np
import matplotlib.pyplot as plt

# =================================================
# b - Events - corrected
# =================================================
b_yield = np.array([
    [-0.25,   41.007889,   1.560346],
    [ 0.25,   26.688442,   1.252819],
    [ 0.75,  268.798218,   4.010045],
    [ 1.25, 1010.887878,   7.747545],
    [ 1.75, 2300.191162,  11.700373],
    [ 2.25, 3494.957031,  14.429359],
    [ 2.75, 3931.815430,  15.309745],
    [ 3.25, 3409.349609,  14.267761],
    [ 3.75, 2327.449219,  11.786798],
    [ 4.25, 1322.819702,   8.882799],
    [ 4.75,  625.590210,   6.121982],
    [ 5.25,  256.905914,   3.928616],
], dtype=float)

# =================================================
# non-b - Events - corrected
# =================================================
nonb_yield = np.array([
    [-0.25, 0.077817, 0.059005],
    [ 0.25, 0.000000, 0.000000],
    [ 0.75, 0.072826, 0.072826],
    [ 1.25, 0.570535, 0.181255],
    [ 1.75, 0.574948, 0.181241],
    [ 2.25, 1.089280, 0.253113],
    [ 2.75, 0.712352, 0.207429],
    [ 3.25, 1.403130, 0.293136],
    [ 3.75, 0.450786, 0.172465],
    [ 4.25, 0.099179, 0.070988],
    [ 4.75, 0.000000, 0.000000],
    [ 5.25, 0.053373, 0.053373],
], dtype=float)

# =================================================
# Total MC = b + non-b
# =================================================
total_yield_mc = np.array([
    [-0.25,   41.085706,   1.561462],
    [ 0.25,   26.688442,   1.252819],
    [ 0.75,  268.871044,   4.010706],
    [ 1.25, 1011.458413,   7.749663],
    [ 1.75, 2300.766110,  11.701777],
    [ 2.25, 3496.046311,  14.431578],
    [ 2.75, 3932.527782,  15.311149],
    [ 3.25, 3410.752739,  14.270772],
    [ 3.75, 2327.900005,  11.788060],
    [ 4.25, 1322.918881,   8.883083],
    [ 4.75,  625.590210,   6.121982],
    [ 5.25,  256.959287,   3.928978],
], dtype=float)

# =================================================
# 2022 real data
# =================================================
data_yield = np.array([
    [-0.25,   32.000000,   5.656854],
    [ 0.25,   29.000000,   5.385165],
    [ 0.75,  243.000000,  15.588457],
    [ 1.25,  936.000000,  30.594117],
    [ 1.75, 2059.000000,  45.376205],
    [ 2.25, 3046.000000,  55.190579],
    [ 2.75, 3379.000000,  58.129167],
    [ 3.25, 2884.000000,  53.702886],
    [ 3.75, 2064.000000,  45.431267],
    [ 4.25, 1137.000000,  33.719431],
    [ 4.75,  556.000000,  23.579652],
    [ 5.25,  237.000000,  15.394804],
], dtype=float)

# =================================================
# Style
# =================================================
plt.rcParams.update({
    "figure.figsize": (9.0, 7.0),
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 17,
    "legend.fontsize": 12,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

COLORS = {
    "b": "#0072B2",
    "nb": "#E69F00",
}

# =================================================
# Variable-width SV-mass binning
# =================================================
bin_edges = np.array([
    -0.5, 0.0, 0.5,
     1.0, 1.5, 2.0,
     2.5, 3.0, 3.5,
     4.0, 4.5, 5.0,
     5.5
], dtype=float)

bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

major_ticks = np.arange(-0.5, 5.5 + 0.5, 0.5)
minor_ticks = np.arange(-0.5, 5.5 + 0.25, 0.25)

# =================================================
# Extract
# =================================================
b = b_yield[:, 1]
nb = nonb_yield[:, 1]

data = data_yield[:, 1]
data_err = data_yield[:, 2]

mc = total_yield_mc[:, 1]
mc_err = total_yield_mc[:, 2]

total = b + nb

# Gaussian propagation for Data / MC
ratio = data / mc
ratio_err = ratio * np.sqrt((data_err / data) ** 2 + (mc_err / mc) ** 2)

# Step arrays
b_step = np.r_[b, b[-1]]
nb_step = np.r_[nb, nb[-1]]
total_step = np.r_[total, total[-1]]
data_step = np.r_[data, data[-1]]
ratio_step = np.r_[ratio, ratio[-1]]

plot_floor = 1e-3

# =================================================
# Figure with ratio panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: data and stacked MC yields
# =================================================
ax.step(
    bin_edges,
    data_step,
    where="post",
    color="black",
    linewidth=1.4,
    zorder=4,
)

ax.errorbar(
    bin_centres,
    data,
    yerr=data_err,
    fmt="o",
    color="black",
    markerfacecolor="white",
    markeredgecolor="black",
    markersize=4,
    markeredgewidth=1.1,
    ecolor="black",
    elinewidth=1.2,
    capsize=3,
    linestyle="none",
    zorder=5,
    label="Data [2022]",
)

# b stacked on top of non-b
ax.fill_between(
    bin_edges,
    np.maximum(nb_step, plot_floor),
    np.maximum(total_step, plot_floor),
    step="post",
    color=COLORS["b"],
    alpha=0.7,
    label="b [MC]",
)

# non-b bottom
ax.fill_between(
    bin_edges,
    plot_floor,
    np.maximum(nb_step, plot_floor),
    step="post",
    color=COLORS["nb"],
    alpha=0.7,
    label="non-b [MC]",
)


ax.set_ylabel("Weighted Yield")
ax.set_title(r"$\chi^2$ - Probe-Jet Flavour Composition - GN2 WP 65%")
ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.set_ylim(plot_floor, 1e4)
ax.set_yscale("log", base=10)
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0.0,
    frameon=False,
)

# =================================================
# Bottom: Data / MC ratio
# =================================================
rax.step(
    bin_edges,
    ratio_step,
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

rax.errorbar(
    bin_centres,
    ratio,
    yerr=ratio_err,
    fmt="o",
    color="black",
    markerfacecolor="white",
    markeredgecolor="black",
    markersize=4,
    markeredgewidth=1.1,
    ecolor="black",
    elinewidth=1.2,
    capsize=3,
    linestyle="none",
    zorder=3,
)

rax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)

rax.set_ylabel("Data / MC")
rax.set_xlabel(r"SV mass [GeV]")
rax.set_ylim(0.5, 1.5)
rax.set_xlim(bin_edges[0], bin_edges[-1])

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{x:.1f}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("DatavMC_probe_jet_flavour_composition_GN2_WP_65.png")
plt.show()
plt.close()