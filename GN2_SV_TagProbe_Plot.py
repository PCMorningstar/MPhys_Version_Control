import numpy as np
import matplotlib.pyplot as plt

b_yield = np.array([
    [-0.25, 1226.807739,  8.499174],
    [ 0.25,  167.533417,  3.145997],
    [ 0.75, 1012.231750,  7.735671],
    [ 1.25, 2440.517578, 12.020829],
    [ 1.75, 4492.460449, 16.329641],
    [ 2.25, 6115.364258, 19.059971],
    [ 2.75, 6570.312012, 19.759647],
    [ 3.25, 5565.834473, 18.206757],
    [ 3.75, 3781.346191, 15.006824],
    [ 4.25, 2141.918457, 11.291224],
    [ 4.75, 1014.655823,  7.788656],
    [ 5.25,  422.212189,  5.028871],
], dtype=float)

nonb_yield = np.array([
    [-0.25, 2.363751, 0.362003],
    [ 0.25, 0.293880, 0.132216],
    [ 0.75, 1.601413, 0.304359],
    [ 1.25, 3.826025, 0.468395],
    [ 1.75, 3.991644, 0.478313],
    [ 2.25, 4.380769, 0.504499],
    [ 2.75, 3.733285, 0.472312],
    [ 3.25, 3.953874, 0.488992],
    [ 3.75, 1.571640, 0.309339],
    [ 4.25, 0.459528, 0.158132],
    [ 4.75, 0.028092, 0.028092],
    [ 5.25, 0.310903, 0.130178],
], dtype=float)

# =================================================
# Total MC = b + non-b
# =================================================
total_yield_mc = np.array([
    [-0.25, 1229.171490,  8.506880],
    [ 0.25,  167.827297,  3.148774],
    [ 0.75, 1013.833163,  7.741656],
    [ 1.25, 2444.343603, 12.029951],
    [ 1.75, 4496.452093, 16.336645],
    [ 2.25, 6119.745027, 19.066647],
    [ 2.75, 6574.045297, 19.765291],
    [ 3.25, 5569.788347, 18.213322],
    [ 3.75, 3782.917831, 15.010012],
    [ 4.25, 2142.377985, 11.292331],
    [ 4.75, 1014.683915,  7.788707],
    [ 5.25,  422.523092,  5.030556],
], dtype=float)

# =================================================
# 2022 real data
# =================================================
data_yield = np.array([
    [-0.25, 1172.000000, 34.234486],
    [ 0.25,  165.000000, 12.845233],
    [ 0.75, 1021.000000, 31.953091],
    [ 1.25, 2333.000000, 48.301139],
    [ 1.75, 4057.000000, 63.694584],
    [ 2.25, 5356.000000, 73.184698],
    [ 2.75, 5684.000000, 75.392307],
    [ 3.25, 4705.000000, 68.593003],
    [ 3.75, 3379.000000, 58.129167],
    [ 4.25, 1834.000000, 42.825226],
    [ 4.75,  897.000000, 29.949958],
    [ 5.25,  386.000000, 19.646883],
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

plot_floor = 1e-2

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
ax.set_ylim(plot_floor, 1e6)
ax.set_yscale("log", base=10)
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()

ax.legend(
    loc="upper right",
    #bbox_to_anchor=(1.02, 0.5),
    #borderaxespad=0.0,
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