import numpy as np
import matplotlib.pyplot as plt

b_yield = np.array([
    [-0.25, 1317.627000, 10.969360],
    [ 0.25,  375.179900,  4.869223],
    [ 0.75, 2340.021000, 16.836670],
    [ 1.25, 4546.752000, 29.105130],
    [ 1.75, 6106.936000, 37.677060],
    [ 2.25, 6272.031000, 38.584250],
    [ 2.75, 5250.030000, 32.978610],
    [ 3.25, 3704.564000, 24.471960],
    [ 3.75, 2186.996000, 15.982940],
    [ 4.25, 1120.495000,  9.773727],
    [ 4.75,  486.522700,  5.706070],
    [ 5.25,  183.900100,  3.266412],
], dtype=float)

nonb_yield = np.array([
    [-0.25, 2.731931, 0.395610],
    [ 0.25, 0.864516, 0.211950],
    [ 0.75, 4.811591, 0.530285],
    [ 1.25, 6.760420, 0.638530],
    [ 1.75, 5.244540, 0.559666],
    [ 2.25, 3.106287, 0.434866],
    [ 2.75, 1.069912, 0.259535],
    [ 3.25, 0.947947, 0.233721],
    [ 3.75, 0.558168, 0.171913],
    [ 4.25, 0.138321, 0.088186],
    [ 4.75, 0.066051, 0.066051],
    [ 5.25, 0.153257, 0.088508],
], dtype=float)

# =================================================
# Total MC = b + non-b
# =================================================
total_yield_mc = np.array([
    [-0.25, 1320.359000, 10.976500],
    [ 0.25,  376.044400,  4.873833],
    [ 0.75, 2344.833000, 16.845020],
    [ 1.25, 4553.513000, 29.112130],
    [ 1.75, 6112.180000, 37.681220],
    [ 2.25, 6275.137000, 38.586700],
    [ 2.75, 5251.100000, 32.979630],
    [ 3.25, 3705.512000, 24.473080],
    [ 3.75, 2187.555000, 15.983870],
    [ 4.25, 1120.634000,  9.774124],
    [ 4.75,  486.588800,  5.706452],
    [ 5.25,  184.053300,  3.267611],
], dtype=float)

# =================================================
# 2022 real data
# =================================================
data_yield = np.array([
    [-0.25, 1291.000000, 35.930490],
    [ 0.25,  415.000000, 20.371550],
    [ 0.75, 2463.000000, 49.628620],
    [ 1.25, 4650.000000, 68.190910],
    [ 1.75, 6209.000000, 78.797210],
    [ 2.25, 6016.000000, 77.562880],
    [ 2.75, 5307.000000, 72.849160],
    [ 3.25, 3613.000000, 60.108240],
    [ 3.75, 2211.000000, 47.021270],
    [ 4.25, 1145.000000, 33.837850],
    [ 4.75,  458.000000, 21.400930],
    [ 5.25,  171.000000, 13.076700],
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
ax.set_title(r"HyPER - Probe-Jet Flavour Composition - GN2 WP 65%")
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

plt.savefig("DatavMC_probe_jet_flavour_composition_GN2_WP_65_Hyper.png")
plt.show()
plt.close()