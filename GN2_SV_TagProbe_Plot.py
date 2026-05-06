import numpy as np
import matplotlib.pyplot as plt

# =================================================
# b - Events - corrected
# =================================================
b_yield = np.array([
    [-0.25,   40.829327,   1.556309],
    [ 0.25,   26.749590,   1.253481],
    [ 0.75,  268.019135,   4.004860],
    [ 1.25, 1010.403015,   7.744958],
    [ 1.75, 2304.259277,  11.712639],
    [ 2.25, 3504.570312,  14.448512],
    [ 2.75, 3945.906494,  15.336808],
    [ 3.25, 3418.974854,  14.287642],
    [ 3.75, 2336.287598,  11.808280],
    [ 4.25, 1326.139038,   8.893599],
    [ 4.75,  626.579956,   6.126587],
    [ 5.25,  257.413513,   3.932517],
], dtype=float)

# =================================================
# non-b - Events - corrected
# =================================================
nonb_yield = np.array([
    [-0.25, 0.130102, 0.078837],
    [ 0.25, 0.035520, 0.035520],
    [ 0.75, 0.187201, 0.108833],
    [ 1.25, 0.594767, 0.182868],
    [ 1.75, 0.778982, 0.208000],
    [ 2.25, 1.379906, 0.286233],
    [ 2.75, 0.837302, 0.226112],
    [ 3.25, 1.499193, 0.301050],
    [ 3.75, 0.583836, 0.196469],
    [ 4.25, 0.163026, 0.095476],
    [ 4.75, 0.000000, 0.000000],
    [ 5.25, 0.053373, 0.053373],
], dtype=float)

# =================================================
# Total MC = b + non-b
# =================================================
total_yield_mc = np.array([
    [-0.25,   40.959429,   1.558305],
    [ 0.25,   26.785110,   1.253984],
    [ 0.75,  268.206336,   4.006339],
    [ 1.25, 1010.997782,   7.747118],
    [ 1.75, 2305.038259,  11.714485],
    [ 2.25, 3505.950218,  14.451348],
    [ 2.75, 3946.743796,  15.338475],
    [ 3.25, 3420.474047,  14.290813],
    [ 3.75, 2336.871434,  11.809915],
    [ 4.25, 1326.302064,   8.894112],
    [ 4.75,  626.579956,   6.126587],
    [ 5.25,  257.466886,   3.932879],
], dtype=float)

# =================================================
# 2022 real data
# =================================================
data_yield = np.array([
    [-0.25,   32.000000,   5.656854],
    [ 0.25,   26.000000,   5.099020],
    [ 0.75,  243.000000,  15.588457],
    [ 1.25,  936.000000,  30.594117],
    [ 1.75, 2064.000000,  45.431267],
    [ 2.25, 3056.000000,  55.281100],
    [ 2.75, 3388.000000,  58.206529],
    [ 3.25, 2898.000000,  53.833075],
    [ 3.75, 2070.000000,  45.497253],
    [ 4.25, 1139.000000,  33.749074],
    [ 4.75,  560.000000,  23.664319],
    [ 5.25,  236.000000,  15.362291],
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


ax.set_ylabel("Weighted Events")
ax.set_title(r"Probe-Jet Flavour Composition - GN2 WP 65%")
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