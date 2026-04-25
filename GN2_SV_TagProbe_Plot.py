import numpy as np
import matplotlib.pyplot as plt

# ================================================= 
# b - Events
# ================================================= 
b_yield = np.array([
    [-0.25, 1795628.0000000000, 12583.4765620000],
    [ 0.25,  243101.5781250000,  4637.9545900000],
    [ 0.75, 1463879.3750000000, 11385.5458980000],
    [ 1.25, 3529116.2500000000, 17687.3007810000],
    [ 1.75, 6532346.5000000000, 24092.0937500000],
    [ 2.25, 8919209.0000000000, 28161.1308590000],
    [ 2.75, 9575216.0000000000, 29190.3222660000],
    [ 3.25, 8117354.0000000000, 26902.5937500000],
    [ 3.75, 5519906.5000000000, 22183.9941410000],
    [ 4.25, 3124403.2500000000, 16690.1562500000],
    [ 4.75, 1476329.3750000000, 11491.6513670000],
    [ 5.25,  616536.5000000000,  7434.7485350000],
], dtype=float)

# ================================================= 
# non-b - Events
# ================================================= 
nonb_yield = np.array([
    [-0.25,  8494.6474610000, 854.3366090000],
    [ 0.25,   877.1328120000, 268.3302310000],
    [ 0.75,  5185.2763670000, 682.6427610000],
    [ 1.25,  7841.1118160000, 816.2751460000],
    [ 1.75, 10966.6796880000, 966.4057010000],
    [ 2.25, 11628.7626950000, 1005.6610110000],
    [ 2.75,  9859.1933590000, 930.8064580000],
    [ 3.25,  9918.4355470000, 937.4350590000],
    [ 3.75,  4363.7104490000, 636.9465940000],
    [ 4.25,  1471.1967770000, 341.1235660000],
    [ 4.75,   227.8755800000, 139.2441710000],
    [ 5.25,   836.2100830000, 271.5788570000],
], dtype=float)

# ================================================= 
# b/non-b ratio
# ================================================= 
b_to_nonb_ratio_data = np.array([
    [-0.25, 211.361267, 21.300404],
    [ 0.25, 277.123807, 84.969193],
    [ 0.75, 282.297119, 38.221882],
    [ 1.25, 450.121643, 47.512066],
    [ 1.75, 595.503418, 57.787998],
    [ 2.25, 766.886353, 72.370926],
    [ 2.75, 971.183594, 93.372314],
    [ 3.25, 818.354004, 79.820076],
    [ 3.75, 1264.610474, 184.039581],
    [ 4.25, 2123.793701, 507.004242],
    [ 4.75, 6480.814941, 3966.520020],
    [ 5.25, 737.340454, 241.010574],
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

MARKERS = {
    "b": "o",
    "nb": "s",
}

BIN_WIDTH = 0.5

# =================================================
# Helpers
# =================================================
def centres_to_edges(x, width=BIN_WIDTH):
    x = np.asarray(x, dtype=float)
    return np.concatenate(([x[0] - width / 2.0], x + width / 2.0))

def plot_with_errorbars(ax, x, y, yerr, key, label):
    edges = centres_to_edges(x, BIN_WIDTH)

    ax.step(
        edges,
        np.r_[y, y[-1]],
        where="post",
        color=COLORS[key],
        linewidth=2.0,
        zorder=2,
    )

    ax.errorbar(
        x, y, yerr=yerr,
        fmt=MARKERS[key],
        color=COLORS[key],
        markerfacecolor="white",
        markeredgecolor=COLORS[key],
        markersize=6,
        markeredgewidth=1.2,
        ecolor=COLORS[key],
        elinewidth=1.2,
        capsize=3,
        linestyle="none",
        label=label,
        zorder=3,
    )

# =================================================
# Extract
# =================================================
pt = b_yield[:, 0]

b = b_yield[:, 1]
nb = nonb_yield[:, 1]

b_err = b_yield[:, 2]
nb_err = nonb_yield[:, 2]

ratio = b_to_nonb_ratio_data[:, 1]
ratio_err = b_to_nonb_ratio_data[:, 2]

plot_edges = centres_to_edges(pt, BIN_WIDTH)
major_ticks = np.arange(-0.5, 5.5 + 0.5, 0.5)
minor_ticks = np.arange(-0.5, 5.5 + BIN_WIDTH, BIN_WIDTH)

# =================================================
# Figure with ratio panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: yields
# =================================================
plot_with_errorbars(ax, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, pt, nb, nb_err, "nb", "non-b")

ax.set_ylabel("Events")
ax.set_title(r"Probe-Jet Flavour Composition - GN2 WP 65%")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b/non-b ratio
# =================================================
rax.step(
    plot_edges,
    np.r_[ratio, ratio[-1]],
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

rax.errorbar(
    pt,
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

rax.set_ylabel("b/non-b")
rax.set_xlabel(r"SV mass [GeV]")
rax.set_ylim(-1000.0, 10000.0)
rax.set_xlim(plot_edges[0], plot_edges[-1])

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{x:.1f}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("probe_jet_flavour_composition_GN2_WP_65.png")
plt.show()
plt.close()