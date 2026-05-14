import numpy as np
import matplotlib.pyplot as plt

################### hyper combined - SV INVARIANT MASS ###################

# =================================================
# b yields
# =================================================
b_yield = np.array([
    [-0.25, 34992.540000, 45.538370],
    [ 0.25,  4505.606000, 16.336590],
    [ 0.75, 19063.130000, 33.630570],
    [ 1.25, 25224.240000, 38.666890],
    [ 1.75, 27951.580000, 40.713420],
    [ 2.25, 26473.230000, 39.626300],
    [ 2.75, 22498.980000, 36.537330],
    [ 3.25, 16444.770000, 31.260760],
    [ 3.75, 10343.720000, 24.789030],
    [ 4.25,  5620.750000, 18.278200],
    [ 4.75,  2690.906000, 12.657100],
    [ 5.25,  1181.594000,  8.385604],
], dtype=float)

# =================================================
# c yields
# =================================================
c_yield = np.array([
    [-0.25, 2225.101000, 11.451620],
    [ 0.25,  111.375500,  2.567530],
    [ 0.75,  343.690400,  4.511803],
    [ 1.25,  269.767300,  3.977071],
    [ 1.75,  153.778800,  3.019043],
    [ 2.25,   51.212160,  1.737043],
    [ 2.75,   29.154790,  1.320740],
    [ 3.25,   15.653320,  0.964116],
    [ 3.75,   10.565430,  0.782011],
    [ 4.25,    6.045898,  0.600458],
    [ 4.75,    5.145020,  0.552910],
    [ 5.25,    4.072510,  0.494523],
], dtype=float)

# =================================================
# light yields
# =================================================
l_yield = np.array([
    [-0.25, 54991.770000, 57.016820],
    [ 0.25,   564.011700,  5.778842],
    [ 0.75,   859.406300,  7.129145],
    [ 1.25,   312.753900,  4.304369],
    [ 1.75,   115.605500,  2.626720],
    [ 2.25,    46.523440,  1.656545],
    [ 2.75,    19.765630,  1.089725],
    [ 3.25,     8.695313,  0.720785],
    [ 3.75,     6.542969,  0.629476],
    [ 4.25,     2.601563,  0.391873],
    [ 4.75,     1.359375,  0.280381],
    [ 5.25,     1.046875,  0.249511],
], dtype=float)

# =================================================
# b fractions
# =================================================
b_frac = np.array([
    [-0.25, 0.380945, 0.000404],
    [ 0.25, 0.869798, 0.001452],
    [ 0.75, 0.939282, 0.000517],
    [ 1.25, 0.976875, 0.000294],
    [ 1.75, 0.991315, 0.000190],
    [ 2.25, 0.998166, 0.000121],
    [ 2.75, 0.997843, 0.000116],
    [ 3.25, 0.998483, 0.000093],
    [ 3.75, 0.998332, 0.000129],
    [ 4.25, 0.998462, 0.000170],
    [ 4.75, 0.997583, 0.000296],
    [ 5.25, 0.995702, 0.000639],
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
    "c": "#56B4E9",
    "l": "#E69F00",
}

MARKERS = {
    "b": "o",
    "c": "s",
    "l": "^",
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
c = c_yield[:, 1]
l = l_yield[:, 1]

b_err = b_yield[:, 2]
c_err = c_yield[:, 2]
l_err = l_yield[:, 2]

frac_b = b_frac[:, 1]
err_b = b_frac[:, 2]

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
# Top: line histogram (yields)
# =================================================
plot_with_errorbars(ax, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, pt, c, c_err, "c", "c")
plot_with_errorbars(ax, pt, l, l_err, "l", "l")

ax.set_ylabel("Weighted Yield")
ax.set_title(r"HyPER Flavour Composition vs. SV Mass - Combined")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(1e-2, 1e6)
ax.set_yscale("log", base=10)
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b fraction (black, step-style)
# =================================================
rax.step(
    plot_edges,
    np.r_[frac_b, frac_b[-1]],
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

rax.errorbar(
    pt,
    frac_b,
    yerr=err_b,
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

rax.set_ylabel("b Fraction")
rax.set_xlabel(r"SV Mass [GeV]")
rax.set_ylim(0.0, 1.1)
rax.set_xlim(plot_edges[0], plot_edges[-1])

# Major ticks every 90 GeV, minor ticks every 30 GeV
rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{float(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("hyper_combined_flavour_vs_SV_invariant_mass.png")
plt.show()
plt.close()