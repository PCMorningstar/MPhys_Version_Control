import numpy as np
import matplotlib.pyplot as plt

b_yield = np.array([
    [-0.25, 4490764.0, 19990.212891],
    [0.25, 856397.75, 8708.93457],
    [0.75, 4224305.5, 19401.419922],
    [1.25, 7104647.5, 25158.882812],
    [1.75, 9775320.0, 29519.195312],
    [2.25, 11098853.0, 31457.304688],
    [2.75, 10713420.0, 30919.421875],
    [3.25, 8550038.0, 27643.728516],
    [3.75, 5597266.5, 22362.837891],
    [4.25, 3123950.5, 16706.958984],
    [4.75, 1475606.0, 11492.398438],
    [5.25, 631929.25, 7532.525879],
], dtype=float)

c_yield = np.array([
    [-0.25, 192160.21875, 4144.179688],
    [0.25, 21766.25, 1394.843506],
    [0.75, 96408.75, 2962.143555],
    [1.25, 111343.625, 3138.117676],
    [1.75, 114386.90625, 3208.688232],
    [2.25, 95494.210938, 2907.963379],
    [2.75, 85187.4375, 2742.323975],
    [3.25, 57307.796875, 2248.375],
    [3.75, 35959.980469, 1782.501587],
    [4.25, 20715.617188, 1373.643433],
    [4.75, 13019.289062, 1100.987793],
    [5.25, 5561.141602, 693.335388],
], dtype=float)

l_yield = np.array([
    [-0.25, 4660101.0, 20394.138672],
    [0.25, 382883.0625, 5844.226562],
    [0.75, 1448286.375, 11385.543945],
    [1.25, 1777492.625, 12582.201172],
    [1.75, 1927078.25, 13108.765625],
    [2.25, 1778849.25, 12594.480469],
    [2.75, 1483273.5, 11499.725586],
    [3.25, 1085628.875, 9850.158203],
    [3.75, 689172.25, 7857.586914],
    [4.25, 383084.9375, 5843.016113],
    [4.75, 195117.90625, 4183.024902],
    [5.25, 106873.421875, 3083.212402],
], dtype=float)

b_frac = np.array([
    [-0.25, 0.480654, 0.001543],
    [0.25, 0.679116, 0.003922],
    [0.75, 0.732242, 0.001744],
    [1.25, 0.789977, 0.001282],
    [1.75, 0.827240, 0.001039],
    [2.25, 0.855522, 0.000922],
    [2.75, 0.872295, 0.000899],
    [3.25, 0.882086, 0.000979],
    [3.75, 0.885307, 0.001199],
    [4.25, 0.885536, 0.001601],
    [4.75, 0.876384, 0.002404],
    [5.25, 0.848952, 0.003915],
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

ax.set_ylabel("Weighted Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. SV Mass - Top2")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(1e3, 1e8)
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
rax.set_ylim(0.0, 1.0)
rax.set_xlim(plot_edges[0], plot_edges[-1])

# Major ticks every 90 GeV, minor ticks every 30 GeV
rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{float(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("chi2_top2_flavour_vs_SV_invariant_mass.png")
plt.show()
plt.close()