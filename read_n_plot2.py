import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style HyPER pair flavour composition vs jet multiplicity
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [Njets, fraction, error]
# -----------------------------
bb_data = np.array([
    [2, 0.997209, 0.001328],
    [3, 0.764500, 0.001239],
    [4, 0.634062, 0.001533],
    [5, 0.554765, 0.002177],
    [6, 0.500675, 0.003402],
    [7, 0.462214, 0.005707],
    [8, 0.438651, 0.010156],
    [9, 0.431781, 0.019102],
    [10, 0.413041, 0.034740],
], dtype=float)

bc_data = np.array([
    [2, 0.001042, 0.000043],
    [3, 0.011336, 0.000150],
    [4, 0.018333, 0.000259],
    [5, 0.022765, 0.000438],
    [6, 0.027239, 0.000788],
    [7, 0.028291, 0.001408],
    [8, 0.032968, 0.002749],
    [9, 0.033283, 0.005312],
    [10, 0.017344, 0.007124],
], dtype=float)

bl_data = np.array([
    [2, 0.001742, 0.000055],
    [3, 0.223916, 0.000669],
    [4, 0.332386, 0.001106],
    [5, 0.388113, 0.001817],
    [6, 0.420516, 0.003094],
    [7, 0.448161, 0.005594],
    [8, 0.454463, 0.010260],
    [9, 0.442593, 0.019138],
    [10, 0.475488, 0.037732],
], dtype=float)

cc_data = np.array([
    [2, 0.000002, 0.000002],
    [3, 0.000012, 0.000005],
    [4, 0.000095, 0.000018],
    [5, 0.000339, 0.000053],
    [6, 0.000394, 0.000097],
    [7, 0.000421, 0.000164],
    [8, 0.000506, 0.000311],
    [9, 0.000784, 0.000784],
    [10, 0.000000, 0.000000],
], dtype=float)

cl_data = np.array([
    [2, 0.000004, 0.000002],
    [3, 0.000112, 0.000015],
    [4, 0.001522, 0.000075],
    [5, 0.003294, 0.000167],
    [6, 0.005069, 0.000338],
    [7, 0.006947, 0.000687],
    [8, 0.006961, 0.001243],
    [9, 0.008240, 0.002536],
    [10, 0.023177, 0.008273],
], dtype=float)

ll_data = np.array([
    [2, 0.000001, 0.000001],
    [3, 0.000153, 0.000017],
    [4, 0.013600, 0.000223],
    [5, 0.030722, 0.000509],
    [6, 0.046106, 0.001023],
    [7, 0.053968, 0.001935],
    [8, 0.066449, 0.003925],
    [9, 0.083319, 0.008087],
    [10, 0.070950, 0.014442],
], dtype=float)

# -----------------------------
# Style
# -----------------------------
plt.rcParams.update({
    "figure.figsize": (9.0, 6.5),
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
    "{b, b}": "#0072B2",
    "{b, c}": "#56B4E9",
    "{b, l}": "#E69F00",
    "{c, c}": "#CC79A7",
    "{c, l}": "#F0E442",
    "{l, l}": "#000000",
}

MARKERS = {
    "{b, b}": "o",
    "{b, c}": "s",
    "{b, l}": "^",
    "{c, c}": "D",
    "{c, l}": "v",
    "{l, l}": "x",
}

BIN_WIDTH = 1.0

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def centres_to_edges(x, width=1.0):
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
        zorder=2
    )

    ax.errorbar(
        x, y, yerr=yerr,
        fmt=MARKERS[key],
        color=COLORS[key],
        markerfacecolor="white" if MARKERS[key] != "x" else COLORS[key],
        markeredgecolor=COLORS[key],
        markersize=6,
        markeredgewidth=1.2,
        ecolor=COLORS[key],
        elinewidth=1.2,
        capsize=3,
        linewidth=0,
        label=label,
        zorder=3
    )

# -----------------------------
# Unpack
# -----------------------------
n_bb, f_bb, e_bb = unpack(bb_data)
n_bc, f_bc, e_bc = unpack(bc_data)
n_bl, f_bl, e_bl = unpack(bl_data)
n_cc, f_cc, e_cc = unpack(cc_data)
n_cl, f_cl, e_cl = unpack(cl_data)
n_ll, f_ll, e_ll = unpack(ll_data)

edges = centres_to_edges(n_bb, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, n_bb, f_bb, e_bb, "{b, b}", "{b, b}")
plot_with_errorbars(ax, n_bc, f_bc, e_bc, "{b, c}", "{b, c}")
plot_with_errorbars(ax, n_bl, f_bl, e_bl, "{b, l}", "{b, l}")
plot_with_errorbars(ax, n_cc, f_cc, e_cc, "{c, c}", "{c, c}")
plot_with_errorbars(ax, n_cl, f_cl, e_cl, "{c, l}", "{c, l}")
plot_with_errorbars(ax, n_ll, f_ll, e_ll, "{l, l}", "{l, l}")

ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title("HyPER Pair Flavour Composition vs Jet Multiplicity")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.05)

ax.set_xticks(np.arange(2, 11, 1))
ax.minorticks_on()

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    borderaxespad=0.0,
    handlelength=1.8
)

plt.tight_layout()
plt.savefig("hyper_pair_flavour_vs_njets.png")
plt.show()
