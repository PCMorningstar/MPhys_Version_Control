import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs jet multiplicity
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [jet multiplicity, fraction, error]
# -----------------------------
bb_data = np.array([
    [2, 0.632237, 0.000503],
    [3, 0.380404, 0.000561],
    [4, 0.261405, 0.000694],
    [5, 0.193209, 0.000950],
    [6, 0.154787, 0.001420],
    [7, 0.124613, 0.002246],
    [8, 0.102189, 0.003705],
    [9, 0.086662, 0.006293],
    [10, 0.080650, 0.011592],
], dtype=float)

bc_data = np.array([
    [2, 0.016645, 0.000133],
    [3, 0.024719, 0.000179],
    [4, 0.027798, 0.000259],
    [5, 0.029951, 0.000409],
    [6, 0.030518, 0.000677],
    [7, 0.028275, 0.001125],
    [8, 0.030139, 0.002088],
    [9, 0.030959, 0.003900],
    [10, 0.038840, 0.008072],
], dtype=float)

bl_data = np.array([
    [2, 0.321095, 0.000487],
    [3, 0.524116, 0.000577],
    [4, 0.562633, 0.000783],
    [5, 0.556167, 0.001193],
    [6, 0.534101, 0.001958],
    [7, 0.508007, 0.003389],
    [8, 0.485846, 0.006097],
    [9, 0.476740, 0.011236],
    [10, 0.471477, 0.021203],
], dtype=float)

cc_data = np.array([
    [2, 0.000247, 0.000017],
    [3, 0.000502, 0.000026],
    [4, 0.001086, 0.000052],
    [5, 0.001564, 0.000095],
    [6, 0.002124, 0.000178],
    [7, 0.002241, 0.000319],
    [8, 0.002019, 0.000522],
    [9, 0.004370, 0.001465],
    [10, 0.005030, 0.002921],
], dtype=float)

cl_data = np.array([
    [2, 0.003105, 0.000058],
    [3, 0.006553, 0.000093],
    [4, 0.013161, 0.000180],
    [5, 0.019848, 0.000334],
    [6, 0.027065, 0.000637],
    [7, 0.031147, 0.001175],
    [8, 0.036895, 0.002277],
    [9, 0.040330, 0.004474],
    [10, 0.035910, 0.007766],
], dtype=float)

ll_data = np.array([
    [2, 0.026671, 0.000169],
    [3, 0.063706, 0.000282],
    [4, 0.133917, 0.000538],
    [5, 0.199262, 0.000959],
    [6, 0.251404, 0.001700],
    [7, 0.305718, 0.003122],
    [8, 0.342911, 0.005802],
    [9, 0.360939, 0.010776],
    [10, 0.368093, 0.020552],
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
jet_bb, f_bb, e_bb = unpack(bb_data)
jet_bc, f_bc, e_bc = unpack(bc_data)
jet_bl, f_bl, e_bl = unpack(bl_data)
jet_cc, f_cc, e_cc = unpack(cc_data)
jet_cl, f_cl, e_cl = unpack(cl_data)
jet_ll, f_ll, e_ll = unpack(ll_data)

plot_edges = centres_to_edges(jet_bb, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, jet_bb, f_bb, e_bb, "{b, b}", "{b, b}")
plot_with_errorbars(ax, jet_bc, f_bc, e_bc, "{b, c}", "{b, c}")
plot_with_errorbars(ax, jet_bl, f_bl, e_bl, "{b, l}", "{b, l}")
plot_with_errorbars(ax, jet_cc, f_cc, e_cc, "{c, c}", "{c, c}")
plot_with_errorbars(ax, jet_cl, f_cl, e_cl, "{c, l}", "{c, l}")
plot_with_errorbars(ax, jet_ll, f_ll, e_ll, "{l, l}", "{l, l}")

ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs Jet Multiplicity")

ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.05)

# Label centres with the jet multiplicity values
ax.set_xticks(jet_bb)
ax.set_xticklabels([str(int(x)) for x in jet_bb])
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
plt.savefig("chi2_pair_flavour_vs_jet_multiplicity.png")
plt.show()
plt.close()

# -----------------------------
# Stacked bar chart
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(jet_bb, dtype=float)

stack_order = [
    ("{b, b}", f_bb, "{b, b}"),
    ("{b, c}", f_bc, "{b, c}"),
    ("{b, l}", f_bl, "{b, l}"),
    ("{c, c}", f_cc, "{c, c}"),
    ("{c, l}", f_cl, "{c, l}"),
    ("{l, l}", f_ll, "{l, l}"),
]

for key, vals, label in stack_order:
    ax.bar(
        jet_bb,
        vals,
        bottom=bottom,
        width=BIN_WIDTH,
        align="center",
        color=COLORS[key],
        edgecolor="white",
        linewidth=0.7,
        label=label
    )
    bottom += vals

ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs Jet Multiplicity")

ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.0)

# Label centres with the jet multiplicity values
ax.set_xticks(jet_bb)
ax.set_xticklabels([str(int(x)) for x in jet_bb])
ax.minorticks_on()

ax.grid(True, axis="y", which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, axis="y", which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    borderaxespad=0.0
)

plt.tight_layout()
plt.savefig("chi2_pair_flavour_stacked_vs_jet_multiplicity.png")
plt.show()
plt.close()