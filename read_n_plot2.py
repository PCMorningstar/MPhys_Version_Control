import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style chi2 pair flavour composition
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data
# -----------------------------
bb_data = np.array([
    [2, 0.632237, 0.000830],
    [3, 0.377607, 0.000711],
    [4, 0.276537, 0.000831],
    [5, 0.210291, 0.001105],
    [6, 0.172765, 0.001630],
    [7, 0.148686, 0.002628],
    [8, 0.104100, 0.003933],
    [9, 0.105969, 0.007385],
    [10, 0.095458, 0.013086],
], dtype=float)

bc_data = np.array([
    [2, 0.016645, 0.000135],
    [3, 0.024797, 0.000182],
    [4, 0.027431, 0.000261],
    [5, 0.029533, 0.000412],
    [6, 0.030127, 0.000683],
    [7, 0.030210, 0.001172],
    [8, 0.026998, 0.001994],
    [9, 0.025652, 0.003561],
    [10, 0.026620, 0.006609],
], dtype=float)

bl_data = np.array([
    [2, 0.321095, 0.000592],
    [3, 0.519539, 0.000832],
    [4, 0.544690, 0.001165],
    [5, 0.542386, 0.001767],
    [6, 0.528069, 0.002855],
    [7, 0.510224, 0.004837],
    [8, 0.464804, 0.008297],
    [9, 0.458021, 0.015147],
    [10, 0.450882, 0.028541],
], dtype=float)

cc_data = np.array([
    [2, 0.000247, 0.000017],
    [3, 0.000499, 0.000026],
    [4, 0.001078, 0.000051],
    [5, 0.001561, 0.000095],
    [6, 0.002109, 0.000179],
    [7, 0.002377, 0.000325],
    [8, 0.002345, 0.000575],
    [9, 0.003035, 0.001193],
    [10, 0.004630, 0.002830],
], dtype=float)

cl_data = np.array([
    [2, 0.003105, 0.000058],
    [3, 0.007306, 0.000099],
    [4, 0.013346, 0.000182],
    [5, 0.020173, 0.000340],
    [6, 0.026439, 0.000637],
    [7, 0.030161, 0.001180],
    [8, 0.036779, 0.002297],
    [9, 0.048475, 0.004985],
    [10, 0.044196, 0.008818],
], dtype=float)

ll_data = np.array([
    [2, 0.026671, 0.000171],
    [3, 0.070252, 0.000306],
    [4, 0.136918, 0.000585],
    [5, 0.196056, 0.001064],
    [6, 0.240492, 0.001924],
    [7, 0.278342, 0.003574],
    [8, 0.364975, 0.007415],
    [9, 0.358848, 0.013520],
    [10, 0.378214, 0.026266],
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

# Okabe-Ito colour-blind-safe palette
COLORS = {
    "{b, b}": "#0072B2",   # blue
    "{b, c}": "#56B4E9",   # sky blue
    "{b, l}": "#E69F00",   # orange
    "{c, c}": "#CC79A7",   # reddish purple
    "{c, l}": "#F0E442",   # yellow
    "{l, l}": "#000000",   # black
}

MARKERS = {
    "{b, b}": "o",
    "{b, c}": "s",
    "{b, l}": "^",
    "{c, c}": "D",
    "{c, l}": "v",
    "{l, l}": "x",
}

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def plot_with_errorbars(ax, x, y, yerr, key, label):
    ax.step(
        x, y, where="mid",
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
nj_bb, frac_bb, err_bb = unpack(bb_data)
nj_bc, frac_bc, err_bc = unpack(bc_data)
nj_bl, frac_bl, err_bl = unpack(bl_data)
nj_cc, frac_cc, err_cc = unpack(cc_data)
nj_cl, frac_cl, err_cl = unpack(cl_data)
nj_ll, frac_ll, err_ll = unpack(ll_data)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, nj_bb, frac_bb, err_bb, "{b, b}", r"{b, b}")
plot_with_errorbars(ax, nj_bc, frac_bc, err_bc, "{b, c}", r"{b, c}")
plot_with_errorbars(ax, nj_bl, frac_bl, err_bl, "{b, l}", r"{b, l}")
plot_with_errorbars(ax, nj_cc, frac_cc, err_cc, "{c, c}", r"{c, c}")
plot_with_errorbars(ax, nj_cl, frac_cl, err_cl, "{c, l}", r"{c, l}")
plot_with_errorbars(ax, nj_ll, frac_ll, err_ll, "{l, l}", r"{l, l}")

# Labels and ranges
ax.set_xlabel(r"Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs. Jet Multiplicity")

ax.set_xlim(1.5, 10.5)
ax.set_ylim(0.0, 0.70)

# Ticks
ax.set_xticks(np.arange(2, 11, 1))
ax.minorticks_on()

# Grid
ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

# Legend
ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    borderaxespad=0.0,
    handlelength=1.8
)

plt.tight_layout()
plt.savefig("chi2_pair_flavour_composition_vs_njets_cb_safe.png")
plt.show()
plt.close()

# -----------------------------
# Optional: stacked plot
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(nj_bb, dtype=float)

stack_order = [
    ("{b, b}", frac_bb, r"{b, b}"),
    ("{b, c}", frac_bc, r"{b, c}"),
    ("{b, l}", frac_bl, r"{b, l}"),
    ("{c, c}", frac_cc, r"{c, c}"),
    ("{c, l}", frac_cl, r"{c, l}"),
    ("{l, l}", frac_ll, r"{l, l}"),
]

for key, vals, label in stack_order:
    ax.bar(
        nj_bb, vals, bottom=bottom, width=0.72,
        color=COLORS[key], edgecolor="white", linewidth=0.7,
        label=label
    )
    bottom += vals

ax.set_xlabel(r"Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition Stacked vs. Jet Multiplicity")

ax.set_xlim(1.4, 10.6)
ax.set_ylim(0.0, 1.0)
ax.set_xticks(np.arange(2, 11, 1))
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
plt.savefig("chi2_pair_flavour_composition_stacked_cb_safe.png")
plt.show()
plt.close()