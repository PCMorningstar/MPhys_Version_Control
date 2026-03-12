import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style chi2 pair flavour composition vs jet pT
# Colour-blind safe (Okabe-Ito palette)
# geq 2j only
# =================================================

# -----------------------------
# Hard-coded data: [pT centre, fraction, error]
# -----------------------------
bb_data = np.array([
    [15,  0.312367, 0.003682],
    [45,  0.650092, 0.002173],
    [75,  0.786641, 0.003964],
    [105, 0.822064, 0.007431],
    [135, 0.853557, 0.013067],
    [165, 0.843059, 0.021721],
    [195, 0.871810, 0.036829],
    [225, 0.836328, 0.053169],
    [255, 0.819067, 0.087678],
    [285, 0.905030, 0.150203],
], dtype=float)

bc_data = np.array([
    [15,  0.023622, 0.001015],
    [45,  0.016313, 0.000345],
    [75,  0.012015, 0.000489],
    [105, 0.010380, 0.000825],
    [135, 0.008829, 0.001381],
    [165, 0.010777, 0.002516],
    [195, 0.004210, 0.002435],
    [225, 0.017780, 0.008094],
    [255, 0.005492, 0.005492],
    [285, 0.000000, 0.000000],
], dtype=float)

bl_data = np.array([
    [15,  0.517566, 0.004750],
    [45,  0.304857, 0.001491],
    [75,  0.191773, 0.001959],
    [105, 0.161868, 0.003294],
    [135, 0.132824, 0.005157],
    [165, 0.143274, 0.008985],
    [195, 0.123980, 0.014006],
    [225, 0.145891, 0.022252],
    [255, 0.158288, 0.036980],
    [285, 0.094970, 0.045431],
], dtype=float)

cc_data = np.array([
    [15,  0.001176, 0.000223],
    [45,  0.000245, 0.000045],
    [75,  0.000070, 0.000037],
    [105, 0.000060, 0.000060],
    [135, 0.000000, 0.000000],
    [165, 0.000000, 0.000000],
    [195, 0.000000, 0.000000],
    [225, 0.000000, 0.000000],
    [255, 0.000000, 0.000000],
    [285, 0.000000, 0.000000],
], dtype=float)

cl_data = np.array([
    [15,  0.012749, 0.000750],
    [45,  0.002670, 0.000139],
    [75,  0.001214, 0.000156],
    [105, 0.000913, 0.000246],
    [135, 0.000687, 0.000366],
    [165, 0.000816, 0.000816],
    [195, 0.000000, 0.000000],
    [225, 0.000000, 0.000000],
    [255, 0.000000, 0.000000],
    [285, 0.000000, 0.000000],
], dtype=float)

ll_data = np.array([
    [15,  0.132520, 0.002418],
    [45,  0.025824, 0.000436],
    [75,  0.008288, 0.000408],
    [105, 0.004715, 0.000565],
    [135, 0.004103, 0.000890],
    [165, 0.002074, 0.001058],
    [195, 0.000000, 0.000000],
    [225, 0.000000, 0.000000],
    [255, 0.000000, 0.000000],
    [285, 0.000000, 0.000000],
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

BIN_WIDTH = 30.0  # GeV

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def centres_to_edges(x, width=30.0):
    x = np.asarray(x, dtype=float)
    return np.concatenate(([x[0] - width/2], x + width/2))

def plot_with_errorbars(ax, x, y, yerr, key, label, bin_width=30.0):
    edges = centres_to_edges(x, width=bin_width)

    # Step histogram spanning full bin width
    ax.step(
        edges,
        np.r_[y, y[-1]],
        where="post",
        color=COLORS[key],
        linewidth=2.0,
        zorder=2
    )

    # Marker + uncertainty at bin centre
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
pt_bb, frac_bb, err_bb = unpack(bb_data)
pt_bc, frac_bc, err_bc = unpack(bc_data)
pt_bl, frac_bl, err_bl = unpack(bl_data)
pt_cc, frac_cc, err_cc = unpack(cc_data)
pt_cl, frac_cl, err_cl = unpack(cl_data)
pt_ll, frac_ll, err_ll = unpack(ll_data)

edges = centres_to_edges(pt_bb, width=BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, pt_bb, frac_bb, err_bb, "{b, b}", r"{b, b}", bin_width=BIN_WIDTH)
plot_with_errorbars(ax, pt_bc, frac_bc, err_bc, "{b, c}", r"{b, c}", bin_width=BIN_WIDTH)
plot_with_errorbars(ax, pt_bl, frac_bl, err_bl, "{b, l}", r"{b, l}", bin_width=BIN_WIDTH)
plot_with_errorbars(ax, pt_cc, frac_cc, err_cc, "{c, c}", r"{c, c}", bin_width=BIN_WIDTH)
plot_with_errorbars(ax, pt_cl, frac_cl, err_cl, "{c, l}", r"{c, l}", bin_width=BIN_WIDTH)
plot_with_errorbars(ax, pt_ll, frac_ll, err_ll, "{l, l}", r"{l, l}", bin_width=BIN_WIDTH)

ax.set_xlabel(r"Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction of events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs. Jet $p_T$")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.0)

ax.set_xticks(edges)
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
plt.savefig("chi2_pair_flavour_composition_vs_jetpt_cb_safe.png")
plt.show()
plt.close()

# -----------------------------
# Optional: stacked plot
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(pt_bb, dtype=float)

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
        pt_bb,
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

ax.set_xlabel(r"Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction of events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs. Jet $p_T$")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.0)

# label bin edges only
ax.set_xticks(edges)
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
plt.savefig("chi2_pair_flavour_composition_stacked_jetpt_cb_safe.png")
plt.show()
plt.close()