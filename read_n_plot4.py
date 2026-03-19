import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs leading jet pT
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [leading jet pT, fraction, error]
# -----------------------------
bb_data = np.array([
    [15, 0.312367, 0.003056],
    [45, 0.557767, 0.000993],
    [75, 0.667784, 0.000879],
    [105, 0.685755, 0.001103],
    [135, 0.686308, 0.001544],
    [165, 0.685588, 0.002210],
    [195, 0.661416, 0.003243],
    [225, 0.629952, 0.004697],
    [255, 0.578906, 0.006663],
    [285, 0.544940, 0.009065],
], dtype=float)

bc_data = np.array([
    [15, 0.023622, 0.001003],
    [45, 0.018285, 0.000268],
    [75, 0.015396, 0.000229],
    [105, 0.015220, 0.000291],
    [135, 0.015462, 0.000411],
    [165, 0.015687, 0.000588],
    [195, 0.016986, 0.000893],
    [225, 0.019103, 0.001319],
    [255, 0.022051, 0.001956],
    [285, 0.024910, 0.002760],
], dtype=float)

bl_data = np.array([
    [15, 0.517566, 0.003300],
    [45, 0.379235, 0.000971],
    [75, 0.296025, 0.000852],
    [105, 0.282680, 0.001070],
    [135, 0.281151, 0.001495],
    [165, 0.277972, 0.002133],
    [195, 0.296195, 0.003128],
    [225, 0.315773, 0.004523],
    [255, 0.354548, 0.006455],
    [285, 0.374020, 0.008802],
], dtype=float)

cc_data = np.array([
    [15, 0.001176, 0.000223],
    [45, 0.000309, 0.000036],
    [75, 0.000180, 0.000025],
    [105, 0.000175, 0.000032],
    [135, 0.000113, 0.000036],
    [165, 0.000250, 0.000076],
    [195, 0.000301, 0.000119],
    [225, 0.000637, 0.000235],
    [255, 0.000431, 0.000305],
    [285, 0.000111, 0.000111],
], dtype=float)

cl_data = np.array([
    [15, 0.012749, 0.000745],
    [45, 0.004183, 0.000129],
    [75, 0.002317, 0.000090],
    [105, 0.001979, 0.000106],
    [135, 0.001794, 0.000142],
    [165, 0.002561, 0.000252],
    [195, 0.003563, 0.000416],
    [225, 0.004251, 0.000633],
    [255, 0.005812, 0.001011],
    [285, 0.005743, 0.001367],
], dtype=float)

ll_data = np.array([
    [15, 0.132520, 0.002250],
    [45, 0.040220, 0.000395],
    [75, 0.018298, 0.000251],
    [105, 0.014190, 0.000282],
    [135, 0.015172, 0.000408],
    [165, 0.017942, 0.000630],
    [195, 0.021538, 0.000990],
    [225, 0.030284, 0.001662],
    [255, 0.038252, 0.002547],
    [285, 0.050276, 0.003926],
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

BIN_WIDTH = 30.0

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

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

ax.set_xlabel("Leading Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs. Leading Jet $p_T$")

ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.05)

# Define bin edges explicitly
edges = np.arange(0.0, 300.0 + BIN_WIDTH, BIN_WIDTH)  # 0 → 5.5 in steps of BIN_WIDTH

# Label centres with the leading jet pT values
ax.set_xticks(edges)
ax.set_xticklabels([f"{e:.0f}" for e in edges])
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
plt.savefig("chi2_pair_flavour_vs_leading_jet_pT.png")
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

ax.set_xlabel("Leading Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs. Leading Jet $p_T$")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.0)

# Label centres with the leading jet pT values
ax.set_xticks(edges)
ax.set_xticklabels([f"{e:.0f}" for e in edges])
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
plt.savefig("chi2_pair_flavour_stacked_vs_leading_jet_pT.png")
plt.show()
plt.close()