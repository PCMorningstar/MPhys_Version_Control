import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs leading SV mass
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [SV mass centre, fraction, error]
# -----------------------------
bb_data = np.array([
    [0.25, 0.229096, 0.001234],
    [0.75, 0.482179, 0.001994],
    [1.25, 0.582427, 0.001576],
    [1.75, 0.656961, 0.001325],
    [2.25, 0.710171, 0.001208],
    [2.75, 0.744505, 0.001195],
    [3.25, 0.764152, 0.001310],
    [3.75, 0.772758, 0.001599],
    [4.25, 0.772337, 0.002143],
    [4.75, 0.752568, 0.003192],
    [5.25, 0.694737, 0.005118],
], dtype=float)

bc_data = np.array([
    [0.25, 0.021591, 0.000428],
    [0.75, 0.026764, 0.000644],
    [1.25, 0.021913, 0.000466],
    [1.75, 0.018010, 0.000372],
    [2.25, 0.014710, 0.000320],
    [2.75, 0.012982, 0.000308],
    [3.25, 0.011855, 0.000332],
    [3.75, 0.011317, 0.000401],
    [4.25, 0.011512, 0.000548],
    [4.75, 0.013441, 0.000860],
    [5.25, 0.015141, 0.001326],
], dtype=float)

bl_data = np.array([
    [0.25, 0.527682, 0.001468],
    [0.75, 0.473158, 0.001993],
    [1.25, 0.390409, 0.001559],
    [1.75, 0.323322, 0.001305],
    [2.25, 0.274399, 0.001188],
    [2.75, 0.242138, 0.001174],
    [3.25, 0.223570, 0.001286],
    [3.75, 0.215577, 0.001570],
    [4.25, 0.215941, 0.002102],
    [4.75, 0.233683, 0.003129],
    [5.25, 0.289399, 0.005042],
], dtype=float)

cc_data = np.array([
    [0.25, 0.001083, 0.000097],
    [0.75, 0.000651, 0.000105],
    [1.25, 0.000382, 0.000063],
    [1.75, 0.000122, 0.000030],
    [2.25, 0.000029, 0.000013],
    [2.75, 0.000013, 0.000009],
    [3.25, 0.000000, 0.000000],
    [3.75, 0.000017, 0.000017],
    [4.25, 0.000000, 0.000000],
    [4.75, 0.000000, 0.000000],
    [5.25, 0.000122, 0.000122],
], dtype=float)

cl_data = np.array([
    [0.25, 0.018268, 0.000394],
    [0.75, 0.005360, 0.000297],
    [1.25, 0.002173, 0.000148],
    [1.75, 0.000815, 0.000080],
    [2.25, 0.000229, 0.000040],
    [2.75, 0.000126, 0.000030],
    [3.25, 0.000167, 0.000039],
    [3.75, 0.000105, 0.000040],
    [4.25, 0.000075, 0.000044],
    [4.75, 0.000283, 0.000128],
    [5.25, 0.000338, 0.000201],
], dtype=float)

ll_data = np.array([
    [0.25, 0.201807, 0.001182],
    [0.75, 0.011651, 0.000429],
    [1.25, 0.002615, 0.000162],
    [1.75, 0.000649, 0.000070],
    [2.25, 0.000366, 0.000052],
    [2.75, 0.000108, 0.000028],
    [3.25, 0.000108, 0.000031],
    [3.75, 0.000087, 0.000036],
    [4.25, 0.000044, 0.000031],
    [4.75, 0.000026, 0.000026],
    [5.25, 0.000152, 0.000152],
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

BIN_WIDTH = 0.5

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def centres_to_edges(x, width=0.5):
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
sv_bb, f_bb, e_bb = unpack(bb_data)
sv_bc, f_bc, e_bc = unpack(bc_data)
sv_bl, f_bl, e_bl = unpack(bl_data)
sv_cc, f_cc, e_cc = unpack(cc_data)
sv_cl, f_cl, e_cl = unpack(cl_data)
sv_ll, f_ll, e_ll = unpack(ll_data)

edges = centres_to_edges(sv_bb, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, sv_bb, f_bb, e_bb, "{b, b}", "{b, b}")
plot_with_errorbars(ax, sv_bc, f_bc, e_bc, "{b, c}", "{b, c}")
plot_with_errorbars(ax, sv_bl, f_bl, e_bl, "{b, l}", "{b, l}")
plot_with_errorbars(ax, sv_cc, f_cc, e_cc, "{c, c}", "{c, c}")
plot_with_errorbars(ax, sv_cl, f_cl, e_cl, "{c, l}", "{c, l}")
plot_with_errorbars(ax, sv_ll, f_ll, e_ll, "{l, l}", "{l, l}")

ax.set_xlabel("Leading SV Mass [GeV]")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs Leading SV Mass")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.05)

# Define bin edges explicitly
edges = np.arange(0.0, 5.5 + 0.5, 0.5)  # 0 → 5.5 in steps of 0.5

ax.set_xticks(edges)
ax.set_xticklabels([f"{e:.1f}" for e in edges])
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
plt.savefig("chi2_pair_flavour_vs_sv_mass.png")
plt.show()
plt.close()

# -----------------------------
# Stacked bar chart
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(sv_bb, dtype=float)

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
        sv_bb,
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

ax.set_xlabel("Leading SV Mass [GeV]")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Pair Flavour Composition vs Leading SV Mass")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.0)

# Define bin edges explicitly
edges = np.arange(0.0, 5.5 + 0.5, 0.5)  # 0 → 5.5 in steps of 0.5

ax.set_xticks(edges)
ax.set_xticklabels([f"{e:.1f}" for e in edges])
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
plt.savefig("chi2_pair_flavour_stacked_vs_sv_mass.png")
plt.show()
plt.close()