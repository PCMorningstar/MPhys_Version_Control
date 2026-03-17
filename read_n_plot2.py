import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style HyPER pair flavour composition vs leading SV mass
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [SV mass centre, fraction, error]
# -----------------------------
bb_data = np.array([
    [0.25, 0.438200, 0.001979],
    [0.75, 0.671267, 0.002781],
    [1.25, 0.751189, 0.002235],
    [1.75, 0.803772, 0.001940],
    [2.25, 0.837975, 0.001838],
    [2.75, 0.859686, 0.001873],
    [3.25, 0.873070, 0.002096],
    [3.75, 0.878860, 0.002569],
    [4.25, 0.882865, 0.003445],
    [4.75, 0.882929, 0.004997],
], dtype=float)

bc_data = np.array([
    [0.25, 0.017915, 0.000397],
    [0.75, 0.017336, 0.000444],
    [1.25, 0.014692, 0.000310],
    [1.75, 0.011792, 0.000233],
    [2.25, 0.009095, 0.000190],
    [2.75, 0.007759, 0.000177],
    [3.25, 0.006986, 0.000187],
    [3.75, 0.006948, 0.000228],
    [4.25, 0.006223, 0.000285],
    [4.75, 0.007519, 0.000455],
], dtype=float)

bl_data = np.array([
    [0.25, 0.448129, 0.001994],
    [0.75, 0.305030, 0.001864],
    [1.25, 0.232326, 0.001236],
    [1.75, 0.183711, 0.000924],
    [2.25, 0.152730, 0.000781],
    [2.75, 0.132488, 0.000731],
    [3.25, 0.119922, 0.000772],
    [3.75, 0.114111, 0.000920],
    [4.25, 0.110866, 0.001215],
    [4.75, 0.109473, 0.001750],
], dtype=float)

cc_data = np.array([
    [0.25, 0.000336, 0.000054],
    [0.75, 0.000201, 0.000047],
    [1.25, 0.000156, 0.000032],
    [1.75, 0.000068, 0.000017],
    [2.25, 0.000005, 0.000004],
    [2.75, 0.000005, 0.000004],
    [3.25, 0.000005, 0.000005],
    [3.75, 0.000007, 0.000007],
    [4.25, 0.000000, 0.000000],
    [4.75, 0.000000, 0.000000],
], dtype=float)

cl_data = np.array([
    [0.25, 0.007378, 0.000254],
    [0.75, 0.001729, 0.000141],
    [1.25, 0.000696, 0.000067],
    [1.75, 0.000388, 0.000042],
    [2.25, 0.000101, 0.000020],
    [2.75, 0.000046, 0.000014],
    [3.25, 0.000007, 0.000005],
    [3.75, 0.000055, 0.000019],
    [4.25, 0.000034, 0.000021],
    [4.75, 0.000053, 0.000038],
], dtype=float)

ll_data = np.array([
    [0.25, 0.088037, 0.000880],
    [0.75, 0.004435, 0.000223],
    [1.25, 0.000939, 0.000079],
    [1.75, 0.000282, 0.000036],
    [2.25, 0.000106, 0.000021],
    [2.75, 0.000021, 0.000010],
    [3.25, 0.000017, 0.000009],
    [3.75, 0.000024, 0.000014],
    [4.25, 0.000013, 0.000013],
    [4.75, 0.000028, 0.000028],
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
ax.set_title("HyPER Pair Flavour Composition vs Leading SV Mass")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.05)

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
plt.savefig("hyper_pair_flavour_vs_sv_mass.png")
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
ax.set_title("HyPER Pair Flavour Composition vs Leading SV Mass")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.0)

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
plt.savefig("hyper_pair_flavour_stacked_vs_sv_mass.png")
plt.show()
plt.close()