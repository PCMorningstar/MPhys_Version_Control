import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style HyPER pair flavour composition vs leading SV mass
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: Jet Multiplicity, efficiency, error
# -----------------------------
chi_gauss = np.array([
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

chi_crystal_ball = np.array([
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

chi_double_gaussian = np.array([
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
    "Guassian": "#0072B2",
    "Crystal Ball": "#56B4E9",
    "Double Gaussian": "#E69F00"
}

MARKERS = {
    "Guassian": "o",
    "Crystal Ball": "s",
    "Double Gaussian": "^"
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
func_gauss, f_gauss, e_gauss = unpack(chi_gauss)
func_crystal_ball, f_crystal_ball, e_crystal_ball = unpack(chi_crystal_ball)
func_double_gaussian, f_double_gaussian, e_double_gaussian = unpack(chi_double_gaussian)

edges = centres_to_edges(func_gauss, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, func_gauss, f_gauss, e_gauss, "Guassian", "Guassian")
plot_with_errorbars(ax, func_crystal_ball, f_crystal_ball, e_crystal_ball, "Crystal Ball", "Crystal Ball")
plot_with_errorbars(ax, func_double_gaussian, f_double_gaussian, e_double_gaussian, "Double Gaussian", "Double Gaussian")

ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Full Reconstruction Efficiency")
ax.set_title(r"$\chi^2$ vs Jet Multiplicity & Truth-Level Fitting Functions")

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
plt.savefig("chi2_vs_jet_multiplicity_and_truth_level_fitting_functions.png")
plt.show()
plt.close()

