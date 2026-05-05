import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style full reconstruction efficiency
# vs jet multiplicity for different truth-level fits
# =================================================

# -----------------------------
# Hard-coded data: Jet Multiplicity, efficiency, error
# -----------------------------

chi_crystal_ball = np.array([ # varying sigma and mu
    [2, 0.8201, 0.0005],
    [3, 0.4083, 0.0007],
    [4, 0.2533, 0.0008],
    [5, 0.1732, 0.0010],
    [6, 0.1304, 0.0015],
    [7, 0.0988, 0.0023],
    [8, 0.0760, 0.0036],
    [9, 0.0628, 0.0059],
    [10, 0.0464, 0.0100],
], dtype=float)
chi_const_sigma_mu = np.array([
    [2, 0.8201, 0.0005],
    [3, 0.4087, 0.0007],
    [4, 0.2538, 0.0008],
    [5, 0.1743, 0.0010],
    [6, 0.1306, 0.0015],
    [7, 0.1003, 0.0023],
    [8, 0.0782, 0.0036],
    [9, 0.0628, 0.0059],
    [10, 0.0430, 0.0097],
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
    r"Constant": "#0072B2",
    r"Non-Constant": "#56B4E9"
}

MARKERS = {
    r"Constant": "o",
    r"Non-Constant": "s"
}

BIN_WIDTH = 1.0  # jet multiplicity spacing is 1

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
        linewidth=0,
        label=label,
        zorder=3,
    )

# -----------------------------
# Unpack
# -----------------------------

jet_chi_crystal_ball, eff_chi_crystal_ball, err_chi_crystal_ball = unpack(chi_crystal_ball)
jet_chi_const_sigma_mu, eff_chi_const_sigma_mu, err_chi_const_sigma_mu = unpack(chi_const_sigma_mu)

edges = centres_to_edges(jet_chi_const_sigma_mu, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, jet_chi_crystal_ball, eff_chi_crystal_ball, err_chi_crystal_ball, r"Non-Constant", r"Non-Constant")
plot_with_errorbars(ax, jet_chi_const_sigma_mu, eff_chi_const_sigma_mu, err_chi_const_sigma_mu, r"Constant", r"Constant")


ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Full Reconstruction Efficiency")
ax.set_title(r"$\chi^2$ Reconstruction Efficiency vs. Fit Parameter Type")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.05)

ax.set_xticks(jet_chi_const_sigma_mu)
ax.minorticks_on()

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="upper right",
    #bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    #borderaxespad=0.0,
    #handlelength=1.8,
)

plt.tight_layout()
plt.savefig("chi2_reconstruction_efficiency_vs_fit_parameter_type.png")
plt.show()
plt.close()