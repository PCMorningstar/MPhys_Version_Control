import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style full reconstruction efficiency
# vs jet multiplicity for different truth-level fits
# =================================================

# -----------------------------
# Hard-coded data: Jet Multiplicity, efficiency, error
# -----------------------------

chi2_constant_sigma_mu = np.array([
    [2, 0.8209, 0.0005],
    [3, 0.4087, 0.0007],
    [4, 0.2537, 0.0008],
    [5, 0.1740, 0.0010],
    [6, 0.1307, 0.0015],
    [7, 0.0995, 0.0023],
    [8, 0.0771, 0.0036],
    [9, 0.0628, 0.0059],
    [10, 0.0458, 0.0100],
], dtype=float)

chi2_non_constant_sigma_mu = np.array([
    [2, 0.8209, 0.0005],
    [3, 0.4083, 0.0007],
    [4, 0.2531, 0.0008],
    [5, 0.1734, 0.0010],
    [6, 0.1299, 0.0015],
    [7, 0.0990, 0.0023],
    [8, 0.0762, 0.0036],
    [9, 0.0624, 0.0059],
    [10, 0.0561, 0.0109],
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
    r"Non-Constant": "#56B4E9",
}

MARKERS = {
    r"Constant": "o",
    r"Non-Constant": "s",
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

jet_constant_sigma_mu, eff_constant_sigma_mu, err_constant_sigma_mu = unpack(chi2_constant_sigma_mu)
jet_non_constant_sigma_mu, eff_non_constant_sigma_mu, err_non_constant_sigma_mu = unpack(chi2_non_constant_sigma_mu)

edges = centres_to_edges(jet_constant_sigma_mu, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, jet_constant_sigma_mu, eff_constant_sigma_mu, err_constant_sigma_mu, r"Constant", r"Constant")
plot_with_errorbars(ax, jet_non_constant_sigma_mu, eff_non_constant_sigma_mu, err_non_constant_sigma_mu, r"Non-Constant", r"Non-Constant")


ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Full Reconstruction Efficiency")
ax.set_title(r"$\chi^2$ Full Reconstruction Efficiency vs. Fit Parameter Type")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(0.0, 1.05)

ax.set_xticks(jet_constant_sigma_mu)
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
plt.savefig("corrected_full_reconstruction_efficiency_vs_fit_parameter_type.png")
plt.show()
plt.close()