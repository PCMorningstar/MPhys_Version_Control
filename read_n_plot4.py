import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs jet multiplicity
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [jet multiplicity, fraction, error]
# -----------------------------
cb_const_data = np.array([
    [2, 0.8228, 0.0005],
    [3, 0.4108, 0.0007],
    [4, 0.2555, 0.0008],
    [5, 0.1753, 0.0010],
    [6, 0.1313, 0.0015],
    [7, 0.1007, 0.0023],
    [8, 0.0789, 0.0037],
    [9, 0.0643, 0.0060],
    [10, 0.0435, 0.0098],
], dtype=float)

cb_data = np.array([
    [2, 0.8228, 0.0005],
    [3, 0.4105, 0.0007],
    [4, 0.2549, 0.0008],
    [5, 0.1742, 0.0010],
    [6, 0.1309, 0.0015],
    [7, 0.0994, 0.0023],
    [8, 0.0766, 0.0036],
    [9, 0.0638, 0.0060],
    [10, 0.0568, 0.0110],
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
    "Constant": "#0072B2",
    "Variable": "#56B4E9"
}

MARKERS = {
    "Constant": "o",
    "Variable": "s"
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
jet_bb, f_bb, e_bb = unpack(cb_const_data)
jet_bc, f_bc, e_bc = unpack(cb_data)


plot_edges = centres_to_edges(jet_bb, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, jet_bb, f_bb, e_bb, "Constant", "Constant")
plot_with_errorbars(ax, jet_bc, f_bc, e_bc, "Variable", "Variable")


ax.set_xlabel("Jet Multiplicity")
ax.set_ylabel("Fraction Of Events")
ax.set_title(r"$\chi^2$ Reconstruction vs. $\mu$ & $\sigma$ Parameters")

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
plt.savefig("chi2_reconstruction_vs_mu_sigma_parameters0.png")
plt.show()
plt.close()

# -----------------------------
# Stacked bar chart
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(jet_bb, dtype=float)

stack_order = [
    ("Constant", f_bb, "Constant"),
    ("Variable", f_bc, "Variable")
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
ax.set_title(r"$\chi^2$ Reconstruction vs $\mu$ & $\sigma$ Parameters")

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
plt.savefig("chi2_reconstruction_vs_mu_sigma_parameters.png")
plt.show()
plt.close()