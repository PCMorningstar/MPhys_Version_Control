import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Cutflow yields (hard-coded)
# Format: [step_index, yield, error]
# =================================================
cutflow_data = np.array([
    [0, 500113664.0, 211947.328125],  # Raw
    [1, 336482720.0, 173903.171875],  # Electron
    [2, 229222096.0, 142321.890625],  # Muon
    [3, 228803616.0, 142192.812500],  # Jet
    [4, 208120864.0, 135615.656250],  # Dilepton
], dtype=float)

labels = ["Raw", "Electron", "Muon", "Jet", "Dilepton"]

# =================================================
# Style (same structure as your plot code)
# =================================================
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

COLOR = "#0072B2"  # Okabe–Ito blue
MARKER = "o"

BIN_WIDTH = 1.0

# =================================================
# Helpers (same logic as your histogram code)
# =================================================
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def centres_to_edges(x, width=1.0):
    x = np.asarray(x, dtype=float)
    return np.concatenate(([x[0] - width/2], x + width/2))

def plot_with_errorbars(ax, x, y, yerr, label, bin_width=1.0):
    edges = centres_to_edges(x, width=bin_width)

    # Step histogram
    ax.step(
        edges,
        np.r_[y, y[-1]],
        where="post",
        color=COLOR,
        linewidth=2.0,
        zorder=2
    )

    # Marker + uncertainty
    ax.errorbar(
        x, y, yerr=yerr,
        fmt=MARKER,
        color=COLOR,
        markerfacecolor="white",
        markeredgecolor=COLOR,
        markersize=6,
        markeredgewidth=1.2,
        ecolor=COLOR,
        elinewidth=1.2,
        capsize=3,
        linewidth=0,
        label=label,
        zorder=3
    )

# =================================================
# Unpack
# =================================================
x, y, yerr = unpack(cutflow_data)
edges = centres_to_edges(x, width=BIN_WIDTH)

# =================================================
# Plot
# =================================================
fig, ax = plt.subplots()

plot_with_errorbars(ax, x, y, yerr, "Cutflow Yield", bin_width=BIN_WIDTH)

ax.set_xlabel("Selection Section Step")
ax.set_ylabel("Weighted Events")
ax.set_title("Selection Cutflow - Specific Sections")

ax.set_xlim(edges[0], edges[-1])
ax.set_xticks(x)
ax.set_xticklabels(labels)

ax.minorticks_on()

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(frameon=False)

plt.tight_layout()
plt.savefig("cutflow_sections_linehist.png")
plt.show()
plt.close()