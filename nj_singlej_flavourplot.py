import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 flavour composition vs jet pT
# Top2, 3-jet region
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# -----------------------------
# Hard-coded data: [jet pT bin centre, fraction, error]
# -----------------------------
b_data = np.array([
    [15, 0.427467, 0.002773],
    [45, 0.585032, 0.000987],
    [75, 0.676344, 0.001000],
    [105, 0.716713, 0.001294],
    [135, 0.734071, 0.001774],
    [165, 0.748685, 0.002450],
    [195, 0.754190, 0.003401],
    [225, 0.762043, 0.004636],
    [255, 0.774155, 0.006311],
    [285, 0.766883, 0.008826],
], dtype=float)

c_data = np.array([
    [15, 0.023944, 0.000851],
    [45, 0.018299, 0.000269],
    [75, 0.015651, 0.000265],
    [105, 0.013792, 0.000335],
    [135, 0.012896, 0.000448],
    [165, 0.013217, 0.000643],
    [195, 0.011297, 0.000838],
    [225, 0.012089, 0.001176],
    [255, 0.011634, 0.001614],
    [285, 0.010458, 0.002016],
], dtype=float)

l_data = np.array([
    [15, 0.548589, 0.002789],
    [45, 0.396670, 0.000980],
    [75, 0.308006, 0.000986],
    [105, 0.269495, 0.001274],
    [135, 0.253033, 0.001747],
    [165, 0.238098, 0.002406],
    [195, 0.234512, 0.003347],
    [225, 0.225868, 0.004555],
    [255, 0.214211, 0.006191],
    [285, 0.222659, 0.008698],
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
    "b": "#0072B2",
    "c": "#56B4E9",
    "l": "#E69F00",
}

MARKERS = {
    "b": "o",
    "c": "s",
    "l": "^",
}

BIN_WIDTH = 30.0  # GeV

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def centres_to_edges(x, width):
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
pt_b, f_b, e_b = unpack(b_data)
pt_c, f_c, e_c = unpack(c_data)
pt_l, f_l, e_l = unpack(l_data)

plot_edges = centres_to_edges(pt_b, BIN_WIDTH)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, pt_b, f_b, e_b, "b", "b")
plot_with_errorbars(ax, pt_c, f_c, e_c, "c", "c")
plot_with_errorbars(ax, pt_l, f_l, e_l, "l", "l")

ax.set_xlabel(r"Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction of events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet $p_T$ [Top2, 3-Jets]")

ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.05)

# Compute bin edges
edges = centres_to_edges(pt_b, BIN_WIDTH)

# Set ticks at edges
ax.set_xticks(edges)

# Label edges (clean integers)
ax.set_xticklabels([str(int(x)) for x in edges])
ax.minorticks_on()

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    borderaxespad=0.0,
    handlelength=1.8,
)

plt.tight_layout()
plt.savefig("chi2_top2_flavour_vs_jet_pT_3jets.png")
plt.show()
plt.close()

# -----------------------------
# Stacked bar chart
# -----------------------------
fig, ax = plt.subplots()

bottom = np.zeros_like(pt_b, dtype=float)

stack_order = [
    ("b", f_b, "b"),
    ("c", f_c, "c"),
    ("l", f_l, "l"),
]

for key, vals, label in stack_order:
    ax.bar(
        pt_b,
        vals,
        bottom=bottom,
        width=BIN_WIDTH,
        align="center",
        color=COLORS[key],
        edgecolor="white",
        linewidth=0.7,
        label=label,
    )
    bottom += vals

ax.set_xlabel(r"Jet $p_T$ [GeV]")
ax.set_ylabel("Fraction of events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet $p_T$ [Top2, 3-Jets]")

ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(0.0, 1.0)

# Compute bin edges - DONT FORGET
edges = centres_to_edges(pt_b, BIN_WIDTH)

# Set ticks at edges
ax.set_xticks(edges)

# Label edges (clean integers)
ax.set_xticklabels([str(int(x)) for x in edges])
ax.minorticks_on()

ax.grid(True, axis="y", which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, axis="y", which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    borderaxespad=0.0,
)

plt.tight_layout()
plt.savefig("chi2_top2_flavour_stacked_vs_jet_pT_3jets.png")
plt.show()
plt.close()