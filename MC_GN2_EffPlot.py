import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs leading jet pT
# Colour-blind safe (Okabe-Ito palette)
# =================================================

b_eff_65 = np.array([
    [15,  0.524621, 0.004161],
    [45,  0.649242, 0.001047],
    [75,  0.698219, 0.000904],
    [105, 0.720919, 0.001119],
    [135, 0.733801, 0.001553],
    [165, 0.736045, 0.002245],
    [195, 0.731648, 0.003311],
    [225, 0.730936, 0.004848],
    [255, 0.719486, 0.007189],
    [285, 0.696439, 0.010360],
    [330, 0.680840, 0.011831],
    [630, 0.667471, 0.018866],
], dtype=float)

b_eff_70 = np.array([
    [15,  0.579578, 0.004113],
    [45,  0.698096, 0.001007],
    [75,  0.742651, 0.000861],
    [105, 0.764602, 0.001059],
    [135, 0.776069, 0.001465],
    [165, 0.777403, 0.002117],
    [195, 0.777896, 0.003109],
    [225, 0.774775, 0.004566],
    [255, 0.768368, 0.006736],
    [285, 0.754475, 0.009710],
    [330, 0.750404, 0.010982],
    [630, 0.723922, 0.017982],
], dtype=float)

b_eff_77 = np.array([
    [15,  0.661183, 0.003942],
    [45,  0.767356, 0.000927],
    [75,  0.803834, 0.000782],
    [105, 0.820794, 0.000957],
    [135, 0.830777, 0.001318],
    [165, 0.832765, 0.001900],
    [195, 0.833886, 0.002787],
    [225, 0.827845, 0.004123],
    [255, 0.823454, 0.006085],
    [285, 0.818409, 0.008709],
    [330, 0.828173, 0.009581],
    [630, 0.803353, 0.015977],
], dtype=float)

b_eff_85 = np.array([
    [15,  0.761995, 0.003552],
    [45,  0.846353, 0.000791],
    [75,  0.869811, 0.000663],
    [105, 0.882621, 0.000804],
    [135, 0.889859, 0.001100],
    [165, 0.889439, 0.001597],
    [195, 0.894255, 0.002303],
    [225, 0.887510, 0.003447],
    [255, 0.883310, 0.005136],
    [285, 0.887889, 0.007214],
    [330, 0.889170, 0.007890],
    [630, 0.882595, 0.013019],
], dtype=float)

b_eff_90 = np.array([
    [15,  0.819365, 0.003207],
    [45,  0.896719, 0.000668],
    [75,  0.913502, 0.000554],
    [105, 0.921746, 0.000671],
    [135, 0.926272, 0.000918],
    [165, 0.926943, 0.001326],
    [195, 0.927833, 0.001929],
    [225, 0.925116, 0.002874],
    [255, 0.931102, 0.004069],
    [285, 0.921725, 0.006163],
    [330, 0.928438, 0.006504],
    [630, 0.915077, 0.011389],
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
    "65%": "#0072B2",
    "70%": "#56B4E9",
    "77%": "#E69F00",
    "85%": "#CC79A7",
    "90%": "#F0E442",
}

MARKERS = {
    "65%": "o",
    "70%": "s",
    "77%": "^",
    "85%": "D",
    "90%": "v",
}

# =================================================
# Actual binning
# =================================================
# Bins are:
# 0-30, 30-60, ..., 270-300, 300-360, 360-900
bin_edges = np.array(
    [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 900],
    dtype=float
)

# -----------------------------
# Helpers
# -----------------------------
def unpack(arr):
    return arr[:, 0], arr[:, 1], arr[:, 2]

def plot_with_errorbars(ax, edges, x, y, yerr, key, label):
    ax.stairs(
        y,
        edges,
        color=COLORS[key],
        linewidth=2.0,
        fill=False,
        zorder=2,
        label=label,
    )

    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt=MARKERS[key],
        color=COLORS[key],
        markerfacecolor="white",
        markeredgecolor=COLORS[key],
        markersize=6,
        markeredgewidth=1.2,
        ecolor=COLORS[key],
        elinewidth=1.2,
        capsize=3,
        linestyle="none",
        zorder=3,
    )

# -----------------------------
# Unpack
# -----------------------------
jet_65, f_65, e_65 = unpack(b_eff_65)
jet_70, f_70, e_70 = unpack(b_eff_70)
jet_77, f_77, e_77 = unpack(b_eff_77)   
jet_85, f_85, e_85 = unpack(b_eff_85)
jet_90, f_90, e_90 = unpack(b_eff_90)

# Consistency check
expected_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
if not np.allclose(jet_65, expected_centres):
    raise ValueError("Bin centres in data do not match supplied bin edges.")

major_ticks = np.arange(0, 901, 90)
minor_ticks = np.array(
    [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 450, 540, 630, 720, 810, 900],
    dtype=float
)

# -----------------------------
# Main plot
# -----------------------------
fig, ax = plt.subplots()

plot_with_errorbars(ax, bin_edges, jet_65, f_65, e_65, "65%", "WP65")
plot_with_errorbars(ax, bin_edges, jet_70, f_70, e_70, "70%", "WP70")
plot_with_errorbars(ax, bin_edges, jet_77, f_77, e_77, "77%", "WP77")
plot_with_errorbars(ax, bin_edges, jet_85, f_85, e_85, "85%", "WP85")
plot_with_errorbars(ax, bin_edges, jet_90, f_90, e_90, "90%", "WP90")

ax.set_xlabel(r"Jet $p_T$ [GeV]")
ax.set_ylabel("b-Tagging Efficiency")
ax.set_title(r"GN2 MC Efficiencies - $\chi^2$ - Leading Jet")

ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.set_ylim(0.4, 1.0)

ax.set_xticks(major_ticks)
ax.set_xticklabels([f"{int(x)}" for x in major_ticks])
ax.set_xticks(minor_ticks, minor=True)

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
plt.savefig("GN2_MC_bTaggingEfficiency_Comparison_Chi2_leading_jet.png")
plt.show()
plt.close()