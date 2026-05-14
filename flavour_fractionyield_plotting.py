import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style b-flavour fraction vs leading vs subleading jet pT
# Colour-blind safe (Okabe-Ito palette)
# =================================================

# Data is updated for the updated selection

b_yield_leading = np.array([
    [15, 866.832764, 7.127416],
    [45, 12639.586914, 27.295240],
    [75, 15869.648438, 30.743673],
    [105, 9926.534180, 24.359756],
    [135, 5014.506836, 17.360830],
    [165, 2387.997070, 12.003186],
    [195, 1107.404541, 8.184131],
    [225, 518.517212, 5.626997],
    [255, 246.061020, 3.902366],
    [285, 124.001411, 2.777023],
    [330, 97.957062, 2.472294],
    [630, 40.591053, 1.625272],
], dtype=float)

nonb_yield_leading = np.array([
    [15, 513.730942, 5.495756],
    [45, 2568.648971, 12.332749],
    [75, 1848.099151, 10.517399],
    [105, 1023.515342, 7.838409],
    [135, 570.466637, 5.854592],
    [165, 335.983661, 4.486757],
    [195, 208.996976, 3.542173],
    [225, 137.174612, 2.873576],
    [255, 95.567451, 2.379880],
    [285, 63.837834, 1.957455],
    [330, 78.560709, 2.171254],
    [630, 97.563583, 2.414945],
], dtype=float)

b_yield_subleading = np.array([
    [15, 741.942993, 6.583130],
    [45, 10381.327148, 24.727894],
    [75, 13325.693359, 28.174053],
    [105, 8366.264648, 22.363771],
    [135, 4313.876953, 16.103434],
    [165, 2148.345215, 11.378210],
    [195, 1047.976929, 7.953327],
    [225, 526.710571, 5.658462],
    [255, 277.989563, 4.116566],
    [285, 155.299957, 3.089904],
    [330, 147.985413, 3.012460],
    [630, 114.669403, 2.651663],
], dtype=float)

nonb_yield_subleading = np.array([
    [15, 638.803496, 6.138078],
    [45, 4826.655746, 16.900718],
    [75, 4391.090270, 16.184774],
    [105, 2582.961998, 12.435903],
    [135, 1270.946503, 8.737989],
    [165, 575.363068, 5.892906],
    [195, 268.310863, 4.032890],
    [225, 128.981239, 2.811056],
    [255, 63.638881, 1.986461],
    [285, 32.419154, 1.410197],
    [330, 28.532356, 1.323478],
    [630, 23.549892, 1.202490],
], dtype=float)

b_frac_leading = np.array([
    [15, 0.627883, 0.003152],
    [45, 0.831101, 0.000739],
    [75, 0.895692, 0.000562],
    [105, 0.906529, 0.000681],
    [135, 0.897857, 0.000993],
    [165, 0.876657, 0.001543],
    [195, 0.841236, 0.002469],
    [225, 0.790794, 0.003903],
    [255, 0.720259, 0.005949],
    [285, 0.660146, 0.008519],
    [330, 0.554942, 0.009244],
    [630, 0.293809, 0.009767],
], dtype=float)

b_frac_subleading = np.array([
    [15, 0.537349, 0.003251],
    [45, 0.682624, 0.000917],
    [75, 0.752151, 0.000792],
    [105, 0.764096, 0.000993],
    [135, 0.772429, 0.001375],
    [165, 0.788757, 0.001921],
    [195, 0.796161, 0.002733],
    [225, 0.803290, 0.003839],
    [255, 0.813719, 0.005237],
    [285, 0.827300, 0.006834],
    [330, 0.838360, 0.006864],
    [630, 0.829619, 0.007923],
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
    "leading": "#0072B2",
    "subleading": "#56B4E9",
    "nonb-leading": "#E69F00",
    "nonb-subleading": "#F0E442",
}

MARKERS = {
    "leading": "o",
    "subleading": "s",
    "nonb-leading": "D",
    "nonb-subleading": "v",
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
# Unpack yields
# -----------------------------
jet_leading, y_b_leading, ey_b_leading = unpack(b_yield_leading)
jet_subleading, y_b_subleading, ey_b_subleading = unpack(b_yield_subleading)

jet_nonb_leading, y_nonb_leading, ey_nonb_leading = unpack(nonb_yield_leading)
jet_nonb_subleading, y_nonb_subleading, ey_nonb_subleading = unpack(nonb_yield_subleading)

# -----------------------------
# Unpack b fractions
# -----------------------------
_, f_b_leading, ef_b_leading = unpack(b_frac_leading)
_, f_b_subleading, ef_b_subleading = unpack(b_frac_subleading)

# Consistency check
expected_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
if not np.allclose(jet_leading, expected_centres):
    raise ValueError("Bin centres in data do not match supplied bin edges.")

major_ticks = np.arange(0, 901, 90)
minor_ticks = np.array(
    [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360,
     450, 540, 630, 720, 810, 900],
    dtype=float
)

# =================================================
# Figure with b-fraction panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# -----------------------------
# Top panel: yields
# -----------------------------
plot_with_errorbars(
    ax, bin_edges, jet_leading, y_b_leading, ey_b_leading,
    "leading", "Leading [b]"
)

plot_with_errorbars(
    ax, bin_edges, jet_subleading, y_b_subleading, ey_b_subleading,
    "subleading", "Sub-leading [b]"
)

plot_with_errorbars(
    ax, bin_edges, jet_nonb_leading, y_nonb_leading, ey_nonb_leading,
    "nonb-leading", "Leading [non-b]"
)

plot_with_errorbars(
    ax, bin_edges, jet_nonb_subleading, y_nonb_subleading, ey_nonb_subleading,
    "nonb-subleading", "Sub-leading [non-b]"
)

ax.set_ylabel("Weighted Yield")
ax.set_title(r"b-Flavour Purity Comparison - $\chi^2$")

ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.set_ylim(0.0, 1.15 * max(
    np.max(y_b_leading),
    np.max(y_b_subleading),
    np.max(y_nonb_leading),
    np.max(y_nonb_subleading),
))

ax.set_xticks(major_ticks)
ax.set_xticklabels([])
ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax.set_xticks(minor_ticks, minor=True)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

ax.legend(
    loc="best",
    frameon=False,
    handlelength=1.8,
    borderaxespad=0.0,
)

# -----------------------------
# Bottom panel: b fractions
# -----------------------------
rax.stairs(
    f_b_leading,
    bin_edges,
    color=COLORS["leading"],
    linewidth=2.0,
    fill=False,
    label="Leading b fraction",
)

rax.stairs(
    f_b_subleading,
    bin_edges,
    color=COLORS["subleading"],
    linewidth=2.0,
    fill=False,
    label="Sub-leading b fraction",
)

rax.errorbar(
    jet_leading,
    f_b_leading,
    yerr=ef_b_leading,
    fmt=MARKERS["leading"],
    color=COLORS["leading"],
    markerfacecolor="white",
    markeredgecolor=COLORS["leading"],
    markersize=5,
    capsize=3,
    linestyle="none",
)

rax.errorbar(
    jet_subleading,
    f_b_subleading,
    yerr=ef_b_subleading,
    fmt=MARKERS["subleading"],
    color=COLORS["subleading"],
    markerfacecolor="white",
    markeredgecolor=COLORS["subleading"],
    markersize=5,
    capsize=3,
    linestyle="none",
)

rax.set_ylabel("b-Fraction")
rax.set_xlabel(r"Jet $p_T$ [GeV]")
rax.set_ylim(0.0, 1.1)

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{int(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.tight_layout()
plt.savefig("bAndNonbYield_BFraction_Comparison_Chi2.png")
plt.show()
plt.close()