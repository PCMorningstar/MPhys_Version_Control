import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style b-flavour fraction vs leading vs subleading jet pT
# Colour-blind safe (Okabe-Ito palette)
# =================================================

b_yield_leading = np.array([
    [15, 1263964.875000, 10530.921875],
    [45, 18349602.000000, 40242.097656],
    [75, 22966610.000000, 45251.671875],
    [105, 14380384.000000, 35876.449219],
    [135, 7293083.000000, 25617.259766],
    [165, 3489088.750000, 17753.044922],
    [195, 1619549.500000, 12109.985352],
    [225, 762270.250000, 8348.328125],
    [255, 361576.500000, 5789.703613],
    [285, 183430.078125, 4129.299316],
    [330, 144699.593750, 3675.663818],
    [630, 59864.078125, 2414.115723],
], dtype=float)

nonb_yield_leading = np.array([
    [15, 750479.312500, 8127.854980],
    [45, 3739814.500000, 18207.851562],
    [75, 2690572.250000, 15527.547852],
    [105, 1494440.125000, 11586.583008],
    [135, 831508.187500, 8645.570312],
    [165, 490262.250000, 6631.447266],
    [195, 306391.437500, 5248.756836],
    [225, 200132.562500, 4246.512207],
    [255, 139956.968750, 3524.711670],
    [285, 93569.984375, 2900.909424],
    [330, 115443.132812, 3219.837891],
    [630, 144494.125000, 3596.121582],
], dtype=float)

b_yield_subleading = np.array([
    [15, 1084860.625000, 9740.370117],
    [45, 15072922.000000, 36458.628906],
    [75, 19289746.000000, 41473.875000],
    [105, 12120180.000000, 32935.507812],
    [135, 6268660.000000, 23750.597656],
    [165, 3135640.500000, 16818.328125],
    [195, 1531205.625000, 11763.998047],
    [225, 772552.437500, 8384.850586],
    [255, 407981.625000, 6103.521484],
    [285, 228780.359375, 4588.052734],
    [330, 218103.578125, 4474.310059],
    [630, 169482.359375, 3944.111572],
], dtype=float)

nonb_yield_subleading = np.array([
    [15, 929860.062500, 9061.609375],
    [45, 7016444.500000, 24934.074219],
    [75, 6367228.000000, 23846.957031],
    [105, 3753354.500000, 18344.455078],
    [135, 1855722.500000, 12918.352539],
    [165, 843311.062500, 8732.306641],
    [195, 394564.312500, 5982.907227],
    [225, 189851.109375, 4173.944824],
    [255, 93551.656250, 2948.103516],
    [285, 48040.164062, 2097.543457],
    [330, 42039.285156, 1964.285400],
    [630, 34972.796875, 1792.599243],
], dtype=float)

b_frac_leading = np.array([
    [15, 0.627449, 0.003194],
    [45, 0.830702, 0.000751],
    [75, 0.895152, 0.000572],
    [105, 0.905861, 0.000695],
    [135, 0.897656, 0.001008],
    [165, 0.876797, 0.001561],
    [195, 0.840915, 0.002501],
    [225, 0.792048, 0.003933],
    [255, 0.720942, 0.006004],
    [285, 0.662202, 0.008570],
    [330, 0.556231, 0.009312],
    [630, 0.292937, 0.009815],
], dtype=float)

b_frac_subleading = np.array([
    [15, 0.538466, 0.003293],
    [45, 0.682372, 0.000932],
    [75, 0.751876, 0.000806],
    [105, 0.763543, 0.001010],
    [135, 0.771588, 0.001397],
    [165, 0.788057, 0.001948],
    [195, 0.795114, 0.002769],
    [225, 0.802732, 0.003883],
    [255, 0.813468, 0.005293],
    [285, 0.826458, 0.006891],
    [330, 0.838399, 0.006914],
    [630, 0.828946, 0.007982],
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

ax.set_ylabel("Events")
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
rax.set_ylim(0.0, 1.0)

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{int(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.tight_layout()
plt.savefig("bAndNonbYield_BFraction_Comparison_Chi2.png")
plt.show()
plt.close()