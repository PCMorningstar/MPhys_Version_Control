import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Publication-style Chi2 pair flavour composition vs leading jet pT
# Colour-blind safe (Okabe-Ito palette)
# =================================================

b_data_65 = np.array([
    [15, 580109.500000, 761.6492],
    [45, 10300460.000000, 3209.4329],
    [75, 14077759.000000, 3752.0339],
    [105, 9112232.000000, 3018.6474],
    [135, 4731268.000000, 2175.1489],
    [165, 2328971.000000, 1526.0968],
    [195, 1107001.000000, 1052.1412],
    [225, 537910.000000, 733.4235],
    [255, 262333.375000, 512.1850],
    [285, 143838.875000, 379.2610],
    [330, 124835.984375, 353.3216],
    [630, 78921.617188, 280.9300],
], dtype=float)

b_data_70 = np.array([
    [15, 649230.312500, 805.7471],
    [45, 11174405.000000, 3342.8130],
    [75, 15098569.000000, 3885.6875],
    [105, 9733567.000000, 3121.4687],
    [135, 5043693.000000, 2245.8168],
    [165, 2475329.750000, 1573.3187],
    [195, 1181951.000000, 1087.1752],
    [225, 575047.937500, 758.3192],
    [255, 281609.093750, 530.6685],
    [285, 154094.046875, 392.5482],
    [330, 134148.656250, 366.2631],
    [630, 84404.296875, 290.5242],
], dtype=float)

b_data_77 = np.array([
    [15, 751593.812500, 866.9566],
    [45, 12453039.000000, 3528.8852],
    [75, 16509936.000000, 4063.2433],
    [105, 10565710.000000, 3250.4946],
    [135, 5462457.500000, 2337.1888],
    [165, 2675158.500000, 1635.5912],
    [195, 1277985.250000, 1130.4571],
    [225, 621761.875000, 788.5809],
    [255, 306527.187500, 553.6490],
    [285, 167401.750000, 409.1490],
    [330, 147669.718750, 384.2847],
    [630, 92918.250000, 304.8233],
], dtype=float)

b_data_85 = np.array([
    [15, 873348.875000, 934.5218],
    [45, 13929894.000000, 3733.6172],
    [75, 18106278.000000, 4255.1459],
    [105, 11496011.000000, 3390.5753],
    [135, 5926773.000000, 2434.2800],
    [165, 2891315.000000, 1700.3867],
    [195, 1385298.000000, 1177.0357],
    [225, 674416.437500, 821.2291],
    [255, 333556.187500, 577.5433],
    [285, 183233.656250, 428.0578],
    [330, 159065.312500, 398.8307],
    [630, 101171.671875, 318.0737],
], dtype=float)

b_data_90 = np.array([
    [15, 945597.312500, 972.4183],
    [45, 14824274.000000, 3850.2304],
    [75, 19127092.000000, 4373.4531],
    [105, 12093795.000000, 3477.6134],
    [135, 6211308.500000, 2492.2497],
    [165, 3034828.000000, 1742.6491],
    [195, 1442333.000000, 1200.9713],
    [225, 706846.062500, 840.7414],
    [255, 350844.156250, 592.3210],
    [285, 191232.218750, 437.3011],
    [330, 166992.156250, 408.6467],
    [630, 105747.328125, 325.1881],
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
jet_65, f_65, e_65 = unpack(b_data_65)
jet_70, f_70, e_70 = unpack(b_data_70)
jet_77, f_77, e_77 = unpack(b_data_77)
jet_85, f_85, e_85 = unpack(b_data_85)
jet_90, f_90, e_90 = unpack(b_data_90)

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
ax.set_ylabel("GN2 b-Tagged Events")
ax.set_title(r"GN2 MC b-Tagging Event Yield Comparison - $\chi^2$")

ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.set_ylim(0.0, 0.2e8)

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
plt.savefig("GN2_MC_bTaggingEventYield_Comparison_Chi2.png")
plt.show()
plt.close()