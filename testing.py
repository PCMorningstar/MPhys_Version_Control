import numpy as np
import matplotlib.pyplot as plt

# =================================================
# T2 flavour analysis - Jet pT region
# =================================================
# Weighted yield and error calculation
# -------------------------------------------------
b_yield = np.array([
    [15, 1172076.0, 10132.2334],
    [45, 16710323.0, 38392.6445],
    [75, 21110552.0, 43391.5742],
    [105, 13256742.0, 34446.0547],
    [135, 6775838.0, 24695.5215],
    [165, 3319361.0, 17307.2070],
    [195, 1578039.75, 11944.0654],
    [225, 769545.9375, 8374.5693],
    [255, 388249.8125, 5976.0981],
    [285, 203232.3281, 4337.1450],
    [330, 179970.3438, 4079.6118],
    [630, 114981.1797, 3280.4287],
], dtype=float)

c_yield = np.array([
    [15, 37934.4805, 1835.1229],
    [45, 250378.7813, 4702.5820],
    [75, 232733.1250, 4555.2705],
    [105, 145486.1250, 3606.9976],
    [135, 70002.1328, 2517.3989],
    [165, 36848.6250, 1831.9004],
    [195, 20601.8105, 1389.9183],
    [225, 10738.8203, 988.4081],
    [255, 7174.3657, 794.9268],
    [285, 5018.2183, 666.7678],
    [330, 4131.8579, 605.5305],
    [630, 7007.4043, 802.3493],
], dtype=float)

l_yield = np.array([
    [15, 804517.6250, 8422.5352],
    [45, 5128971.0, 21327.8340],
    [75, 4314777.5, 19629.1016],
    [105, 2471826.5, 14891.5127],
    [135, 1278569.5, 10713.1875],
    [165, 622784.0, 7497.8398],
    [195, 327300.25, 5441.4087],
    [225, 182117.5938, 4076.4072],
    [255, 106109.0938, 3098.1531],
    [285, 68569.8438, 2488.8782],
    [330, 76040.6250, 2620.6802],
    [630, 82369.6250, 2712.0237],
], dtype=float)

# -------------------------------------------------
# Weighted event fraction and error calculation
# -------------------------------------------------
b_frac = np.array([
    [15, 0.5818117261, 0.0032591186],
    [45, 0.7564767599, 0.0008593851],
    [75, 0.8227648139, 0.0007122966],
    [105, 0.8351202011, 0.0008819182],
    [135, 0.8340098262, 0.0012372674],
    [165, 0.8342212439, 0.0017715958],
    [195, 0.8193599582, 0.0026388799],
    [225, 0.7996094227, 0.0038969154],
    [255, 0.7741257548, 0.0056229348],
    [285, 0.7341668010, 0.0080028529],
    [330, 0.6918135285, 0.0086326723],
    [630, 0.5626452565, 0.0104843600],
], dtype=float)

# =================================================
# Style
# =================================================
plt.rcParams.update({
    "figure.figsize": (9.0, 7.0),
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

# =================================================
# Actual binning
# =================================================
# Bins are:
# 0-30, 30-60, ..., 270-300, 300-360, 360-900
bin_edges = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 900], dtype=float)

# Optional consistency check
expected_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
pt = b_yield[:, 0]
if not np.allclose(pt, expected_centres):
    raise ValueError("Bin centres in data do not match the supplied bin edges.")

# =================================================
# Helpers
# =================================================
def plot_with_errorbars(ax, edges, x, y, yerr, key, label):
    # Proper histogram drawing (no vertical seams or spikes)
    ax.stairs(
        y,
        edges,
        color=COLORS[key],
        linewidth=2.0,
        fill=False,
        label=label,
        zorder=2,
    )

    # error bars
    xerr = np.vstack((x - edges[:-1], edges[1:] - x))

    ax.errorbar(
        x, y,
        xerr=xerr,
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

# =================================================
# Extract
# =================================================
b = b_yield[:, 1]
c = c_yield[:, 1]
l = l_yield[:, 1]

b_err = b_yield[:, 2]
c_err = c_yield[:, 2]
l_err = l_yield[:, 2]

frac_b = b_frac[:, 1]
err_b = b_frac[:, 2]

major_ticks = np.arange(0, 901, 90)
minor_ticks = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 450, 540, 630, 720, 810, 900])

# =================================================
# Figure with ratio panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: line histogram (yields)
# =================================================
plot_with_errorbars(ax, bin_edges, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, bin_edges, pt, c, c_err, "c", "c")
plot_with_errorbars(ax, bin_edges, pt, l, l_err, "l", "l")

ax.set_ylabel("Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet $p_T$ - Top2")
ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b fraction
# =================================================
rax.step(
    bin_edges,
    np.r_[frac_b, frac_b[-1]],
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

xerr_frac = np.vstack((pt - bin_edges[:-1], bin_edges[1:] - pt))

rax.errorbar(
    pt,
    frac_b,
    xerr=xerr_frac,
    yerr=err_b,
    fmt="o",
    color="black",
    markerfacecolor="white",
    markeredgecolor="black",
    markersize=4,
    markeredgewidth=1.1,
    ecolor="black",
    elinewidth=1.2,
    capsize=3,
    linestyle="none",
    zorder=3,
)

rax.set_ylabel("b Fraction")
rax.set_xlabel(r"Jet $p_T$ [GeV]")
rax.set_ylim(0.0, 1.0)
rax.set_xlim(bin_edges[0], bin_edges[-1])

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{int(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("chi2_top2_flavour_vs_jet_pT.png")
plt.show()
plt.close()