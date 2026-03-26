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
    [75, 21110552.0, 43391.5703],
    [105, 13256744.0, 34446.0547],
    [135, 6775839.0, 24695.5215],
    [165, 3319361.0, 17307.2070],
    [195, 1578039.75, 11944.0654],
    [225, 769545.9375, 8374.5703],
    [255, 388249.8125, 5976.0986],
    [285, 203232.3438, 4337.1450],
    [315, 115425.8984, 3260.8491],
    [345, 64544.4492, 2451.5496],
    [375, 35956.3359, 1851.0824],
    [405, 24487.5820, 1505.5728],
    [435, 17759.7539, 1292.5221],
    [465, 10389.0371, 993.5731],
    [495, 6751.3574, 803.5013],
    [525, 5443.6250, 686.3103],
    [555, 2996.3906, 541.4762],
    [585, 3318.4653, 556.1020],
    [615, 1433.1758, 337.8076],
    [645, 1540.4990, 369.8357],
    [675, 1529.6481, 369.9860],
    [705, 1078.0552, 314.2697],
    [735, 1056.6492, 307.2703],
    [765, 810.5961, 269.8098],
    [795, 169.2023, 119.7465],
    [825, 0.0, 0.0],
    [855, 260.8063, 151.5731],
    [885, 0.0, 0.0],
], dtype=float)
c_yield = np.array([
    [15, 37934.4766, 1835.1229],
    [45, 250378.7813, 4702.5820],
    [75, 232733.1250, 4555.2705],
    [105, 145486.1250, 3606.9976],
    [135, 70002.1250, 2517.3989],
    [165, 36848.6289, 1831.9004],
    [195, 20601.8086, 1389.9183],
    [225, 10738.8203, 988.4081],
    [255, 7174.3652, 794.9268],
    [285, 5018.2183, 666.7678],
    [315, 2062.7432, 421.1903],
    [345, 2069.1150, 435.0470],
    [375, 1738.4304, 402.6272],
    [405, 1215.6693, 343.0108],
    [435, 827.3414, 267.2363],
    [465, 603.4117, 232.4057],
    [495, 938.9663, 300.3517],
    [525, 81.0169, 81.0169],
    [555, 377.6390, 191.3372],
    [585, 205.0949, 145.0254],
    [615, 161.7507, 114.3842],
    [645, 245.8441, 129.5389],
    [675, 159.1696, 112.6263],
    [705, 80.7937, 80.7937],
    [735, 0.0, 0.0],
    [765, 110.6735, 110.6735],
    [795, 0.0, 0.0],
    [825, 176.5603, 124.8895],
    [855, 85.0423, 85.0423],
    [885, 0.0, 0.0],
], dtype=float)
l_yield = np.array([
    [15, 804517.625, 8422.5352],
    [45, 5128971.5, 21327.8320],
    [75, 4314778.0, 19629.1016],
    [105, 2471826.75, 14891.5127],
    [135, 1278569.5, 10713.1875],
    [165, 622784.0, 7497.8398],
    [195, 327300.25, 5441.4087],
    [225, 182117.5938, 4076.4072],
    [255, 106109.0938, 3098.1531],
    [285, 68569.8438, 2488.8784],
    [315, 46340.9375, 2045.0286],
    [345, 29699.6953, 1638.8478],
    [375, 22139.0234, 1393.8093],
    [405, 12580.6543, 1062.5563],
    [435, 12444.9961, 1052.2178],
    [465, 9636.4043, 936.9900],
    [495, 6758.2070, 789.5706],
    [525, 4481.6655, 623.5134],
    [555, 3911.2793, 581.1094],
    [585, 2142.4856, 439.1186],
    [615, 2187.7134, 440.2750],
    [645, 1319.6367, 341.7672],
    [675, 1503.9375, 391.0549],
    [705, 533.5271, 220.1659],
    [735, 1117.3851, 313.6956],
    [765, 77.2006, 77.2006],
    [795, 679.4915, 249.7397],
    [825, 411.4770, 191.2211],
    [855, 240.7136, 139.2006],
    [885, 203.8207, 144.2401],
], dtype=float)

# -------------------------------------------------
# Weighted event fraction and error calculation
# -------------------------------------------------
b_frac = np.array([
    [15, 0.5818116665, 0.0032591184],
    [45, 0.7564768195, 0.0008593852],
    [75, 0.8227649331, 0.0007122969],
    [105, 0.8351202607, 0.0008819179],
    [135, 0.8340100050, 0.0012372678],
    [165, 0.8342213035, 0.0017715959],
    [195, 0.8193600178, 0.0026388801],
    [225, 0.7996094227, 0.0038969147],
    [255, 0.7741257548, 0.0056229340],
    [285, 0.7341668606, 0.0080028521],
    [315, 0.7045485973, 0.0107335216],
    [345, 0.6701512337, 0.0144805775],
    [375, 0.6009369493, 0.0190980011],
    [405, 0.6396312118, 0.0234275413],
    [435, 0.5723028183, 0.0267991917],
    [465, 0.5036167502, 0.0335714652],
    [495, 0.4672694802, 0.0402998016],
    [525, 0.5440193415, 0.0463319223],
    [555, 0.4112921357, 0.0557447907],
    [585, 0.5856757164, 0.0627579303],
    [615, 0.3788824081, 0.0717830708],
    [645, 0.4959784150, 0.0837145179],
    [675, 0.4790997207, 0.0858652740],
    [705, 0.6370068789, 0.1110672002],
    [735, 0.4860314727, 0.1009713036],
    [765, 0.8118380308, 0.1209255569],
    [795, 0.1993678808, 0.1272907450],
    [825, 0.0000000000, 0.0000000000],
    [855, 0.4446353316, 0.1894348081],
    [885, 0.0000000000, 0.0000000000],
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

BIN_WIDTH = 30.0

# =================================================
# Helpers
# =================================================
def centres_to_edges(x, width=BIN_WIDTH):
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
        linestyle="none",
        label=label,
        zorder=3,
    )

# =================================================
# Extract
# =================================================
pt = b_yield[:, 0]

b = b_yield[:, 1]
c = c_yield[:, 1]
l = l_yield[:, 1]

b_err = b_yield[:, 2]
c_err = c_yield[:, 2]
l_err = l_yield[:, 2]

frac_b = b_frac[:, 1]
err_b = b_frac[:, 2]

plot_edges = centres_to_edges(pt, BIN_WIDTH)
major_ticks = np.arange(0, 900 + 90, 90)
minor_ticks = np.arange(0, 900 + BIN_WIDTH, BIN_WIDTH)

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
plot_with_errorbars(ax, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, pt, c, c_err, "c", "c")
plot_with_errorbars(ax, pt, l, l_err, "l", "l")

ax.set_ylabel("Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet $p_T$ - Top2")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b fraction (black, step-style)
# =================================================
rax.step(
    plot_edges,
    np.r_[frac_b, frac_b[-1]],
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

rax.errorbar(
    pt,
    frac_b,
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
rax.set_xlim(plot_edges[0], plot_edges[-1])

# Major ticks every 90 GeV, minor ticks every 30 GeV
rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{int(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("chi2_top2_flavour_vs_jet_pT.png")
plt.show()
plt.close()