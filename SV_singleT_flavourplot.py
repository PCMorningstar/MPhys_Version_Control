import numpy as np
import matplotlib.pyplot as plt

################### chi2 combined - SV INVARIANT MASS ###################

# =================================================
# Combined b yields
# =================================================
b_yield = np.array([
    [-0.25, 18346.723633, 32.913587],
    [ 0.25,  2884.421997, 13.047619],
    [ 0.75, 13409.761230, 28.176313],
    [ 1.25, 20989.498047, 35.262175],
    [ 1.75, 27564.194336, 40.412556],
    [ 2.25, 30050.157226, 42.214343],
    [ 2.75, 28436.375000, 41.075650],
    [ 3.25, 22307.510742, 36.411997],
    [ 3.75, 14536.535156, 29.395959],
    [ 4.25,  8023.555176, 21.838867],
    [ 4.75,  3793.313476, 15.044420],
    [ 5.25,  1623.820190,  9.840460],
], dtype=float)

# =================================================
# Combined c yields
# =================================================
c_yield = np.array([
    [-0.25, 1380.454590, 9.010381],
    [ 0.25,  126.965889, 2.746112],
    [ 0.75,  501.097305, 5.437906],
    [ 1.25,  580.969238, 5.838259],
    [ 1.75,  568.441864, 5.791357],
    [ 2.25,  458.939728, 5.204337],
    [ 2.75,  372.710541, 4.680701],
    [ 3.25,  265.138993, 3.954210],
    [ 3.75,  167.208573, 3.132159],
    [ 4.25,   91.480133, 2.334776],
    [ 4.75,   49.303593, 1.717279],
    [ 5.25,   21.343060, 1.102722],
], dtype=float)

# =================================================
# Combined light yields
# =================================================
l_yield = np.array([
    [-0.25, 34356.345703, 44.996290],
    [ 0.25,  1996.076904, 10.856088],
    [ 0.75,  7235.364014, 20.679374],
    [ 1.25,  8726.052246, 22.716846],
    [ 1.75,  9297.580566, 23.432197],
    [ 2.25,  8548.253906, 22.490627],
    [ 2.75,  7115.196045, 20.515272],
    [ 3.25,  5126.545410, 17.436331],
    [ 3.75,  3181.104614, 13.732486],
    [ 4.25,  1724.376465, 10.111237],
    [ 4.75,   835.989746,  7.042694],
    [ 5.25,   397.331207,  4.855377],
], dtype=float)

# =================================================
# Combined b fractions
# =================================================
b_frac = np.array([
    [-0.25, 0.339229, 0.000495],
    [ 0.25, 0.576024, 0.001697],
    [ 0.75, 0.634145, 0.000805],
    [ 1.25, 0.692802, 0.000645],
    [ 1.75, 0.736416, 0.000554],
    [ 2.25, 0.769385, 0.000519],
    [ 2.75, 0.791564, 0.000521],
    [ 3.25, 0.805349, 0.000579],
    [ 3.75, 0.812785, 0.000710],
    [ 4.25, 0.815451, 0.000953],
    [ 4.75, 0.810778, 0.001396],
    [ 5.25, 0.795018, 0.002175],
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

BIN_WIDTH = 0.5

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
major_ticks = np.arange(-0.5, 5.5 + 0.5, 0.5)
minor_ticks = np.arange(-0.5, 5.5 + BIN_WIDTH, BIN_WIDTH)

# =================================================
# Figure with ratio panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: line histogram (yields)
# =================================================
plot_with_errorbars(ax, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, pt, c, c_err, "c", "c")
plot_with_errorbars(ax, pt, l, l_err, "l", "l")

ax.set_ylabel("Weighted Yield")
ax.set_title(r"$\chi^2$ Flavour Composition vs. SV Mass - Combined")
ax.set_xlim(plot_edges[0], plot_edges[-1])
ax.set_ylim(1e-2, 1e6)
ax.set_yscale("log", base=10)
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
rax.set_xlabel(r"SV Mass [GeV]")
rax.set_ylim(0.0, 1.1)
rax.set_xlim(plot_edges[0], plot_edges[-1])

# Major ticks every 90 GeV, minor ticks every 30 GeV
rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{float(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("chi2_combined_flavour_vs_SV_invariant_mass.png")
plt.show()
plt.close()