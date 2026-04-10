import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Data (yields + b fraction)
# =================================================
b_yield = np.array([
    [2, 65579652.0, 76486.640625],
    [3, 43396952.0, 62030.56640625],
    [4, 19657534.0, 41645.37890625],
    [5, 7366859.5, 25449.76171875],
    [6, 2474501.75, 14714.4072265625],
    [7, 734081.375, 7991.5068359375],
    [8, 202954.375, 4186.7260742188],
    [9, 57400.90625, 2247.84375],
    [10, 15510.7099609375, 1153.5701904297],
], dtype=float)
c_yield = np.array([
    [2, 828055.75, 8598.359375],
    [3, 1065230.125, 9703.0703125],
    [4, 755854.5625, 8136.8315429688],
    [5, 399746.6875, 5912.3569335938],
    [6, 174852.265625, 3904.9409179688],
    [7, 60902.375, 2303.5249023438],
    [8, 22543.9921875, 1391.8453369141],
    [9, 6581.0786132812, 760.9968261719],
    [10, 1713.2377929688, 379.3201293945],
], dtype=float)
l_yield = np.array([
    [2, 15464859.0, 37150.8046875],
    [3, 21936354.0, 44075.0703125],
    [4, 14961321.0, 36311.765625],
    [5, 7446050.5, 25547.6953125],
    [6, 3031623.0, 16299.40234375],
    [7, 1094413.75, 9737.62109375],
    [8, 360017.53125, 5618.1669921875],
    [9, 107794.8125, 3048.0871582031],
    [10, 29751.240234375, 1583.9177246094],
], dtype=float)

# -------------------------------------------------
# Weighted event fraction and error calculation
# -------------------------------------------------
b_frac = np.array([
    [2, 0.8009967208, 0.0004168279],
    [3, 0.6535828710, 0.0005496177],
    [4, 0.5556945205, 0.0007844154],
    [5, 0.4842585623, 0.0012005092],
    [6, 0.4355767965, 0.0019464498],
    [7, 0.3885267377, 0.0033050032],
    [8, 0.3466248512, 0.0057937907],
    [9, 0.3341598213, 0.0106427476],
    [10, 0.3301894069, 0.0200403550],
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
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
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

BIN_WIDTH = 1.0

# =================================================
# Helpers
# =================================================
def centres_to_edges(x, width):
    x = np.asarray(x, dtype=float)
    return np.concatenate(([x[0] - width/2], x + width/2))

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
        markerfacecolor="white",
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

# =================================================
# Extract
# =================================================
nj = b_yield[:, 0]

b = b_yield[:, 1]
c = c_yield[:, 1]
l = l_yield[:, 1]

b_err = b_yield[:, 2]
c_err = c_yield[:, 2]
l_err = l_yield[:, 2]

frac_b = b_frac[:, 1]
err_b = b_frac[:, 2]

edges = centres_to_edges(nj, BIN_WIDTH)

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
plot_with_errorbars(ax, nj, b, b_err, "b", "b")
plot_with_errorbars(ax, nj, c, c_err, "c", "c")
plot_with_errorbars(ax, nj, l, l_err, "l", "l")

ax.set_ylabel("Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet Multiplicity - Top2")

ax.set_xlim(edges[0], edges[-1])

# REMOVE x labels on top plot
ax.tick_params(labelbottom=False)

ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b fraction (black)
# =================================================
# Step line (histogram style)
rax.step(
    edges,
    np.r_[frac_b, frac_b[-1]],
    where="post",
    color="black",
    linewidth=2.0,
)

# Error bars (points only, no connecting line)
rax.errorbar(
    nj,
    frac_b,
    yerr=err_b,
    fmt='o',
    color='black',
    markersize=4,
    capsize=3,
    linewidth=0,
)

rax.set_ylabel("b Fraction")
rax.set_xlabel("Jet Multiplicity")
rax.set_ylim(0, 1)

rax.set_xticks(nj)

rax.grid(True, linestyle=":", alpha=0.5)

# =================================================
# Final
# =================================================
#plt.tight_layout()

plt.savefig("chi2_top2_flavour_vs_jet_multiplicity.png")
plt.show()
plt.close()