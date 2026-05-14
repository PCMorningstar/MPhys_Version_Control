import numpy as np
import matplotlib.pyplot as plt

################### hyper combined - Jet Multiplicity ###################

# =================================================
# b yields
# =================================================
b_yield = np.array([
    [2, 96844.870000, 75.899160],
    [3, 55081.480000, 57.149950],
    [4, 21646.160000, 35.775160],
    [5,  7204.078000, 20.624690],
    [6,  2170.109000, 11.324060],
    [7,   581.796900,  5.857521],
    [8,   151.765600,  3.004066],
    [9,    38.375000,  1.496089],
    [10,    9.218750,  0.724164],
], dtype=float)

# =================================================
# c yields
# =================================================
c_yield = np.array([
    [2, 1160.876000, 8.290800],
    [3,  959.765300, 7.532715],
    [4,  509.360800, 5.468034],
    [5,  207.069300, 3.476372],
    [6,   71.518070, 2.050971],
    [7,   19.728760, 1.065891],
    [8,    6.385742, 0.609438],
    [9,    1.868408, 0.328613],
    [10,   0.259521, 0.108886],
], dtype=float)

# =================================================
# light yields
# =================================================
l_yield = np.array([
    [2, 20775.900000, 35.121010],
    [3, 16694.460000, 31.405730],
    [4,  8333.195000, 22.162460],
    [5,  3229.832000, 13.767130],
    [6,  1062.867000,  7.887948],
    [7,   307.968800,  4.243101],
    [8,    81.843750,  2.189787],
    [9,    21.585940,  1.107172],
    [10,    5.570313,  0.559672],
], dtype=float)

# =================================================
# b fractions
# =================================================
b_frac = np.array([
    [2, 0.814944, 0.000429],
    [3, 0.757147, 0.000602],
    [4, 0.710163, 0.000930],
    [5, 0.676906, 0.001546],
    [6, 0.657513, 0.002699],
    [7, 0.638245, 0.004875],
    [8, 0.633481, 0.008931],
    [9, 0.619376, 0.016534],
    [10, 0.612058, 0.031867],
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
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: line histogram (yields)
# =================================================
plot_with_errorbars(ax, nj, b, b_err, "b", "b")
plot_with_errorbars(ax, nj, c, c_err, "c", "c")
plot_with_errorbars(ax, nj, l, l_err, "l", "l")

ax.set_ylabel("Weighted Yield")
ax.set_title(r"HyPER Flavour Composition vs. Jet Multiplicity - Combined")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(1e-2, 1e6)
ax.set_yscale("log", base=10)

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
rax.set_ylim(0, 1.1)

rax.set_xticks(nj)

rax.grid(True, linestyle=":", alpha=0.5)

# =================================================
# Final
# =================================================
#plt.tight_layout()

plt.savefig("hyper_combined_flavour_vs_jet_multiplicity.png")
plt.show()
plt.close()