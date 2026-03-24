import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Data (yields + b fraction)
# =================================================
b_yield = np.array([
    [2, 65597776.0, 76494.34375],
    [3, 43462832.0, 62070.24609375],
    [4, 19623906.0, 41603.78515625],
    [5, 7380257.0, 25465.98828125],
    [6, 2477211.0, 14760.00390625],
    [7, 748341.125, 8064.1650390625],
    [8, 215984.515625, 4351.6538085938],
    [9, 59451.7734375, 2265.2238769531],
    [10, 16522.642578125, 1171.0747070312],
], dtype=float)

c_yield = np.array([
    [2, 829367.375, 8586.912109375],
    [3, 1079626.75, 9766.158203125],
    [4, 768266.375, 8216.3125],
    [5, 407414.28125, 5966.9692382812],
    [6, 175589.34375, 3926.6940917969],
    [7, 60320.8203125, 2274.6611328125],
    [8, 19154.72265625, 1277.7174072266],
    [9, 7318.7626953125, 799.4478759766],
    [10, 2290.0251464844, 429.6172790527],
], dtype=float)

l_yield = np.array([
    [2, 15445648.0, 37137.96484375],
    [3, 21855014.0, 44004.28125],
    [4, 14981479.0, 36340.31640625],
    [5, 7424835.5, 25518.6015625],
    [6, 3028300.25, 16253.1171875],
    [7, 1080631.375, 9683.814453125],
    [8, 350376.6875, 5519.0810546875],
    [9, 105006.265625, 3025.29296875],
    [10, 28162.521484375, 1558.0169677734],
], dtype=float)

b_frac = np.array([
    [2, 0.8012156487, 0.0004167029],
    [3, 0.6545857787, 0.0005493056],
    [4, 0.5547605753, 0.0007846206],
    [5, 0.4851440191, 0.0012005592],
    [6, 0.4360441864, 0.0019477886],
    [7, 0.3960957825, 0.0033156750],
    [8, 0.3688789606, 0.0058940331],
    [9, 0.3460989892, 0.0106819842],
    [10, 0.3517312706, 0.0201896060],
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
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet p_T - Top1")

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

plt.savefig("chi2_top1_flavour_vs_jet_pT.png")
plt.show()
plt.close()