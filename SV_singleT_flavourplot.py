import numpy as np
import matplotlib.pyplot as plt

################### Top1 - SV INVARIANT MASS ###################

# =================================================
# Combined b yields
# =================================================
b_yield = np.array([
    [-0.25, 5999.538330, 18.882497],
    [ 0.25, 1145.986145,  8.237840],
    [ 0.75, 5647.082031, 18.342605],
    [ 1.25, 9481.313476, 23.751360],
    [ 1.75,13073.138672, 27.898680],
    [ 2.25,14820.346679, 29.707186],
    [ 2.75,14316.632324, 29.207980],
    [ 3.25,11424.116700, 26.113555],
    [ 3.75, 7488.186280, 21.144381],
    [ 4.25, 4178.453613, 15.789378],
    [ 4.75, 1970.867249, 10.859252],
    [ 5.25,  841.976593,  7.107015],
], dtype=float)

# =================================================
# Combined c yields
# =================================================
c_yield = np.array([
    [-0.25, 265.521607,  3.978835],
    [ 0.25,  31.259203,  1.367419],
    [ 0.75, 128.244426,  2.777355],
    [ 1.25, 149.327256,  2.973750],
    [ 1.75, 149.537621,  2.992203],
    [ 2.25, 129.569481,  2.775463],
    [ 2.75, 107.207969,  2.512585],
    [ 3.25,  77.894836,  2.145470],
    [ 3.75,  48.554190,  1.694569],
    [ 4.25,  27.131262,  1.280203],
    [ 4.75,  15.507302,  0.971417],
    [ 5.25,   7.853784,  0.668558],
], dtype=float)

# =================================================
# Combined light yields
# =================================================
l_yield = np.array([
    [-0.25, 6216.005371, 19.251798],
    [ 0.25,  507.331680,  5.492514],
    [ 0.75, 1931.749573, 10.734381],
    [ 1.25, 2383.779908, 11.913478],
    [ 1.75, 2563.450805, 12.362791],
    [ 2.25, 2381.116944, 11.910453],
    [ 2.75, 1983.589478, 10.882947],
    [ 3.25, 1447.239136,  9.298738],
    [ 3.75,  908.239044,  7.363358],
    [ 4.25,  507.112915,  5.503574],
    [ 4.75,  262.946839,  3.960775],
    [ 5.25,  144.515327,  2.932768],
], dtype=float)

# =================================================
# Combined b fractions
# =================================================
b_frac = np.array([
    [-0.25, 0.4801940370, 0.0010917182],
    [ 0.25, 0.6780781887, 0.0027722780],
    [ 0.75, 0.7323645218, 0.0012334135],
    [ 1.25, 0.7892193875, 0.0009081336],
    [ 1.75, 0.8277945785, 0.0007353628],
    [ 2.25, 0.8551403915, 0.0006533566],
    [ 2.75, 0.8725686428, 0.0006374352],
    [ 3.25, 0.8819554912, 0.0006931165],
    [ 3.75, 0.8866997354, 0.0008474897],
    [ 4.25, 0.8847711466, 0.0011320208],
    [ 4.75, 0.8762469308, 0.0017001119],
    [ 5.25, 0.8474544022, 0.0027711955],
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
ax.set_ylim(1e0, 1e5)
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
rax.set_ylim(0.0, 1.0)
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