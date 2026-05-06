import numpy as np
import matplotlib.pyplot as plt

################### Top2 - SV INVARIANT MASS ###################

b_yield = np.array([
    [-0.25, 2999.617432, 13.352517],
    [0.25, 572.033081, 5.817156],
    [0.75, 2821.635742, 12.959229],
    [1.25, 4745.568359, 16.804943],
    [1.75, 6529.450195, 19.717424],
    [2.25, 7413.509277, 21.011988],
    [2.75, 7156.057617, 20.652710],
    [3.25, 5711.020020, 18.464701],
    [3.75, 3738.709229, 14.937316],
    [4.25, 2086.650879, 11.159458],
    [4.75, 985.634888, 7.676378],
    [5.25, 422.098816, 5.031370],
], dtype=float)

c_yield = np.array([
    [-0.25, 128.353928, 2.768116],
    [0.25, 14.538824, 0.931689],
    [0.75, 64.396484, 1.978572],
    [1.25, 74.372269, 2.096114],
    [1.75, 76.405037, 2.143252],
    [2.25, 63.785614, 1.942382],
    [2.75, 56.901176, 1.831743],
    [3.25, 38.278896, 1.501808],
    [3.75, 24.019564, 1.190627],
    [4.25, 13.837051, 0.917529],
    [4.75, 8.696269, 0.735408],
    [5.25, 3.714579, 0.463115],
], dtype=float)

l_yield = np.array([
    [-0.25, 3112.727051, 13.622320],
    [0.25, 255.747772, 3.903667],
    [0.75, 967.386780, 7.605004],
    [1.25, 1187.280762, 8.404314],
    [1.75, 1287.196899, 8.756034],
    [2.25, 1188.187012, 8.412517],
    [2.75, 990.756470, 7.681273],
    [3.25, 725.148682, 6.579440],
    [3.75, 460.334381, 5.248496],
    [4.25, 255.882584, 3.902858],
    [4.75, 130.329529, 2.794063],
    [5.25, 71.386383, 2.059440],
], dtype=float)

b_frac = np.array([
    [-0.25, 0.480654, 0.001543],
    [0.25, 0.679116, 0.003922],
    [0.75, 0.732242, 0.001744],
    [1.25, 0.789977, 0.001282],
    [1.75, 0.827240, 0.001039],
    [2.25, 0.855522, 0.000922],
    [2.75, 0.872295, 0.000899],
    [3.25, 0.882086, 0.000979],
    [3.75, 0.885307, 0.001199],
    [4.25, 0.885536, 0.001601],
    [4.75, 0.876384, 0.002404],
    [5.25, 0.848952, 0.003915],
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

ax.set_ylabel("Weighted Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. SV Mass - Top2")
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

plt.savefig("chi2_top2_flavour_vs_SV_invariant_mass.png")
plt.show()
plt.close()