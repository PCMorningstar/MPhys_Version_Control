import numpy as np
import matplotlib.pyplot as plt

################### Top1 - CORRECTED ###################

b_yield = np.array([
    [2, 67642488.0, 77678.84375],
    [3, 44680776.0, 62938.46875],
    [4, 20200750.0, 42213.26953125],
    [5, 7565092.5, 25784.783203125],
    [6, 2530426.0, 14881.904296875],
    [7, 750748.9375, 8084.57421875],
    [8, 208508.703125, 4244.4213867188],
    [9, 58892.671875, 2274.7270507812],
    [10, 15808.669921875, 1155.8424072266],
], dtype=float)

c_yield = np.array([
    [2, 849311.25, 8707.837890625],
    [3, 1097127.125, 9845.3125],
    [4, 782518.6875, 8279.9853515625],
    [5, 410817.5625, 5994.2739257812],
    [6, 179672.3125, 3956.7790527344],
    [7, 62202.0078125, 2326.9797363281],
    [8, 22674.7421875, 1394.3541259766],
    [9, 6652.6225585938, 764.5280151367],
    [10, 1735.7578125, 384.7969360352],
], dtype=float)

l_yield = np.array([
    [2, 15917842.0, 37691.3046875],
    [3, 22623738.0, 44761.1640625],
    [4, 15404992.0, 36843.55859375],
    [5, 7652224.5, 25897.6171875],
    [6, 3111496.5, 16508.87890625],
    [7, 1123668.875, 9865.087890625],
    [8, 366564.9375, 5667.79296875],
    [9, 110425.4375, 3085.5478515625],
    [10, 29997.390625, 1594.6661376953],
], dtype=float)

b_frac = np.array([
    [2, 0.8013594747, 0.0004102347],
    [3, 0.6532118320, 0.0005416408],
    [4, 0.5551447272, 0.0007734524],
    [5, 0.4840688407, 0.0011842830],
    [6, 0.4346620440, 0.0019222875],
    [7, 0.3876594603, 0.0032635877],
    [8, 0.3488235474, 0.0057429288],
    [9, 0.3346731067, 0.0105128437],
    [10, 0.3325213492, 0.0198742946],
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

ax.set_ylabel("Weighted Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet Multiplicity - Top2")

ax.set_xlim(edges[0], edges[-1])
ax.set_ylim(1e3, 1e8)
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