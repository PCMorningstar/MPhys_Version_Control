import numpy as np
import matplotlib.pyplot as plt

################### Top2 - CORRECTED ###################

b_yield = np.array([
    [2, 45181.9804687500, 51.8857955933],
    [3, 29844.6445312500, 42.0399246216],
    [4, 13493.1406250000, 28.1964626312],
    [5, 5053.1230468750, 17.2230129242],
    [6, 1690.2047119141, 9.9404077530],
    [7, 501.4647216797, 5.4001126289],
    [8, 139.2739410400, 2.8350725174],
    [9, 39.3375167847, 1.5194100142],
    [10, 10.5594425201, 0.7720479965],
], dtype=float)

c_yield = np.array([
    [2, 567.2996826172, 5.8164238930],
    [3, 732.8289794922, 6.5762023926],
    [4, 522.6854858398, 5.5306377411],
    [5, 274.4066772461, 4.0038909912],
    [6, 120.0125885010, 2.6429409981],
    [7, 41.5479965210, 1.5543123484],
    [8, 15.1456537247, 0.9313625097],
    [9, 4.4436368942, 0.5106685162],
    [10, 1.1594040394, 0.2570261359],
], dtype=float)

l_yield = np.array([
    [2, 10632.3652343750, 25.1760082245],
    [3, 15111.5869140625, 29.8983364105],
    [4, 10289.8046875000, 24.6097526550],
    [5, 5111.3232421875, 17.2983798981],
    [6, 2078.3320312500, 11.0271501541],
    [7, 750.5576171875, 6.5894112587],
    [8, 244.8480072021, 3.7858173847],
    [9, 73.7589645386, 2.0609998703],
    [10, 20.0368366241, 1.0651613474],
], dtype=float)

b_frac = np.array([
    [2, 0.8013598323, 0.0004102350],
    [3, 0.6532121301, 0.0005416409],
    [4, 0.5551446080, 0.0007734525],
    [5, 0.4840687811, 0.0011842830],
    [6, 0.4346620142, 0.0019222874],
    [7, 0.3876594305, 0.0032635874],
    [8, 0.3488235772, 0.0057429294],
    [9, 0.3346731067, 0.0105128441],
    [10, 0.3325213492, 0.0198742949],
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
ax.set_ylim(1e0, 1e5)
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