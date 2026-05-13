import numpy as np
import matplotlib.pyplot as plt

# =================================================
# Combined - Jet pT region
# =================================================
# =================================================
# Combined b yields
# =================================================
b_yield = np.array([
    [2, 90388.6289062500, 73.3868731835],
    [3, 59728.4082031250, 59.4709739458],
    [4, 26956.4189453125, 39.8509634984],
    [5, 10109.8955078125, 24.3602055627],
    [6, 3381.6485595703, 14.0749022842],
    [7, 1013.3573608399, 7.6714533313],
    [8, 285.0007781982, 4.0778645300],
    [9, 79.9205093384, 2.1559697411],
    [10, 21.5790929794, 1.0995295338],
], dtype=float)

# =================================================
# Combined c yields
# =================================================
c_yield = np.array([
    [2, 1137.6088867188, 8.2260415515],
    [3, 1477.6850585938, 9.3375647306],
    [4, 1050.7073364257, 7.8453365054],
    [5, 554.4685668945, 5.6899718346],
    [6, 241.0328826905, 3.7540412920],
    [7, 82.5691528320, 2.1849104656],
    [8, 28.1573505402, 1.2681832551],
    [9, 9.1815319061, 0.7322370432],
    [10, 2.6427168846, 0.3830029966],
], dtype=float)

# =================================================
# Combined light yields
# =================================================
l_yield = np.array([
    [2, 21237.0781250000, 35.5854931265],
    [3, 30171.1982421875, 42.2502631472],
    [4, 20603.4218750000, 34.8252933413],
    [5, 10213.2304687500, 24.4547443817],
    [6, 4154.5002441406, 15.5747466818],
    [7, 1491.1445312500, 9.2932482768],
    [8, 484.3770599365, 5.3142288738],
    [9, 145.9781875611, 2.9069248056],
    [10, 39.2895565033, 1.4957648092],
], dtype=float)

# =================================================
# Combined b fractions
# =================================================
b_frac = np.array([
    [2, 0.8008047444, 0.0002899770],
    [3, 0.6525722894, 0.0003831738],
    [4, 0.5545384882, 0.0005471070],
    [5, 0.4842462225, 0.0008374130],
    [6, 0.4348167337, 0.0013594541],
    [7, 0.3916982735, 0.0023107618],
    [8, 0.3628539473, 0.0041141126],
    [9, 0.3400529832, 0.0074579249],
    [10, 0.3403027405, 0.0141366203],
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
bin_edges = np.array(
    [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 900],
    dtype=float
)

# Optional consistency check
expected_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
pt = b_yield[:, 0]
if not np.allclose(pt, expected_centres):
    raise ValueError("Bin centres in data do not match the supplied bin edges.")

# =================================================
# Helpers
# =================================================
def plot_with_errorbars(ax, edges, x, y, yerr, key, label):
    ax.stairs(
        y,
        edges,
        color=COLORS[key],
        linewidth=2.0,
        fill=False,
        label=label,
        zorder=2,
    )

    ax.errorbar(
        x,
        y,
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
minor_ticks = np.array(
    [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 450, 540, 630, 720, 810, 900],
    dtype=float
)

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
plot_with_errorbars(ax, bin_edges, pt, b, b_err, "b", "b")
plot_with_errorbars(ax, bin_edges, pt, c, c_err, "c", "c")
plot_with_errorbars(ax, bin_edges, pt, l, l_err, "l", "l")

ax.set_ylabel("Events")
ax.set_title(r"$\chi^2$ Flavour Composition vs. Jet $p_T$ - Combined")
ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()
ax.legend(loc="upper right", frameon=False)

# =================================================
# Bottom: b fraction
# =================================================
rax.stairs(
    frac_b,
    bin_edges,
    color="black",
    linewidth=2.0,
    fill=False,
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
rax.set_xlim(bin_edges[0], bin_edges[-1])

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{int(x)}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("chi2_combined_flavour_vs_jet_pT.png")
plt.show()
plt.close()