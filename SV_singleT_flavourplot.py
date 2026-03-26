import numpy as np
import matplotlib.pyplot as plt

# ================================================= 
# T2 flavour analysis - SV invariant mass region
# ================================================= 
# Weighted yield and error calculation
# -------------------------------------------------
b_yield = np.array([
    [0.25, 15438121.0000000000, 36933.0781250000],
    [0.75,  9748399.0000000000, 29388.4042968750],
    [1.25, 15264166.0000000000, 36806.3125000000],
    [1.75, 20020576.0000000000, 42145.2343750000],
    [2.25, 21836464.0000000000, 44030.0234375000],
    [2.75, 20657042.0000000000, 42842.4375000000],
    [3.25, 16197140.0000000000, 37968.8320312500],
    [3.75, 10556579.0000000000, 30648.3613281250],
    [4.25,  5824061.0000000000, 22769.4160156250],
    [4.75,  2754722.0000000000, 15683.3906250000],
    [5.25,  1188169.2500000000, 10304.0673828125],
], dtype=float)

c_yield = np.array([
    [0.25, 1084951.7500000000,  9778.2265625000],
    [0.75,  365692.3437500000,  5692.8535156250],
    [1.25,  418672.3750000000,  6070.1074218750],
    [1.75,  413598.6562500000,  6039.4238281250],
    [2.25,  330130.0937500000,  5397.0434570312],
    [2.75,  271796.6875000000,  4894.1005859375],
    [3.25,  191849.0000000000,  4110.6689453125],
    [3.75,  120138.3281250000,  3250.6762695312],
    [4.25,   68267.7109375000,  2469.4719238281],
    [4.75,   35315.4531250000,  1786.8688964844],
    [5.25,   15067.5937500000,  1131.9222412109],
], dtype=float)

l_yield = np.array([
    [0.25, 26458062.0000000000, 48327.9765625000],
    [0.75,  5274568.0000000000, 21604.4042968750],
    [1.25,  6329407.0000000000, 23651.8691406250],
    [1.75,  6760142.5000000000, 24445.1054687500],
    [2.25,  6221953.5000000000, 23480.4218750000],
    [2.75,  5163587.0000000000, 21380.5605468750],
    [3.25,  3745751.5000000000, 18226.3828125000],
    [3.75,  2328120.5000000000, 14386.6806640625],
    [4.25,  1259192.5000000000, 10576.2089843750],
    [4.75,   607370.1250000000,  7341.7832031250],
    [5.25,   284026.6875000000,  5016.5776367188],
], dtype=float)

# -------------------------------------------------
# Weighted event fraction and error calculation
# -------------------------------------------------
b_frac = np.array([
    [0.25, 0.3591836393, 0.0006877458],
    [0.75, 0.6334794164, 0.0011557733],
    [1.25, 0.6934397817, 0.0009243817],
    [1.75, 0.7362045050, 0.0007948714],
    [2.25, 0.7691996098, 0.0007445064],
    [2.75, 0.7916872501, 0.0007482519],
    [3.25, 0.8044373989, 0.0008326081],
    [3.75, 0.8117423654, 0.0010219604],
    [4.25, 0.8143807650, 0.0013707097],
    [4.75, 0.8108306527, 0.0020036612],
    [5.25, 0.7988962531, 0.0030939195],
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
major_ticks = np.arange(0, 5.5 + 0.5, 0.5)
minor_ticks = np.arange(0, 5.5 + BIN_WIDTH, BIN_WIDTH)

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
ax.set_title(r"$\chi^2$ Flavour Composition vs. SV Mass - Top2")
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