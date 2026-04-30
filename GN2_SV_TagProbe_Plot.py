import numpy as np
import matplotlib.pyplot as plt

# ================================================= 
# b - Events
# ================================================= 
b_yield = np.array([
    [-0.25,   59613.984375,   2300.355957],
    [ 0.25,   38874.960938,   1847.597656],
    [ 0.75,  389292.562500,   5906.194824],
    [ 1.25, 1463274.500000,  11407.054688],
    [ 1.75, 3339577.500000,  17253.572266],
    [ 2.25, 5080530.000000,  21285.535156],
    [ 2.75, 5715139.000000,  22584.531250],
    [ 3.25, 4958822.000000,  21053.396484],
    [ 3.75, 3384584.250000,  17387.451172],
    [ 4.25, 1927614.250000,  13121.601562],
    [ 4.75,  909799.500000,   9031.028320],
    [ 5.25,  374631.187500,   5805.943848],
], dtype=float)


# ================================================= 
# non-b - Events
# ================================================= 
nonb_yield = np.array([
    [-0.25, 194.777802, 118.028427],
    [ 0.25,  53.177876,  53.177876],
    [ 0.75, 280.261353, 162.934464],
    [ 1.25, 890.433167, 273.773407],
    [ 1.75,1166.223755, 311.399109],
    [ 2.25,2065.873535, 428.523071],
    [ 2.75,1253.535034, 338.515320],
    [ 3.25,2162.624756, 443.213593],
    [ 3.75, 874.067261, 294.135834],
    [ 4.25, 244.068115, 142.938339],
    [ 4.75,   0.000000,   0.000000],   # undefined ratio
    [ 5.25,  79.905853,  79.905853],
], dtype=float)

# ================================================= 
# 2022 real data (scaled to MC normalisation)
# Columns: SV centre, N_data_scaled, err_scaled
# ================================================= 
data_yield_scaled = np.array([
    [-0.25,    54659.784816,   9662.576125],
    [ 0.25,    42702.956888,   8540.591378],
    [ 0.75,   406532.149569,  26351.565309],
    [ 1.25,  1549263.275880,  51442.442740],
    [ 1.75,  3440150.206860,  76656.268099],
    [ 2.25,  5071403.159964,  93072.855441],
    [ 2.75,  5614584.771574,  97930.459295],
    [ 3.25,  4789563.644507,  90449.661099],
    [ 3.75,  3441858.325136,  76675.296588],
    [ 4.25,  1907968.113735,  57087.960238],
    [ 4.75,   932632.578424,  39912.989759],
    [ 5.25,   399699.676467,  26129.185255],
], dtype=float)

# ================================================= 
# Data / MC ratio
# ================================================= 
data_mc_ratio = np.array([
    [-0.25, 0.913999, 0.163525],
    [ 0.25, 1.096990, 0.220483],
    [ 0.75, 1.043602, 0.067755],
    [ 1.25, 1.058281, 0.035272],
    [ 1.75, 1.029720, 0.023207],
    [ 2.25, 0.998184, 0.018344],
    [ 2.75, 0.982369, 0.017300],
    [ 3.25, 0.965482, 0.018335],
    [ 3.75, 1.016665, 0.024498],
    [ 4.25, 0.989606, 0.029882],
    [ 4.75, 1.025084, 0.043972],
    [ 5.25, 1.066725, 0.071075],
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
    "nb": "#E69F00",
}

# =================================================
# Variable-width SV-mass binning
# =================================================
bin_edges = np.array([
    -0.5, 0.0, 0.5,
     1.0,
     1.5,
     2.0,
     2.5,
     3.0,
     3.5,
     4.0,
     4.5,
     5.0,
     5.5
], dtype=float)

bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# =================================================
# Extract
# =================================================
b = b_yield[:, 1]
nb = nonb_yield[:, 1]
b_err = b_yield[:, 2]
nb_err = nonb_yield[:, 2]

ratio = data_mc_ratio[:, 1]
ratio_err = data_mc_ratio[:, 2]

total = b + nb

nb_step = np.r_[nb, nb[-1]]
total_step = np.r_[total, total[-1]]
ratio_step = np.r_[ratio, ratio[-1]]

major_ticks = np.arange(-0.5, 5.5 + 0.5, 0.5)
minor_ticks = np.arange(-0.5, 5.5 + 0.25, 0.25)

# =================================================
# Figure with ratio panel
# =================================================
fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)

ax = fig.add_subplot(gs[0])
rax = fig.add_subplot(gs[1], sharex=ax)

# =================================================
# Top: stacked yields
# =================================================

data_step = np.r_[data_yield_scaled[:, 1], data_yield_scaled[-1, 1]]

ax.step(
    bin_edges,
    data_step,
    where="post",
    color="black",
    linewidth=1.4,
    zorder=2,
)

ax.errorbar(
    bin_centres,
    data_yield_scaled[:, 1],
    yerr=data_yield_scaled[:, 2],
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
    label="Data",
)

# b stacked on top of non-b
ax.fill_between(
    bin_edges,
    nb_step,
    total_step,
    step="post",
    color=COLORS["b"],
    alpha=0.7,
    label="b",
)

# non-b bottom
ax.fill_between(
    bin_edges,
    1e0,
    nb_step,
    step="post",
    color=COLORS["nb"],
    alpha=0.7,
    label="non-b",
)

ax.set_ylabel("Events")
ax.set_title(r"Probe-Jet Flavour Composition - GN2 WP 65%")
ax.set_xlim(bin_edges[0], bin_edges[-1])
ax.set_ylim(1e0, 1e8)
ax.set_yscale("log", base=10)
ax.tick_params(labelbottom=False)

ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)
ax.minorticks_on()

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0.0,
    frameon=False,
)

# =================================================
# Bottom: Data / MC ratio
# =================================================
rax.step(
    bin_edges,
    ratio_step,
    where="post",
    color="black",
    linewidth=2.0,
    zorder=2,
)

rax.errorbar(
    bin_centres,
    ratio,
    yerr=ratio_err,
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

rax.set_ylabel("Data / MC")
rax.set_xlabel(r"SV mass [GeV]")
rax.set_ylim(0.5, 1.5)
#rax.set_yscale("log", base=10)
rax.set_xlim(bin_edges[0], bin_edges[-1])

rax.set_xticks(major_ticks)
rax.set_xticklabels([f"{x:.1f}" for x in major_ticks])
rax.set_xticks(minor_ticks, minor=True)

rax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.7)
rax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.35)

plt.savefig("DatavMC_probe_jet_flavour_composition_GN2_WP_65.png")
plt.show()
plt.close()