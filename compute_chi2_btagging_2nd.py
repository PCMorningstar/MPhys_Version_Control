import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Config
# -------------------------
fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

MEV_TO_GEV = 1.0 / 1000.0
bin_edges = np.linspace(0, 6, 60)  # GeV
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_widths  = np.diff(bin_edges)

# -------------------------
# Load branches (once)
# -------------------------
branches = [
    "sv1_ordered_jets_NOSYS",
    "ordered_jet_truth_flavour_NOSYS",
    "raw_chi2_minval_truthall_NOSYS",
    "selection_cuts_NOSYS",
    "jet_size_NOSYS",
    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

# -------------------------
# Region: selection_cuts==1 and >=2 jets
# -------------------------
evt_region = (arr["selection_cuts_NOSYS"] == 1) & (arr["jet_size_NOSYS"] >= 2)

sv_evt    = arr["sv1_ordered_jets_NOSYS"][evt_region]            # events -> jets (MeV)
truth_evt = arr["ordered_jet_truth_flavour_NOSYS"][evt_region]    # events -> jets
chi_evt   = arr["raw_chi2_minval_truthall_NOSYS"][evt_region]     # events -> [12] floats

w_evt = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[evt_region]

# -------------------------
# Denominator: total true b-jets (truth flavour only)
# -------------------------
w_jet = ak.broadcast_arrays(w_evt, sv_evt)[0]

sv_flat_mev = ak.to_numpy(ak.flatten(sv_evt))
truth_flat  = ak.to_numpy(ak.flatten(truth_evt))
w_flat      = ak.to_numpy(ak.flatten(w_jet))

mask_trueb = (np.abs(truth_flat) == 5) & np.isfinite(sv_flat_mev) & (sv_flat_mev >= 0.0)

sv_trueb = sv_flat_mev[mask_trueb] * MEV_TO_GEV
w_trueb  = w_flat[mask_trueb]

# -------------------------
# Numerator: χ²-selected jets that are true b (selected jet is b)
# -------------------------
chi_i = ak.values_astype(chi_evt[:, 0], int)  # chosen jet index for l+
chi_j = ak.values_astype(chi_evt[:, 1], int)  # chosen jet index for l-
evt_valid_chi = (chi_i >= 0) & (chi_j >= 0)

sv_v    = sv_evt[evt_valid_chi]
truth_v = truth_evt[evt_valid_chi]
w_v     = w_evt[evt_valid_chi]
chi_i_v = chi_i[evt_valid_chi]
chi_j_v = chi_j[evt_valid_chi]

# guard against out-of-range indices
njet_v = ak.num(sv_v)
evt_inrange = (
    (chi_i_v >= 0) & (chi_i_v < njet_v) &
    (chi_j_v >= 0) & (chi_j_v < njet_v)
)

sv_v    = sv_v[evt_inrange]
truth_v = truth_v[evt_inrange]
w_v     = w_v[evt_inrange]
chi_i_v = chi_i_v[evt_inrange]
chi_j_v = chi_j_v[evt_inrange]

evt_idx = np.arange(len(sv_v))

sv_sel_i_mev = ak.to_numpy(sv_v[evt_idx, chi_i_v])
sv_sel_j_mev = ak.to_numpy(sv_v[evt_idx, chi_j_v])

truth_sel_i  = ak.to_numpy(truth_v[evt_idx, chi_i_v])
truth_sel_j  = ak.to_numpy(truth_v[evt_idx, chi_j_v])

mask_i = (np.abs(truth_sel_i) == 5) & np.isfinite(sv_sel_i_mev) & (sv_sel_i_mev >= 0.0)
mask_j = (np.abs(truth_sel_j) == 5) & np.isfinite(sv_sel_j_mev) & (sv_sel_j_mev >= 0.0)

sv_corr = np.concatenate([sv_sel_i_mev[mask_i], sv_sel_j_mev[mask_j]]) * MEV_TO_GEV
w_corr  = np.concatenate([ak.to_numpy(w_v)[mask_i], ak.to_numpy(w_v)[mask_j]])

# -------------------------
# Weighted histograms + stat errors (sqrt(sumw2))
# -------------------------
# Top panel inputs
counts_true, _ = np.histogram(sv_trueb, bins=bin_edges, weights=w_trueb)
sumw2_true, _  = np.histogram(sv_trueb, bins=bin_edges, weights=w_trueb**2)

counts_corr, _ = np.histogram(sv_corr,  bins=bin_edges, weights=w_corr)
sumw2_corr, _  = np.histogram(sv_corr,  bins=bin_edges, weights=w_corr**2)

err_true = np.sqrt(sumw2_true)
err_corr = np.sqrt(sumw2_corr)

# Bottom panel efficiency
num = counts_corr
den = counts_true
num_err = err_corr
den_err = err_true

eff = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=(den > 0))

eff_err = np.zeros_like(eff, dtype=float)
m = (den > 0) & (num > 0)
eff_err[m] = eff[m] * np.sqrt((num_err[m] / num[m])**2 + (den_err[m] / den[m])**2)

# Remove trailing empty bins for efficiency panel
valid_bins = den > 0
x_eff   = bin_centers[valid_bins]
y_eff   = eff[valid_bins]
yerr_eff = eff_err[valid_bins]

# Cut x-axis at the last populated bin
xmax = x_eff[-1] + bin_widths[0] / 2.0

# -------------------------
# Combined figure (top + bottom)
# -------------------------
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(9, 7),
    gridspec_kw={"height_ratios": [3, 1]},
    sharex=True
)

# ---- Top: entry counts (log scale)
ax_top.bar(
    bin_centers,
    counts_true,
    width=bin_widths,
    color="lightsteelblue",
    label="Total Truth b-jets",
    align="center",
)

ax_top.step(
    bin_centers,
    counts_corr,
    where="mid",
    color="red",
    linewidth=2,
    label=r"Correctly Tagged by $\chi^2$ (selected jet is b)",
)

ax_top.set_yscale("log")
ax_top.set_ylabel("Entry Count")
ax_top.set_title(r"$\chi^2$ b-tagging Performance vs. SV Mass")
ax_top.legend()
ax_top.set_xlim(0, xmax)

# ---- Bottom: efficiency (black) with errors
ax_bot.errorbar(
    x_eff,
    y_eff,
    yerr=yerr_eff,
    fmt='-',
    color='black',
    linewidth=1.6,
    capsize=2
)

ax_bot.set_ylim(0.0, 1.05)
ax_bot.set_xlabel("SV1 Mass [GeV]")
ax_bot.set_ylabel("Efficiency\n(Correct/Truth)")
ax_bot.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig("chi2_btagging_performance_sv1mass_top_and_eff_geq2jets.png", dpi=300)
plt.show()
plt.close(fig)