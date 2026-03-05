import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# SV1 is in MeV in your ntuples; we convert to GeV before histogramming
bins = np.linspace(0, 6, 60)

with uproot.open(fname) as f:
    arr = f[tree].arrays([
        "sv1_ordered_jets_NOSYS",
        "ordered_jet_truth_flavour_NOSYS",
        "raw_chi2_minval_truthall_NOSYS",
        "selection_cuts_NOSYS",
        "jet_size_NOSYS",
        "weight_mc_NOSYS",
        "weight_pileup_NOSYS",
        "weight_leptonSF_tight_NOSYS",
        "weight_jvt_effSF_NOSYS",
    ], library="ak")

# -------------------------
# Region selection: selection_cuts==1 and exactly 2 jets
# -------------------------
region_mask = (arr["selection_cuts_NOSYS"] == 1) & (arr["jet_size_NOSYS"] >= 0)

sv    = arr["sv1_ordered_jets_NOSYS"][region_mask]          # events -> jets (MeV)
truth = arr["ordered_jet_truth_flavour_NOSYS"][region_mask]  # events -> jets
chi   = arr["raw_chi2_minval_truthall_NOSYS"][region_mask]   # events -> 12 floats

w_event = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[region_mask]

# -------------------------
# Total true b-jets distribution (jet-level, truth flavour only)
# -------------------------
w_jet = ak.broadcast_arrays(w_event, sv)[0]

sv_flat    = ak.to_numpy(ak.flatten(sv))
truth_flat = ak.to_numpy(ak.flatten(truth))
w_flat     = ak.to_numpy(ak.flatten(w_jet))

mask_trueb = (np.abs(truth_flat) == 5) & np.isfinite(sv_flat) & (sv_flat >= 0)

sv_trueb = (sv_flat[mask_trueb] / 1000.0)   # MeV -> GeV
w_trueb  = w_flat[mask_trueb]

# -------------------------
# χ² selected jets; "correctly tagged" defined using truth flavour only
# (selected jet is b/bbar => abs(truth)==5)
# -------------------------
chi_i = ak.values_astype(chi[:, 0], int)  # best jet index for l+
chi_j = ak.values_astype(chi[:, 1], int)  # best jet index for l-
valid = (chi_i >= 0) & (chi_j >= 0)

sv_v    = sv[valid]
truth_v = truth[valid]
w_v     = w_event[valid]
chi_i_v = chi_i[valid]
chi_j_v = chi_j[valid]

# safety: ensure indices in range
jet_counts = ak.num(sv_v)
valid_idx = (
    (chi_i_v >= 0) & (chi_i_v < jet_counts) &
    (chi_j_v >= 0) & (chi_j_v < jet_counts)
)

sv_v    = sv_v[valid_idx]
truth_v = truth_v[valid_idx]
w_v     = w_v[valid_idx]
chi_i_v = chi_i_v[valid_idx]
chi_j_v = chi_j_v[valid_idx]

evt_idx = np.arange(len(sv_v))

sv_sel_i = sv_v[evt_idx, chi_i_v]
sv_sel_j = sv_v[evt_idx, chi_j_v]

truth_sel_i = truth_v[evt_idx, chi_i_v]
truth_sel_j = truth_v[evt_idx, chi_j_v]

mask_i = (np.abs(ak.to_numpy(truth_sel_i)) == 5)
mask_j = (np.abs(ak.to_numpy(truth_sel_j)) == 5)

sv_corr = ak.to_numpy(ak.concatenate([sv_sel_i[mask_i], sv_sel_j[mask_j]])) / 1000.0  # MeV -> GeV
w_corr  = ak.to_numpy(ak.concatenate([w_v[mask_i],     w_v[mask_j]]))

m_corr = np.isfinite(sv_corr) & (sv_corr >= 0)
sv_corr = sv_corr[m_corr]
w_corr  = w_corr[m_corr]

# -------------------------
# Histograms with errors (sqrt(sumw2))
# -------------------------
counts_true, edges = np.histogram(sv_trueb, bins=bins, weights=w_trueb)
sumw2_true, _      = np.histogram(sv_trueb, bins=bins, weights=w_trueb**2)

counts_corr, _     = np.histogram(sv_corr,  bins=bins, weights=w_corr)
sumw2_corr, _      = np.histogram(sv_corr,  bins=bins, weights=w_corr**2)

err_true = np.sqrt(sumw2_true)
err_corr = np.sqrt(sumw2_corr)

centers = 0.5 * (edges[1:] + edges[:-1])
widths  = np.diff(edges)

# -------------------------
# Plot (log y-axis)
# -------------------------
plt.figure(figsize=(7, 5))

plt.bar(
    centers,
    counts_true,
    width=widths,
    color="lightsteelblue",
    label="Total Truth b-jets",
    align="center",
)

plt.step(
    centers,
    counts_corr,
    where="mid",
    color="red",
    linewidth=2,
    label="Correctly Tagged by χ² (selected jet is b)",
)


plt.yscale("log")
plt.xlabel("SV1 Mass [GeV]")
plt.ylabel("Entry Count")
plt.title(r"$\chi^2$ b-tagging Performance vs SV1 Mass")
plt.legend()
plt.tight_layout()
plt.savefig("chi2_btagging_performance_v_sv1_mass_geq2jets.png", dpi=300)
plt.show()
plt.close()