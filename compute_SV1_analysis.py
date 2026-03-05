import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

fname = "your_ntuple.root"
tree  = "nominal"

bins = np.linspace(0,6,60)

# -------------------------
# Load branches
# -------------------------
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
        "weight_jvt_effSF_NOSYS"
    ], library="ak")

# -------------------------
# Region selection
# -------------------------
region_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] >= 2)
)

sv     = arr["sv1_ordered_jets_NOSYS"][region_mask]
truth  = arr["ordered_jet_truth_flavour_NOSYS"][region_mask]
chi    = arr["raw_chi2_minval_truthall_NOSYS"][region_mask]

# -------------------------
# Event weights
# -------------------------
w_event = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[region_mask]

# broadcast event weights to jets
w_jet = ak.broadcast_arrays(w_event, sv)[0]

# -------------------------
# TRUE B-JETS DISTRIBUTION
# -------------------------
sv_flat   = ak.to_numpy(ak.flatten(sv))
truth_flat = ak.to_numpy(ak.flatten(truth))
w_flat     = ak.to_numpy(ak.flatten(w_jet))

mask_trueb = (np.abs(truth_flat) == 5)

sv_trueb = sv_flat[mask_trueb]
w_trueb  = w_flat[mask_trueb]

# -------------------------
# CHI2 SELECTED JETS
# -------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:,0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:,1], int))

valid = (chi_i >= 0) & (chi_j >= 0)

sv = sv[valid]
truth = truth[valid]
w_event = w_event[valid]

sv_i = ak.to_numpy(sv[ak.local_index(sv)[:,0]*0 + np.arange(len(sv)), chi_i[valid]])
sv_j = ak.to_numpy(sv[ak.local_index(sv)[:,0]*0 + np.arange(len(sv)), chi_j[valid]])

truth_i = ak.to_numpy(truth[ak.local_index(truth)[:,0]*0 + np.arange(len(truth)), chi_i[valid]])
truth_j = ak.to_numpy(truth[ak.local_index(truth)[:,0]*0 + np.arange(len(truth)), chi_j[valid]])

mask_i = (np.abs(truth_i) == 5)
mask_j = (np.abs(truth_j) == 5)

sv_corr = np.concatenate([sv_i[mask_i], sv_j[mask_j]])
w_corr  = np.concatenate([w_event[mask_i], w_event[mask_j]])

# -------------------------
# HISTOGRAMS (weighted)
# -------------------------
counts_true, edges = np.histogram(
    sv_trueb,
    bins=bins,
    weights=w_trueb
)

sumw2_true, _ = np.histogram(
    sv_trueb,
    bins=bins,
    weights=w_trueb**2
)

counts_corr, _ = np.histogram(
    sv_corr,
    bins=bins,
    weights=w_corr
)

sumw2_corr, _ = np.histogram(
    sv_corr,
    bins=bins,
    weights=w_corr**2
)

err_true = np.sqrt(sumw2_true)
err_corr = np.sqrt(sumw2_corr)

centers = 0.5*(edges[1:] + edges[:-1])

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(7,5))

plt.bar(
    centers,
    counts_true,
    width=np.diff(edges),
    color="lightsteelblue",
    label="Total Truth b-jets",
    align="center"
)

plt.step(
    centers,
    counts_corr,
    where="mid",
    color="red",
    linewidth=2,
    label="Correctly Tagged by χ²"
)

plt.errorbar(
    centers,
    counts_corr,
    yerr=err_corr,
    fmt="none",
    color="black",
    capsize=2
)

plt.yscale("log")

plt.xlabel("SV1 Mass [GeV]")
plt.ylabel("Entry Count")
plt.title("Chi2 b-tagging Performance vs SV1 Mass")

plt.legend()
plt.tight_layout()
plt.savefig("sv1_analysis_performance_v_sv1_mass.png", dpi=300)
plt.show()
plt.close()