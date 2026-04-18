import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# Stored leading-jet pT region flags from ntuple
# NOTE:
# These are still leading-jet pT bins, not probe-jet pT bins.
# Keep only if this is intentionally what you want.
# -------------------------------------------------
pt_regions = [
    ("0to30",    "jet_pt_region_0to30_GeV_NOSYS"),
    ("30to60",   "jet_pt_region_30to60_GeV_NOSYS"),
    ("60to90",   "jet_pt_region_60to90_GeV_NOSYS"),
    ("90to120",  "jet_pt_region_90to120_GeV_NOSYS"),
    ("120to150", "jet_pt_region_120to150_GeV_NOSYS"),
    ("150to180", "jet_pt_region_150to180_GeV_NOSYS"),
    ("180to210", "jet_pt_region_180to210_GeV_NOSYS"),
    ("210to240", "jet_pt_region_210to240_GeV_NOSYS"),
    ("240to270", "jet_pt_region_240to270_GeV_NOSYS"),
    ("270to300", "jet_pt_region_270to300_GeV_NOSYS"),
    ("300to360", "jet_pt_region_300to360_GeV_NOSYS"),
    ("360to900", "jet_pt_region_360to900_GeV_NOSYS"),
]

# -------------------------------------------------
# Helpers
# -------------------------------------------------
def is_b_flavour(f):
    try:
        return abs(int(f)) == 5
    except Exception:
        return False

def efficiency_and_error(sumw_num, sumw2_num, sumw_den, sumw2_den):
    if sumw_den <= 0.0:
        return 0.0, 0.0

    eff = sumw_num / sumw_den
    var = (
        ((1.0 - eff) ** 2) * sumw2_num
        + (eff ** 2) * (sumw2_den - sumw2_num)
    ) / (sumw_den ** 2)

    if var < 0.0:
        var = 0.0

    return eff, np.sqrt(var)

# -------------------------------------------------
# Load branches
# -------------------------------------------------
branches = [
    "selection_cuts_NOSYS",
    "jet_size_NOSYS",
    "ordered_jet_truth_flavour_NOSYS",
    "raw_chi2_minval_truthall_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_65_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_70_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_77_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_85_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_90_NOSYS",
    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
] + [b for _, b in pt_regions]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

# -------------------------------------------------
# Base event selection
# -------------------------------------------------
base_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] == 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi   = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]

# Choose the WP branch you actually want here
wp_probe = arr["jet_select_GN2v01_FixedCutBEff_65_NOSYS"][base_mask]

w_event = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[base_mask]

region_flags = {
    label: arr[branch][base_mask]
    for label, branch in pt_regions
}

# -------------------------------------------------
# Extract chi2-selected pair
# Convention:
#   chi[:, 0] -> Top1
#   chi[:, 1] -> Top2 = PROBE
# -------------------------------------------------
chi_top1 = ak.to_numpy(ak.values_astype(chi[:, 0], int))
chi_top2 = ak.to_numpy(ak.values_astype(chi[:, 1], int))
n_jets   = ak.to_numpy(ak.num(truth, axis=1))

valid = (
    (chi_top1 >= 0) & (chi_top2 >= 0)
    & (chi_top1 < n_jets) & (chi_top2 < n_jets)
    & (chi_top1 != chi_top2)
)

truth    = truth[valid]
wp_probe = wp_probe[valid]
w_event  = ak.to_numpy(w_event[valid])
chi_top1 = chi_top1[valid]
chi_top2 = chi_top2[valid]

for label in region_flags:
    region_flags[label] = ak.to_numpy(region_flags[label][valid])

idx = np.arange(len(truth))

# Top1 = probe jet
probe_truth = ak.to_numpy(truth[idx, chi_top1])
probe_wp    = ak.to_numpy(ak.values_astype(wp_probe[idx, chi_top1], np.int32))

probe_is_b = np.array([is_b_flavour(x) for x in probe_truth], dtype=bool)
probe_pass = (probe_wp == 1)

# -------------------------------------------------
# Compute probe-jet tagging efficiency by region
# -------------------------------------------------
print("\n" + "=" * 90)
print("leading_jet_pt_bin, numerator_sumw, denominator_sumw, efficiency, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)

    if not np.any(bin_mask):
        print(f"{label}, 0.000000, 0.000000, 0.000000, 0.000000")
        continue

    w_bin          = w_event[bin_mask]
    probe_is_b_bin = probe_is_b[bin_mask]
    probe_pass_bin = probe_pass[bin_mask]

    # Denominator: true b probe jets only
    den_weights = w_bin[probe_is_b_bin]

    # Numerator: true b probe jets passing the chosen WP
    num_weights = w_bin[probe_is_b_bin & probe_pass_bin]

    sumw_den  = np.sum(den_weights)
    sumw2_den = np.sum(den_weights ** 2)

    sumw_num  = np.sum(num_weights)
    sumw2_num = np.sum(num_weights ** 2)

    eff, err = efficiency_and_error(
        sumw_num=sumw_num,
        sumw2_num=sumw2_num,
        sumw_den=sumw_den,
        sumw2_den=sumw2_den,
    )

    print(f"{label}, {sumw_num:.6f}, {sumw_den:.6f}, {eff:.6f}, {err:.6f}")

print("=" * 90)

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------
print("\nDiagnostics:")
for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)
    count = np.count_nonzero(bin_mask)
    sumw  = np.sum(w_event[bin_mask]) if count > 0 else 0.0
    n_b_probe = np.count_nonzero(probe_is_b[bin_mask]) if count > 0 else 0
    print(f"{label}: events={count}, event_sumw={sumw:.6f}, true_b_probes={n_b_probe}")