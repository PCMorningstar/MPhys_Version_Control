import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# Stored leading-jet pT region flags from ntuple
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
    ("300to330", "jet_pt_region_300to330_GeV_NOSYS"),
    ("330to360", "jet_pt_region_330to360_GeV_NOSYS"),
    ("360to390", "jet_pt_region_360to390_GeV_NOSYS"),
    ("390to420", "jet_pt_region_390to420_GeV_NOSYS"),
    ("420to450", "jet_pt_region_420to450_GeV_NOSYS"),
    ("450to480", "jet_pt_region_450to480_GeV_NOSYS"),
    ("480to510", "jet_pt_region_480to510_GeV_NOSYS"),
    ("510to540", "jet_pt_region_510to540_GeV_NOSYS"),
    ("540to570", "jet_pt_region_540to570_GeV_NOSYS"),
    ("570to600", "jet_pt_region_570to600_GeV_NOSYS"),

]

# -------------------------------------------------
# Helpers
# -------------------------------------------------
def is_b_flavour(f):
    """Return True if jet truth flavour corresponds to a b-jet."""
    try:
        return abs(int(f)) == 5
    except Exception:
        return False

def efficiency_and_error(sumw_num, sumw2_num, sumw_den, sumw2_den):
    """
    Weighted efficiency and approximate binomial-like uncertainty:
        eff = sumw_num / sumw_den

    Variance expression for weighted pass/fail counting.
    """
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
    "jet_select_GN2v01_FixedCutBEff_77_NOSYS",
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
wp77  = arr["jet_select_GN2v01_FixedCutBEff_77_NOSYS"][base_mask]

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
# Assumes raw_chi2_minval_truthall_NOSYS[:,0:2] are the selected jet indices
# -------------------------------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:, 0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:, 1], int))
n_jets = ak.to_numpy(ak.num(truth, axis=1))

valid = (
    (chi_i >= 0) & (chi_j >= 0)
    & (chi_i < n_jets) & (chi_j < n_jets)
    & (chi_i != chi_j)
)

truth   = truth[valid]
wp77    = wp77[valid]
w_event = ak.to_numpy(w_event[valid])
chi_i   = chi_i[valid]
chi_j   = chi_j[valid]

for label in region_flags:
    region_flags[label] = ak.to_numpy(region_flags[label][valid])

idx = np.arange(len(truth))

truth_i = ak.to_numpy(truth[idx, chi_i])
truth_j = ak.to_numpy(truth[idx, chi_j])

wp_i = ak.to_numpy(ak.values_astype(wp77[idx, chi_i], np.int32))
wp_j = ak.to_numpy(ak.values_astype(wp77[idx, chi_j], np.int32))

# -------------------------------------------------
# Compute inclusive tagging efficiency by leading-jet pT region
# Each selected true b-jet contributes separately
# -------------------------------------------------
print("\n" + "=" * 90)
print("leading_jet_pt_bin, numerator_sumw, denominator_sumw, efficiency, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)

    if not np.any(bin_mask):
        print(f"{label}, 0.000000, 0.000000, 0.000000, 0.000000")
        continue

    # Event quantities in this pT region
    w_bin       = w_event[bin_mask]
    truth_i_bin = truth_i[bin_mask]
    truth_j_bin = truth_j[bin_mask]
    wp_i_bin    = wp_i[bin_mask]
    wp_j_bin    = wp_j[bin_mask]

    # Per-selected-jet truth b flags
    jet1_is_b = np.array([is_b_flavour(x) for x in truth_i_bin], dtype=bool)
    jet2_is_b = np.array([is_b_flavour(x) for x in truth_j_bin], dtype=bool)

    # Per-selected-jet pass flags at WP77
    jet1_pass = (wp_i_bin == 1)
    jet2_pass = (wp_j_bin == 1)

    # Denominator: all selected true b-jets
    den_weights = np.concatenate([
        w_bin[jet1_is_b],
        w_bin[jet2_is_b],
    ])

    # Numerator: all selected true b-jets passing WP77
    num_weights = np.concatenate([
        w_bin[jet1_is_b & jet1_pass],
        w_bin[jet2_is_b & jet2_pass],
    ])

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
    sumw = np.sum(w_event[bin_mask]) if count > 0 else 0.0
    print(f"{label}: events={count}, event_sumw={sumw:.6f}")