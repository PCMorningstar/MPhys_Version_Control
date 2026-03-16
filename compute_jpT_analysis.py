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
]

flavour_cats = ["bb", "bc", "bl", "cc", "cl", "ll"]

# -------------------------------------------------
# Helpers
# -------------------------------------------------
def flavour_label(f):
    f = abs(int(f))
    if f == 5:
        return "b"
    elif f == 4:
        return "c"
    elif f == 0:
        return "l"
    else:
        return None

def pair_category(f1, f2):
    l1 = flavour_label(f1)
    l2 = flavour_label(f2)
    if l1 is None or l2 is None:
        return None
    return "".join(sorted([l1, l2]))

def weighted_fraction_and_error(sumw_cat, sumw2_cat, sumw_tot, sumw2_tot):
    """
    Weighted fraction and uncertainty for an exclusive category:
        p = sumw_cat / sumw_tot

    Uses the standard variance expression for a weighted binomial-like fraction.
    """
    if sumw_tot <= 0.0:
        return 0.0, 0.0

    frac = sumw_cat / sumw_tot

    var = (
        ((1.0 - frac) ** 2) * sumw2_cat
        + (frac ** 2) * (sumw2_tot - sumw2_cat)
    ) / (sumw_tot ** 2)

    if var < 0.0:
        var = 0.0

    return frac, np.sqrt(var)

# -------------------------------------------------
# Load branches
# -------------------------------------------------
branches = [
    "selection_cuts_NOSYS",
    "jet_size_NOSYS",
    "ordered_jet_truth_flavour_NOSYS",
    "raw_chi2_minval_truthall_NOSYS",
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
    & (arr["jet_size_NOSYS"] >= 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi   = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]

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
# Extract χ2-selected pair
# -------------------------------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:, 0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:, 1], int))
n_jets = ak.to_numpy(ak.num(truth, axis=1))

valid = (
    (chi_i >= 0) & (chi_j >= 0)
    & (chi_i < n_jets) & (chi_j < n_jets)
    & (chi_i != chi_j)
)

truth = truth[valid]
w_event = ak.to_numpy(w_event[valid])
chi_i = chi_i[valid]
chi_j = chi_j[valid]

for label in region_flags:
    region_flags[label] = ak.to_numpy(region_flags[label][valid])

idx = np.arange(len(truth))
truth_i = ak.to_numpy(truth[idx, chi_i])
truth_j = ak.to_numpy(truth[idx, chi_j])

pair_cat = np.array(
    [pair_category(f1, f2) for f1, f2 in zip(truth_i, truth_j)],
    dtype=object
)

# keep only events that map to one of the six categories
valid_cat = np.isin(pair_cat, flavour_cats)

pair_cat = pair_cat[valid_cat]
w_event = w_event[valid_cat]

for label in region_flags:
    region_flags[label] = region_flags[label][valid_cat]

# -------------------------------------------------
# Output weighted event fractions by pT region
# -------------------------------------------------
print("\n" + "=" * 90)
print("leading_jet_pt_bin, flavour, fraction, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)

    if not np.any(bin_mask):
        for cat in flavour_cats:
            print(f"{label}, {cat}, 0.000000, 0.000000")
        continue

    w_bin = w_event[bin_mask]
    cats_bin = pair_cat[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin ** 2)

    for cat in flavour_cats:
        cat_mask = (cats_bin == cat)
        w_cat = w_bin[cat_mask]

        sumw_cat = np.sum(w_cat)
        sumw2_cat = np.sum(w_cat ** 2)

        frac, err = weighted_fraction_and_error(
            sumw_cat=sumw_cat,
            sumw2_cat=sumw2_cat,
            sumw_tot=sumw_tot,
            sumw2_tot=sumw2_tot,
        )

        print(f"{label}, {cat}, {frac:.6f}, {err:.6f}")

print("=" * 90)

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------
print("\nDiagnostics:")
for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)
    count = np.count_nonzero(bin_mask)
    sumw = np.sum(w_event[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_event[bin_mask] ** 2) if count > 0 else 0.0
    print(f"{label}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")