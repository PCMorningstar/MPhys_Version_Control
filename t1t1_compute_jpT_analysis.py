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

# -------------------------------------------------
# Single-jet flavour mapping
# -------------------------------------------------
single_flavours = ["b", "c", "l"]

def single_flavour_label(f):
    f = abs(int(f))
    if f == 5:
        return "b"
    elif f == 4:
        return "c"
    elif f == 0:
        return "l"
    else:
        return None

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
    & (arr["jet_size_NOSYS"] == 3)
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

# -------------------------------------------------
# Extract flavours separately
# -------------------------------------------------
flav_i = np.array([single_flavour_label(f) for f in truth_i], dtype=object)
flav_j = np.array([single_flavour_label(f) for f in truth_j], dtype=object)

valid_i = np.isin(flav_i, single_flavours)
valid_j = np.isin(flav_j, single_flavours)

# apply masks separately
flav_i = flav_i[valid_i]
w_i = w_event[valid_i]

flav_j = flav_j[valid_j]
w_j = w_event[valid_j]

region_flags_i = {k: v[valid_i] for k, v in region_flags.items()}
region_flags_j = {k: v[valid_j] for k, v in region_flags.items()}

# -------------------------------------------------
# Output: Top1 (χ² index 0)
# -------------------------------------------------
print("\n" + "=" * 90)
print("Top1 (chi2 index 0): leading_jet_pt_bin, flavour, fraction, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags_i[label] == 1)

    if not np.any(bin_mask):
        for f in single_flavours:
            print(f"{label}, {f}, 0.000000, 0.000000")
        continue

    w_bin = w_i[bin_mask]
    flav_bin = flav_i[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin ** 2)

    for f in single_flavours:
        mask = (flav_bin == f)
        w_cat = w_bin[mask]

        sumw_cat = np.sum(w_cat)
        sumw2_cat = np.sum(w_cat ** 2)

        frac, err = weighted_fraction_and_error(
            sumw_cat=sumw_cat,
            sumw2_cat=sumw2_cat,
            sumw_tot=sumw_tot,
            sumw2_tot=sumw2_tot,
        )

        print(f"{label}, {f}, {frac:.6f}, {err:.6f}")

# -------------------------------------------------
# Output: Top2 (χ² index 1)
# -------------------------------------------------
print("\n" + "=" * 90)
print("Top2 (chi2 index 1): leading_jet_pt_bin, flavour, fraction, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags_j[label] == 1)

    if not np.any(bin_mask):
        for f in single_flavours:
            print(f"{label}, {f}, 0.000000, 0.000000")
        continue

    w_bin = w_j[bin_mask]
    flav_bin = flav_j[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin ** 2)

    for f in single_flavours:
        mask = (flav_bin == f)
        w_cat = w_bin[mask]

        sumw_cat = np.sum(w_cat)
        sumw2_cat = np.sum(w_cat ** 2)

        frac, err = weighted_fraction_and_error(
            sumw_cat=sumw_cat,
            sumw2_cat=sumw2_cat,
            sumw_tot=sumw_tot,
            sumw2_tot=sumw2_tot,
        )

        print(f"{label}, {f}, {frac:.6f}, {err:.6f}")

print("=" * 90)

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------
print("\nDiagnostics Top1:")
for label, _ in pt_regions:
    bin_mask = (region_flags_i[label] == 1)
    count = np.count_nonzero(bin_mask)
    sumw = np.sum(w_i[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_i[bin_mask] ** 2) if count > 0 else 0.0
    print(f"{label}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")

print("\nDiagnostics Top2:")
for label, _ in pt_regions:
    bin_mask = (region_flags_j[label] == 1)
    count = np.count_nonzero(bin_mask)
    sumw = np.sum(w_j[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_j[bin_mask] ** 2) if count > 0 else 0.0
    print(f"{label}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")