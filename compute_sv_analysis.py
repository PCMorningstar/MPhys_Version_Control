import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# Stored SV invariant mass region flags
# -------------------------------------------------
sv_regions = [
    ("0to0point5",   "sv_invariant_mass_region_0to0point5_GeV_NOSYS"),
    ("0point5to1",   "sv_invariant_mass_region_0point5to1_GeV_NOSYS"),
    ("1to1point5",   "sv_invariant_mass_region_1to1point5_GeV_NOSYS"),
    ("1point5to2",   "sv_invariant_mass_region_1point5to2_GeV_NOSYS"),
    ("2to2point5",   "sv_invariant_mass_region_2to2point5_GeV_NOSYS"),
    ("2point5to3",   "sv_invariant_mass_region_2point5to3_GeV_NOSYS"),
    ("3to3point5",   "sv_invariant_mass_region_3to3point5_GeV_NOSYS"),
    ("3point5to4",   "sv_invariant_mass_region_3point5to4_GeV_NOSYS"),
    ("4to4point5",   "sv_invariant_mass_region_4to4point5_GeV_NOSYS"),
    ("4point5to5",   "sv_invariant_mass_region_4point5to5_GeV_NOSYS"),
]

flavour_cats = ["bb", "bc", "bl", "cc", "cl", "ll"]

# --------------------------------------------------
# Helper
# --------------------------------------------------
def flavour_label(f, b_code=5, c_code=4, l_code=0):
    if abs(f) == b_code:
        return "b"
    elif abs(f) == c_code:
        return "c"
    elif f == l_code:
        return "l"
    else:
        return None

def pair_category(f1, f2):
    a = flavour_label(f1)
    b = flavour_label(f2)

    if a is None or b is None:
        return None

    return "".join(sorted([a, b]))

def weighted_fraction_and_error(sumw_cat, sumw2_cat, sumw_tot, sumw2_tot):

    if sumw_tot <= 0.0:
        return 0.0, 0.0

    frac = sumw_cat / sumw_tot

    var = (
        ((1.0 - frac) ** 2) * sumw2_cat +
        (frac ** 2) * (sumw2_tot - sumw2_cat)
    ) / (sumw_tot ** 2)

    var = max(var, 0.0)

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
] + [b for _, b in sv_regions]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

# -------------------------------------------------
# Base event selection
# -------------------------------------------------
base_mask = (
    (arr["selection_cuts_NOSYS"] == 1) &
    (arr["jet_size_NOSYS"] >= 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi   = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]

w_event = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[base_mask]

# region flags
region_flags = {
    label: arr[branch][base_mask]
    for label, branch in sv_regions
}

# -------------------------------------------------
# Extract χ2 pair
# -------------------------------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:,0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:,1], int))

n_jets = ak.to_numpy(ak.num(truth, axis=1))
w_event = ak.to_numpy(w_event)

valid = (
    np.isfinite(w_event)
    & (chi_i >= 0)
    & (chi_j >= 0)
    & (chi_i < n_jets)
    & (chi_j < n_jets)
    & (chi_i != chi_j)
)

truth = truth[valid]
w_event = w_event[valid]

chi_i = chi_i[valid]
chi_j = chi_j[valid]

for label in region_flags:
    region_flags[label] = ak.to_numpy(region_flags[label][valid])

idx = np.arange(len(truth))

truth_i = ak.to_numpy(truth[idx, chi_i])
truth_j = ak.to_numpy(truth[idx, chi_j])

pair_cat = np.array([
    pair_category(f1, f2)
    for f1, f2 in zip(truth_i, truth_j)
])

# -------------------------------------------------
# Output weighted event fractions by SV mass region
# -------------------------------------------------
print("\n" + "="*90)
print("sv_invariant_mass_region, flavour, fraction, error")
print("="*90)

for label, _ in sv_regions:

    bin_mask = (region_flags[label] == 1)

    if not np.any(bin_mask):
        for cat in flavour_cats:
            print(f"{label}, {cat}, 0.000000, 0.000000")
        continue

    w_bin = w_event[bin_mask]
    cats_bin = pair_cat[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin**2)

    for cat in flavour_cats:

        cat_mask = (cats_bin == cat)

        w_cat = w_bin[cat_mask]

        sumw_cat = np.sum(w_cat)
        sumw2_cat = np.sum(w_cat**2)

        frac, err = weighted_fraction_and_error(
            sumw_cat,
            sumw2_cat,
            sumw_tot,
            sumw2_tot
        )

        print(f"{label}, {cat}, {frac:.6f}, {err:.6f}")

print("="*90)

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------
print("\nDiagnostics:")

for label, _ in sv_regions:

    bin_mask = (region_flags[label] == 1)

    count = np.count_nonzero(bin_mask)
    sumw  = np.sum(w_event[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_event[bin_mask]**2) if count > 0 else 0.0

    print(f"{label}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")