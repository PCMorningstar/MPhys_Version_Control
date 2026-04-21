import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# YAML-defined SV invariant mass regions
# -------------------------------------------------
sv_regions = [
    ("sv_invariant_mass_region_neg0point5upto0_GeV_region",   "sv_invariant_mass_region_neg0point5upto0_GeV_NOSYS"),
    ("sv_invariant_mass_region_0to0point5_GeV_region",   "sv_invariant_mass_region_0to0point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_0point5to1_GeV_region",   "sv_invariant_mass_region_0point5to1_GeV_NOSYS"),
    ("sv_invariant_mass_region_1to1point5_GeV_region",   "sv_invariant_mass_region_1to1point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_1point5to2_GeV_region",   "sv_invariant_mass_region_1point5to2_GeV_NOSYS"),
    ("sv_invariant_mass_region_2to2point5_GeV_region",   "sv_invariant_mass_region_2to2point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_2point5to3_GeV_region",   "sv_invariant_mass_region_2point5to3_GeV_NOSYS"),
    ("sv_invariant_mass_region_3to3point5_GeV_region",   "sv_invariant_mass_region_3to3point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_3point5to4_GeV_region",   "sv_invariant_mass_region_3point5to4_GeV_NOSYS"),
    ("sv_invariant_mass_region_4to4point5_GeV_region",   "sv_invariant_mass_region_4to4point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_4point5to5_GeV_region",   "sv_invariant_mass_region_4point5to5_GeV_NOSYS"),
    ("sv_invariant_mass_region_5to5point5_GeV_region",   "sv_invariant_mass_region_5to5point5_GeV_NOSYS"),
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

    Variance used:
        var(p) = [ (1-p)^2 * sumw2_cat + p^2 * (sumw2_tot - sumw2_cat) ] / (sumw_tot^2)
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

def weighted_yield_and_error(w):
    """
    Weighted event yield and uncertainty:
        N = sum w
        sigma_N = sqrt(sum w^2)
    """
    sumw = np.sum(w)
    sumw2 = np.sum(w ** 2)
    err = np.sqrt(sumw2)
    return sumw, err

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
] + [branch for _, branch in sv_regions]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

# -------------------------------------------------
# Base event selection: match YAML baseline for SV regions
# selection_cuts_NOSYS == 1 && jet_size_NOSYS >= 2
# -------------------------------------------------
base_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] >= 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi   = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]

w_event = ak.to_numpy((
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[base_mask])

# Event-level SV region flags from YAML
sv_region_flags = {
    region_name: ak.to_numpy(arr[branch_name][base_mask])
    for region_name, branch_name in sv_regions
}

# -------------------------------------------------
# Extract chi2-selected pair
# Assumes chi[:,0] and chi[:,1] are the selected jet indices
# -------------------------------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:, 0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:, 1], int))
n_jets_truth = ak.to_numpy(ak.num(truth, axis=1))

valid = (
    (chi_i >= 0) & (chi_j >= 0)
    & (chi_i < n_jets_truth) & (chi_j < n_jets_truth)
    & (chi_i != chi_j)
)

truth = truth[valid]
w_event = w_event[valid]
chi_i = chi_i[valid]
chi_j = chi_j[valid]

for key in sv_region_flags:
    sv_region_flags[key] = sv_region_flags[key][valid]

idx = np.arange(len(truth))
truth_i = ak.to_numpy(truth[idx, chi_i])
truth_j = ak.to_numpy(truth[idx, chi_j])

# -------------------------------------------------
# Extract flavours separately for Top1 and Top2
# -------------------------------------------------
flav_i = np.array([single_flavour_label(f) for f in truth_i], dtype=object)
flav_j = np.array([single_flavour_label(f) for f in truth_j], dtype=object)

valid_i = np.isin(flav_i, single_flavours)
valid_j = np.isin(flav_j, single_flavours)

flav_i = flav_i[valid_i]
w_i = w_event[valid_i]

flav_j = flav_j[valid_j]
w_j = w_event[valid_j]

sv_region_flags_i = {k: v[valid_i] for k, v in sv_region_flags.items()}
sv_region_flags_j = {k: v[valid_j] for k, v in sv_region_flags.items()}

# -------------------------------------------------
# Lists to store results
# -------------------------------------------------
top1_fraction_list = []
top1_yield_list = []

top2_fraction_list = []
top2_yield_list = []

# -------------------------------------------------
# Top1 (chi2 index 0), using YAML SV invariant-mass regions
# -------------------------------------------------
for region_name, _ in sv_regions:
    bin_mask = (sv_region_flags_i[region_name] == 1)

    if not np.any(bin_mask):
        for f in single_flavours:
            top1_fraction_list.append((region_name, f, 0.0, 0.0))
            top1_yield_list.append((region_name, f, 0.0, 0.0))
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

        frac, err_frac = weighted_fraction_and_error(
            sumw_cat=sumw_cat,
            sumw2_cat=sumw2_cat,
            sumw_tot=sumw_tot,
            sumw2_tot=sumw2_tot,
        )

        yield_val, err_yield = weighted_yield_and_error(w_cat)

        top1_fraction_list.append((region_name, f, frac, err_frac))
        top1_yield_list.append((region_name, f, yield_val, err_yield))

# -------------------------------------------------
# Top2 (chi2 index 1), using YAML SV invariant-mass regions
# -------------------------------------------------
for region_name, _ in sv_regions:
    bin_mask = (sv_region_flags_j[region_name] == 1)

    if not np.any(bin_mask):
        for f in single_flavours:
            top2_fraction_list.append((region_name, f, 0.0, 0.0))
            top2_yield_list.append((region_name, f, 0.0, 0.0))
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

        frac, err_frac = weighted_fraction_and_error(
            sumw_cat=sumw_cat,
            sumw2_cat=sumw2_cat,
            sumw_tot=sumw_tot,
            sumw2_tot=sumw2_tot,
        )

        yield_val, err_yield = weighted_yield_and_error(w_cat)

        top2_fraction_list.append((region_name, f, frac, err_frac))
        top2_yield_list.append((region_name, f, yield_val, err_yield))

# -------------------------------------------------
# Print results
# -------------------------------------------------
print("\n" + "=" * 110)
print("ERROR FORMULAE USED")
print("=" * 110)
print("Weighted fraction p = sumw_cat / sumw_tot")
print("Variance on fraction:")
print("  var(p) = [ (1-p)^2 * sumw2_cat + p^2 * (sumw2_tot - sumw2_cat) ] / (sumw_tot^2)")
print("  sigma_p = sqrt(var(p))")
print("")
print("Weighted event yield N = sumw")
print("Uncertainty on yield:")
print("  sigma_N = sqrt(sumw2) = sqrt(sum(w^2))")
print("=" * 110)

print("\n" + "=" * 110)
print("Top1 (chi2 index 0): weighted event FRACTIONS vs YAML SV invariant-mass regions")
print("Format: region_name, flavour, fraction, error")
print("=" * 110)
for row in top1_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top1 (chi2 index 0): weighted NUMBER OF EVENTS vs YAML SV invariant-mass regions")
print("Format: region_name, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in top1_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top2 (chi2 index 1): weighted event FRACTIONS vs YAML SV invariant-mass regions")
print("Format: region_name, flavour, fraction, error")
print("=" * 110)
for row in top2_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top2 (chi2 index 1): weighted NUMBER OF EVENTS vs YAML SV invariant-mass regions")
print("Format: region_name, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in top2_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("=" * 110)

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------
print("\nDiagnostics Top1:")
for region_name, _ in sv_regions:
    bin_mask = (sv_region_flags_i[region_name] == 1)
    count = np.count_nonzero(bin_mask)
    sumw = np.sum(w_i[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_i[bin_mask] ** 2) if count > 0 else 0.0
    print(f"{region_name}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")

print("\nDiagnostics Top2:")
for region_name, _ in sv_regions:
    bin_mask = (sv_region_flags_j[region_name] == 1)
    count = np.count_nonzero(bin_mask)
    sumw = np.sum(w_j[bin_mask]) if count > 0 else 0.0
    sumw2 = np.sum(w_j[bin_mask] ** 2) if count > 0 else 0.0
    print(f"{region_name}: count={count}, sumw={sumw:.6f}, sumw2={sumw2:.6f}")

# -------------------------------------------------
# Optional: lists only, for easy storage/copying
# -------------------------------------------------
print("\n" + "=" * 110)
print("PYTHON LISTS FOR EASY STORAGE")
print("=" * 110)

print("\ntop1_fraction_list = [")
for row in top1_fraction_list:
    print(f'    ("{row[0]}", "{row[1]}", {row[2]:.10f}, {row[3]:.10f}),')
print("]")

print("\ntop1_yield_list = [")
for row in top1_yield_list:
    print(f'    ("{row[0]}", "{row[1]}", {row[2]:.10f}, {row[3]:.10f}),')
print("]")

print("\ntop2_fraction_list = [")
for row in top2_fraction_list:
    print(f'    ("{row[0]}", "{row[1]}", {row[2]:.10f}, {row[3]:.10f}),')
print("]")

print("\ntop2_yield_list = [")
for row in top2_yield_list:
    print(f'    ("{row[0]}", "{row[1]}", {row[2]:.10f}, {row[3]:.10f}),')
print("]")