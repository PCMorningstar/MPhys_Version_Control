import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree = "reco"

# -------------------------------------------------
# MC normalisation
# -------------------------------------------------
luminosity = 29300.0
xsec = 85.482
filter_eff = 1.0
kfactor = 1.138433852
sum_of_weights = 4268786417.0

norm = luminosity * xsec * filter_eff * kfactor / sum_of_weights

# -------------------------------------------------
# Stored SV invariant mass region flags from ntuple
# -------------------------------------------------
pt_regions = [
    ("neg0point5upto0", "sv_invariant_mass_region_neg0point5upto0_GeV_NOSYS"),
    ("0to0point5", "sv_invariant_mass_region_0to0point5_GeV_NOSYS"),
    ("0point5to1", "sv_invariant_mass_region_0point5to1_GeV_NOSYS"),
    ("1to1point5", "sv_invariant_mass_region_1to1point5_GeV_NOSYS"),
    ("1point5to2", "sv_invariant_mass_region_1point5to2_GeV_NOSYS"),
    ("2to2point5", "sv_invariant_mass_region_2to2point5_GeV_NOSYS"),
    ("2point5to3", "sv_invariant_mass_region_2point5to3_GeV_NOSYS"),
    ("3to3point5", "sv_invariant_mass_region_3to3point5_GeV_NOSYS"),
    ("3point5to4", "sv_invariant_mass_region_3point5to4_GeV_NOSYS"),
    ("4to4point5", "sv_invariant_mass_region_4to4point5_GeV_NOSYS"),
    ("4point5to5", "sv_invariant_mass_region_4point5to5_GeV_NOSYS"),
    ("5to5point5", "sv_invariant_mass_region_5to5point5_GeV_NOSYS"),
]

# -------------------------------------------------
# Single-jet flavour mapping
# -------------------------------------------------
single_flavours = ["b", "c", "l"]

def single_flavour_label(f):
    f = abs(int(f))
    if f == 5:
        return "b"
    if f == 4:
        return "c"
    if f == 0:
        return "l"
    return None

def weighted_fraction_and_error(sumw_cat, sumw2_cat, sumw_tot, sumw2_tot):
    if sumw_tot <= 0.0:
        return 0.0, 0.0

    frac = sumw_cat / sumw_tot

    var = (
        ((1.0 - frac) ** 2) * sumw2_cat
        + (frac ** 2) * (sumw2_tot - sumw2_cat)
    ) / (sumw_tot ** 2)

    return frac, np.sqrt(max(var, 0.0))

def weighted_yield_and_error(w):
    sumw = np.sum(w)
    sumw2 = np.sum(w ** 2)
    return sumw, np.sqrt(sumw2)

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
# Base event selection: exactly 2 jets
# -------------------------------------------------
base_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] >= 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]

# -------------------------------------------------
# Total FastFrames-like event weight
# -------------------------------------------------
w_event = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
    * norm
)[base_mask]

region_flags = {
    label: arr[branch][base_mask]
    for label, branch in pt_regions
}

# -------------------------------------------------
# Extract chi2-selected pair
# Assumes chi[:,0] and chi[:,1] are selected jet indices
# -------------------------------------------------
chi_i = ak.to_numpy(ak.values_astype(chi[:, 0], int))
chi_j = ak.to_numpy(ak.values_astype(chi[:, 1], int))
n_jets = ak.to_numpy(ak.num(truth, axis=1))

valid = (
    (chi_i >= 0)
    & (chi_j >= 0)
    & (chi_i < n_jets)
    & (chi_j < n_jets)
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

region_flags_i = {k: v[valid_i] for k, v in region_flags.items()}
region_flags_j = {k: v[valid_j] for k, v in region_flags.items()}

# -------------------------------------------------
# Lists to store results
# -------------------------------------------------
top1_fraction_list = []
top1_yield_list = []

top2_fraction_list = []
top2_yield_list = []

# -------------------------------------------------
# Top1 (chi2 index 0)
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = region_flags_i[label] == 1

    if not np.any(bin_mask):
        for f in single_flavours:
            top1_fraction_list.append((label, f, 0.0, 0.0))
            top1_yield_list.append((label, f, 0.0, 0.0))
        continue

    w_bin = w_i[bin_mask]
    flav_bin = flav_i[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin ** 2)

    for f in single_flavours:
        mask = flav_bin == f
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

        top1_fraction_list.append((label, f, frac, err_frac))
        top1_yield_list.append((label, f, yield_val, err_yield))

# -------------------------------------------------
# Top2 (chi2 index 1)
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = region_flags_j[label] == 1

    if not np.any(bin_mask):
        for f in single_flavours:
            top2_fraction_list.append((label, f, 0.0, 0.0))
            top2_yield_list.append((label, f, 0.0, 0.0))
        continue

    w_bin = w_j[bin_mask]
    flav_bin = flav_j[bin_mask]

    sumw_tot = np.sum(w_bin)
    sumw2_tot = np.sum(w_bin ** 2)

    for f in single_flavours:
        mask = flav_bin == f
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

        top2_fraction_list.append((label, f, frac, err_frac))
        top2_yield_list.append((label, f, yield_val, err_yield))



print("\n" + "=" * 110)
print("Top1 (chi2 index 0): weighted event FRACTIONS")
print("Format: sv_invariant_mass_region, flavour, fraction, error")
print("=" * 110)
for row in top1_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top1 (chi2 index 0): weighted NUMBER OF EVENTS")
print("Format: sv_invariant_mass_region, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in top1_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top2 (chi2 index 1): weighted event FRACTIONS")
print("Format: sv_invariant_mass_region, flavour, fraction, error")
print("=" * 110)
for row in top2_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Top2 (chi2 index 1): weighted NUMBER OF EVENTS")
print("Format: sv_invariant_mass_region, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in top2_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("=" * 110)