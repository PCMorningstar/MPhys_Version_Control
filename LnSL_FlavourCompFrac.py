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
    ("300to360", "jet_pt_region_300to360_GeV_NOSYS"),
    ("360to900", "jet_pt_region_360to900_GeV_NOSYS"),
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

    This is the same weighted binomial-like expression you were already using.
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

    This is the standard weighted counting uncertainty.
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
    "jet_pt_new_NOSYS",
] + [b for _, b in pt_regions]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

# -------------------------------------------------
# Base event selection: NOW 2 JETS
# -------------------------------------------------
base_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] == 2)
)

truth = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
chi   = arr["raw_chi2_minval_truthall_NOSYS"][base_mask]
jet_pt = arr["jet_pt_new_NOSYS"][base_mask]

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
# Extract chi2-selected pair - leading and subleading jets
# -------------------------------------------------

n_truth_flavours = ak.to_numpy(ak.num(truth, axis=1))
n_jet_pt         = ak.to_numpy(ak.num(jet_pt, axis=1))

valid = (
    (n_truth_flavours == 2)
    & (n_jet_pt == 2)
)

truth = truth[valid]
jet_pt = jet_pt[valid]
w_event = w_event[valid]

# -------------------------------------------------
# Separate jets into leading and subleading by pT
# -------------------------------------------------
idx = np.arange(len(truth))

pt0 = ak.to_numpy(jet_pt[:, 0])
pt1 = ak.to_numpy(jet_pt[:, 1])

lead_idx    = np.where(pt0 >= pt1, 0, 1)
sublead_idx = np.where(pt0 >= pt1, 1, 0)

truth_i = ak.to_numpy(truth[idx, lead_idx])
truth_j = ak.to_numpy(truth[idx, sublead_idx])

# -------------------------------------------------
# Extract flavours separately for leading and subleading jets
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
# Lists to store results
# -------------------------------------------------
top1_fraction_list = []
top1_yield_list = []

top2_fraction_list = []
top2_yield_list = []

# -------------------------------------------------
# Leading jet
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = (region_flags_i[label] == 1)

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

        top1_fraction_list.append((label, f, frac, err_frac))
        top1_yield_list.append((label, f, yield_val, err_yield))

# -------------------------------------------------
# Subleading jet
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = (region_flags_j[label] == 1)

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

        top2_fraction_list.append((label, f, frac, err_frac))
        top2_yield_list.append((label, f, yield_val, err_yield))

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

# -------------------------------------------------
# Leading jet fractions
# -------------------------------------------------
print("\n" + "=" * 110)
print("Leading jet: weighted event FRACTIONS")
print("Format: leading_jet_pt_bin, flavour, fraction, error")
print("=" * 110)

for row in top1_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

# -------------------------------------------------
# Leading jet yields
# -------------------------------------------------
print("\n" + "=" * 110)
print("Leading jet: weighted NUMBER OF EVENTS")
print("Format: leading_jet_pt_bin, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)

for row in top1_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

# -------------------------------------------------
# Subleading jet fractions
# -------------------------------------------------
print("\n" + "=" * 110)
print("Subleading jet: weighted event FRACTIONS")
print("Format: leading_jet_pt_bin, flavour, fraction, error")
print("=" * 110)

for row in top2_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

# -------------------------------------------------
# Subleading jet yields
# -------------------------------------------------
print("\n" + "=" * 110)
print("Subleading jet: weighted NUMBER OF EVENTS")
print("Format: leading_jet_pt_bin, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)

for row in top2_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("=" * 110)
