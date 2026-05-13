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
# Stored pT region flags from ntuple
# -------------------------------------------------
pt_regions = [
    ("0to30", "jet_pt_region_0to30_GeV_NOSYS"),
    ("30to60", "jet_pt_region_30to60_GeV_NOSYS"),
    ("60to90", "jet_pt_region_60to90_GeV_NOSYS"),
    ("90to120", "jet_pt_region_90to120_GeV_NOSYS"),
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
    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
    "jet_pt_new_NOSYS",
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
jet_pt = arr["jet_pt_new_NOSYS"][base_mask]

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
# Require exactly 2 truth flavours and 2 jet pT entries
# -------------------------------------------------
n_truth_flavours = ak.to_numpy(ak.num(truth, axis=1))
n_jet_pt = ak.to_numpy(ak.num(jet_pt, axis=1))

valid = (
    (n_truth_flavours == 2)
    & (n_jet_pt == 2)
)

truth = truth[valid]
jet_pt = jet_pt[valid]
w_event = ak.to_numpy(w_event[valid])

for label in region_flags:
    region_flags[label] = ak.to_numpy(region_flags[label][valid])

# -------------------------------------------------
# Separate jets into leading and subleading by pT
# -------------------------------------------------
idx = np.arange(len(truth))

pt0 = ak.to_numpy(jet_pt[:, 0])
pt1 = ak.to_numpy(jet_pt[:, 1])

lead_idx = np.where(pt0 >= pt1, 0, 1)
sublead_idx = np.where(pt0 >= pt1, 1, 0)

truth_lead = ak.to_numpy(truth[idx, lead_idx])
truth_sublead = ak.to_numpy(truth[idx, sublead_idx])

# -------------------------------------------------
# Extract flavours separately for leading and subleading jets
# -------------------------------------------------
flav_lead = np.array([single_flavour_label(f) for f in truth_lead], dtype=object)
flav_sublead = np.array([single_flavour_label(f) for f in truth_sublead], dtype=object)

valid_lead = np.isin(flav_lead, single_flavours)
valid_sublead = np.isin(flav_sublead, single_flavours)

flav_lead = flav_lead[valid_lead]
w_lead = w_event[valid_lead]

flav_sublead = flav_sublead[valid_sublead]
w_sublead = w_event[valid_sublead]

region_flags_lead = {k: v[valid_lead] for k, v in region_flags.items()}
region_flags_sublead = {k: v[valid_sublead] for k, v in region_flags.items()}

# -------------------------------------------------
# Lists to store results
# -------------------------------------------------
leading_fraction_list = []
leading_yield_list = []

subleading_fraction_list = []
subleading_yield_list = []

# -------------------------------------------------
# Leading jet
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = region_flags_lead[label] == 1

    if not np.any(bin_mask):
        for f in single_flavours:
            leading_fraction_list.append((label, f, 0.0, 0.0))
            leading_yield_list.append((label, f, 0.0, 0.0))
        continue

    w_bin = w_lead[bin_mask]
    flav_bin = flav_lead[bin_mask]

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

        leading_fraction_list.append((label, f, frac, err_frac))
        leading_yield_list.append((label, f, yield_val, err_yield))

# -------------------------------------------------
# Subleading jet
# -------------------------------------------------
for label, _ in pt_regions:
    bin_mask = region_flags_sublead[label] == 1

    if not np.any(bin_mask):
        for f in single_flavours:
            subleading_fraction_list.append((label, f, 0.0, 0.0))
            subleading_yield_list.append((label, f, 0.0, 0.0))
        continue

    w_bin = w_sublead[bin_mask]
    flav_bin = flav_sublead[bin_mask]

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

        subleading_fraction_list.append((label, f, frac, err_frac))
        subleading_yield_list.append((label, f, yield_val, err_yield))

# -------------------------------------------------
# Print results
# -------------------------------------------------

print("\n" + "=" * 110)
print("Leading jet: weighted event FRACTIONS")
print("Format: leading_jet_pt_bin, flavour, fraction, error")
print("=" * 110)
for row in leading_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Leading jet: weighted NUMBER OF EVENTS")
print("Format: leading_jet_pt_bin, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in leading_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Subleading jet: weighted event FRACTIONS")
print("Format: subleading_jet_pt_bin, flavour, fraction, error")
print("=" * 110)
for row in subleading_fraction_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("\n" + "=" * 110)
print("Subleading jet: weighted NUMBER OF EVENTS")
print("Format: subleading_jet_pt_bin, flavour, yield(sumw), error(sqrt(sumw2))")
print("=" * 110)
for row in subleading_yield_list:
    print(f"{row[0]}, {row[1]}, {row[2]:.6f}, {row[3]:.6f}")

print("=" * 110)