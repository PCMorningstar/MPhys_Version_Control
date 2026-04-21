import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# Stored leading-jet pT region flags from ntuple
# These are event-level flags based on the leading jet pT
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

def to_event_scalar(x, default=0):
    """
    Convert a branch to one scalar per event.
    If already flat, return as numpy.
    If jagged, take first entry per event and fill missing with default.
    """
    t = str(ak.type(x))
    if "var *" not in t:
        return ak.to_numpy(x)
    return ak.to_numpy(ak.fill_none(ak.firsts(x), default))

def to_event_weight(x, default=1.0):
    """
    Convert a weight branch to one scalar per event.
    If flat, return as awkward/numpy-compatible.
    If jagged, multiply entries per event.
    """
    t = str(ak.type(x))
    if "var *" not in t:
        return x
    return ak.fill_none(ak.prod(x, axis=1), default)

# -------------------------------------------------
# Load branches
# -------------------------------------------------
branches = [
    "selection_cuts_NOSYS",
    "jet_size_NOSYS",
    "ordered_jet_truth_flavour_NOSYS",
    "jet_pt_new_NOSYS",
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
# Event-level mask
# -------------------------------------------------
selection_cuts = to_event_scalar(arr["selection_cuts_NOSYS"], default=0)
jet_size       = to_event_scalar(arr["jet_size_NOSYS"], default=-1)

base_mask = (
    (selection_cuts == 1)
    & (jet_size == 2)
)

# -------------------------------------------------
# Object branches after event selection
# -------------------------------------------------
truth  = arr["ordered_jet_truth_flavour_NOSYS"][base_mask]
jet_pt = arr["jet_pt_new_NOSYS"][base_mask]

# Choose the WP branch you actually want here
wp_probe = arr["jet_select_GN2v01_FixedCutBEff_90_NOSYS"][base_mask]

# -------------------------------------------------
# Event weights
# -------------------------------------------------
w_mc     = to_event_weight(arr["weight_mc_NOSYS"])[base_mask]
w_pileup = to_event_weight(arr["weight_pileup_NOSYS"])[base_mask]
w_lepSF  = to_event_weight(arr["weight_leptonSF_tight_NOSYS"])[base_mask]
w_jvtSF  = to_event_weight(arr["weight_jvt_effSF_NOSYS"])[base_mask]

w_event = ak.to_numpy(w_mc * w_pileup * w_lepSF * w_jvtSF)

# -------------------------------------------------
# Region flags after event selection
# -------------------------------------------------
region_flags = {
    label: to_event_scalar(arr[branch], default=0)[base_mask]
    for label, branch in pt_regions
}

# -------------------------------------------------
# Require exactly 2 jets in the selected branches as well
# -------------------------------------------------
n_truth = ak.to_numpy(ak.num(truth, axis=1))
n_pt    = ak.to_numpy(ak.num(jet_pt, axis=1))
n_wp    = ak.to_numpy(ak.num(wp_probe, axis=1))

valid = (
    (n_truth == 2)
    & (n_pt == 2)
    & (n_wp == 2)
)

truth    = truth[valid]
jet_pt   = jet_pt[valid]
wp_probe = wp_probe[valid]
w_event  = w_event[valid]

for label in region_flags:
    region_flags[label] = region_flags[label][valid]

# -------------------------------------------------
# Separate jets into leading and subleading by pT
# -------------------------------------------------
idx = np.arange(len(truth))

pt0 = ak.to_numpy(jet_pt[:, 0])
pt1 = ak.to_numpy(jet_pt[:, 1])

lead_idx    = np.where(pt0 >= pt1, 0, 1)
sublead_idx = np.where(pt0 >= pt1, 1, 0)

# Leading jet = probe jet
probe_truth = ak.to_numpy(truth[idx, sublead_idx])
probe_wp    = ak.to_numpy(ak.values_astype(wp_probe[idx, sublead_idx], np.int32))
probe_pt    = ak.to_numpy(jet_pt[idx, sublead_idx])

# Optional diagnostics
sublead_truth = ak.to_numpy(truth[idx, lead_idx])
sublead_wp    = ak.to_numpy(ak.values_astype(wp_probe[idx, lead_idx], np.int32))
sublead_pt    = ak.to_numpy(jet_pt[idx, lead_idx])

probe_is_b = np.array([is_b_flavour(x) for x in probe_truth], dtype=bool)
probe_pass = (probe_wp == 1)

# -------------------------------------------------
# Compute leading-jet probe tagging efficiency by region
# -------------------------------------------------
print("\n" + "=" * 90)
print("subleading_jet_pt_bin, numerator_sumw, denominator_sumw, efficiency, error")
print("=" * 90)

for label, _ in pt_regions:
    bin_mask = (region_flags[label] == 1)

    if not np.any(bin_mask):
        print(f"{label}, 0.000000, 0.000000, 0.000000, 0.000000")
        continue

    w_bin          = w_event[bin_mask]
    probe_is_b_bin = probe_is_b[bin_mask]
    probe_pass_bin = probe_pass[bin_mask]

    den_weights = w_bin[probe_is_b_bin]
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
