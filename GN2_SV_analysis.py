import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree  = "reco"

# -------------------------------------------------
# YAML-defined SV invariant mass regions
# -------------------------------------------------
sv_regions = [
    ("sv_invariant_mass_region_neg0point5upto0_GeV_region", "sv_invariant_mass_region_neg0point5upto0_GeV_NOSYS"),
    ("sv_invariant_mass_region_0to0point5_GeV_region",      "sv_invariant_mass_region_0to0point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_0point5to1_GeV_region",      "sv_invariant_mass_region_0point5to1_GeV_NOSYS"),
    ("sv_invariant_mass_region_1to1point5_GeV_region",      "sv_invariant_mass_region_1to1point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_1point5to2_GeV_region",      "sv_invariant_mass_region_1point5to2_GeV_NOSYS"),
    ("sv_invariant_mass_region_2to2point5_GeV_region",      "sv_invariant_mass_region_2to2point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_2point5to3_GeV_region",      "sv_invariant_mass_region_2point5to3_GeV_NOSYS"),
    ("sv_invariant_mass_region_3to3point5_GeV_region",      "sv_invariant_mass_region_3to3point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_3point5to4_GeV_region",      "sv_invariant_mass_region_3point5to4_GeV_NOSYS"),
    ("sv_invariant_mass_region_4to4point5_GeV_region",      "sv_invariant_mass_region_4to4point5_GeV_NOSYS"),
    ("sv_invariant_mass_region_4point5to5_GeV_region",      "sv_invariant_mass_region_4point5to5_GeV_NOSYS"),
    ("sv_invariant_mass_region_5to5point5_GeV_region",      "sv_invariant_mass_region_5to5point5_GeV_NOSYS"),
]

# -------------------------------------------------
# Helpers
# -------------------------------------------------
def flavour_label(f):
    f = abs(int(f))
    if f == 5:
        return "b"
    elif f == 4 or f == 0:
        return "nonb"
    else:
        return None

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
    "jet_pt_new_NOSYS",

    # GN2 tag/probe selections
    "jet_select_GN2v01_FixedCutBEff_65_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_77_NOSYS",

    # weights
    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
] + [branch for _, branch in sv_regions]

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
jet_pt = arr["jet_pt_new_NOSYS"][base_mask]

wp65 = arr["jet_select_GN2v01_FixedCutBEff_65_NOSYS"][base_mask]
wp77 = arr["jet_select_GN2v01_FixedCutBEff_77_NOSYS"][base_mask]

w_event = ak.to_numpy((
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)[base_mask])

sv_region_flags = {
    region_name: ak.to_numpy(arr[branch_name][base_mask])
    for region_name, branch_name in sv_regions
}

# -------------------------------------------------
# Require at least 2 jets in required jet-level branches
# -------------------------------------------------
valid = (
    (ak.num(truth, axis=1) >= 2)
    & (ak.num(jet_pt, axis=1) >= 2)
    & (ak.num(wp65, axis=1) >= 2)
    & (ak.num(wp77, axis=1) >= 2)
)

truth = truth[valid]
jet_pt = jet_pt[valid]
wp65 = wp65[valid]
wp77 = wp77[valid]
w_event = w_event[valid]

for region_name in sv_region_flags:
    sv_region_flags[region_name] = sv_region_flags[region_name][valid]

# -------------------------------------------------
# Define leading/subleading by pT
# probe = leading
# tag   = subleading
# -------------------------------------------------
idx = np.arange(len(truth))

pt0 = ak.to_numpy(jet_pt[:, 0])
pt1 = ak.to_numpy(jet_pt[:, 1])

lead_idx = np.where(pt0 >= pt1, 0, 1)
sublead_idx = np.where(pt0 >= pt1, 1, 0)

probe_idx = lead_idx
tag_idx = sublead_idx

probe_truth = ak.to_numpy(truth[idx, probe_idx])
probe_wp65 = ak.to_numpy(ak.values_astype(wp65[idx, probe_idx], np.int32))
tag_wp77 = ak.to_numpy(ak.values_astype(wp77[idx, tag_idx], np.int32))

probe_flavour = np.array([flavour_label(f) for f in probe_truth], dtype=object)

# -------------------------------------------------
# Calibration-region selection:
# tag passes WP77 and probe passes WP65
# -------------------------------------------------
tag_pass = tag_wp77 == 1
probe_pass = probe_wp65 == 1

fit_mask = (
    tag_pass
    & probe_pass
    & np.isin(probe_flavour, ["b", "nonb"])
)

probe_flavour = probe_flavour[fit_mask]
w_event = w_event[fit_mask]

for region_name in sv_region_flags:
    sv_region_flags[region_name] = sv_region_flags[region_name][fit_mask]

# -------------------------------------------------
# Extract N_b_MC and N_nonb_MC per SV region
# -------------------------------------------------
mc_sv_template_list = []

for region_name, _ in sv_regions:
    region_mask = sv_region_flags[region_name] == 1

    if not np.any(region_mask):
        mc_sv_template_list.append((region_name, 0.0, 0.0, 0.0, 0.0))
        continue

    w_bin = w_event[region_mask]
    flav_bin = probe_flavour[region_mask]

    w_b = w_bin[flav_bin == "b"]
    w_nonb = w_bin[flav_bin == "nonb"]

    N_b, err_b = weighted_yield_and_error(w_b)
    N_nonb, err_nonb = weighted_yield_and_error(w_nonb)

    mc_sv_template_list.append((region_name, N_b, err_b, N_nonb, err_nonb))

# -------------------------------------------------
# Print results
# -------------------------------------------------
print("\n" + "=" * 120)
print("MC SV MASS TEMPLATES FOR PROBE JET")
print("Selection: tag WP77 && probe WP65")
print("Probe = leading jet, tag = subleading jet")
print("Format: region_name, N_b_MC, err_b_MC, N_nonb_MC, err_nonb_MC")
print("=" * 120)

for row in mc_sv_template_list:
    print(f"{row[0]}, {row[1]:.6f}, {row[2]:.6f}, {row[3]:.6f}, {row[4]:.6f}")

print("=" * 120)

# -------------------------------------------------
# Python list for copy/paste into fit
# -------------------------------------------------
print("\nmc_sv_template_list = [")
for row in mc_sv_template_list:
    print(
        f'    ("{row[0]}", {row[1]:.10f}, {row[2]:.10f}, '
        f'{row[3]:.10f}, {row[4]:.10f}),'
    )
print("]")