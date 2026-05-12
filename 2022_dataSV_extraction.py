import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/data_0_2022_data.root"
tree = "reco"

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

branches = [
    "selection_cuts_NOSYS",
    "jet_size_NOSYS",
    "jet_pt_new_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_65_pt_ordered_NOSYS",
    "jet_select_GN2v01_FixedCutBEff_77_pt_ordered_NOSYS",
] + [branch for _, branch in sv_regions]

with uproot.open(fname) as f:
    arr = f[tree].arrays(branches, library="ak")

base_mask = (
    (arr["selection_cuts_NOSYS"] == 1)
    & (arr["jet_size_NOSYS"] >= 2)
)

jet_pt = arr["jet_pt_new_NOSYS"][base_mask]
wp65 = arr["jet_select_GN2v01_FixedCutBEff_65_pt_ordered_NOSYS"][base_mask]
wp77 = arr["jet_select_GN2v01_FixedCutBEff_77_pt_ordered_NOSYS"][base_mask]

sv_region_flags = {
    region_name: ak.to_numpy(arr[branch_name][base_mask])
    for region_name, branch_name in sv_regions
}

valid = (
    (ak.num(jet_pt, axis=1) >= 2)
    & (ak.num(wp65, axis=1) >= 2)
    & (ak.num(wp77, axis=1) >= 2)
)

valid_np = ak.to_numpy(valid)

jet_pt = jet_pt[valid]
wp65 = wp65[valid]
wp77 = wp77[valid]

for region_name in sv_region_flags:
    sv_region_flags[region_name] = sv_region_flags[region_name][valid_np]

idx = np.arange(len(jet_pt))

pt0 = ak.to_numpy(jet_pt[:, 0])
pt1 = ak.to_numpy(jet_pt[:, 1])

lead_idx = np.where(pt0 >= pt1, 0, 1)
sublead_idx = np.where(pt0 >= pt1, 1, 0)

probe_idx = lead_idx
tag_idx = sublead_idx

probe_wp65 = ak.to_numpy(ak.values_astype(wp65[idx, probe_idx], np.int32))
tag_wp77 = ak.to_numpy(ak.values_astype(wp77[idx, tag_idx], np.int32))

fit_mask = (
    (tag_wp77 == 1)
    & (probe_wp65 == 1)
)

for region_name in sv_region_flags:
    sv_region_flags[region_name] = sv_region_flags[region_name][fit_mask]

print("SV region, N_data, err_data")

for region_name, _ in sv_regions:
    region_mask = sv_region_flags[region_name] == 1

    N_data = int(np.sum(region_mask))
    err_data = np.sqrt(N_data)

    print(f"{region_name}, {N_data}, {err_data:.6f}")