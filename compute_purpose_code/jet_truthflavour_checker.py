import uproot
import numpy as np

ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree = uproot.open(ntuple_path)["reco"]

# branches listed in your YAML
selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
jet_size       = tree["jet_size_NOSYS"].array(library="np")
ordered_truth_flav = tree["ordered_jet_truth_flavour_extended_NOSYS"].array(library="np")

# weights listed in your YAML
weights_all = (
    tree["weight_mc_NOSYS"].array(library="np") *
    tree["weight_pileup_NOSYS"].array(library="np") *
    tree["weight_leptonSF_tight_NOSYS"].array(library="np") *
    tree["weight_jvt_effSF_NOSYS"].array(library="np")
)

# regions as in your ntuples:regions list
regions = {
    "sec_region": (selection_cuts == 1),  # matches your comment: use same as selection
    **{f"{nj}jets_region": (selection_cuts == 1) & (jet_size == nj) for nj in range(2, 11)}
}

def print_flavour_codes_for_region(region_name, mask):
    if not np.any(mask):
        print(f"{region_name}: no events")
        return

    # flatten only selected events in this region
    flavs = np.concatenate(ordered_truth_flav[mask])
    uniq, cnt = np.unique(flavs, return_counts=True)

    print(f"{region_name}: events={int(np.sum(mask))}, jets={int(len(flavs))}")
    print("  flavour codes and counts:")
    for u, c in zip(uniq, cnt):
        print(f"    {int(u)}: {int(c)}")
    print("  unique codes:", uniq)
    print()

# print for all ntuple regions
for rname, rmask in regions.items():
    print_flavour_codes_for_region(rname, rmask)