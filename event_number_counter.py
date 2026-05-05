import uproot

file = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = file["reco"]

n_events = tree.num_entries
print(n_events)

import uproot
import awkward as ak
import numpy as np

fname = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree = "reco"

branches = ["selection_cuts_NOSYS"]

# load data
arr = uproot.open(fname)[tree].arrays(branches)

# apply selection
mask = arr["selection_cuts_NOSYS"] == 1

# count events
n_selected = ak.sum(mask)

print(int(n_selected))