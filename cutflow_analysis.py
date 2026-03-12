import uproot
import awkward as ak
import numpy as np

filename = "output_ntuples/ttll_601230_mc23a_fullsim.root"
treename = "reco"

branches = [
    "selection_cuts_NOSYS",

    "electron_selections_paper_NOSYS",
    "muon_selections_paper_NOSYS",
    "jet_selections_paper_NOSYS",
    "dilepton_selections_paper_NOSYS",

    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
]

with uproot.open(filename) as f:
    tree = f[treename]
    arr = tree.arrays(branches, library="ak")

# ------------------------------------------------
# Event weight
# ------------------------------------------------
w = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
)

# ------------------------------------------------
# Masks
# ------------------------------------------------
raw_mask = (arr["selection_cuts_NOSYS"] == 1) | (arr["selection_cuts_NOSYS"] == 0)

electron_mask = arr["electron_selections_paper_NOSYS"] == 1
muon_mask     = arr["muon_selections_paper_NOSYS"] == 1
jet_mask      = arr["jet_selections_paper_NOSYS"] == 1
dilepton_mask = arr["dilepton_selections_paper_NOSYS"] == 1

# Your order
mask_raw      = raw_mask
mask_electron = mask_raw & electron_mask
mask_muon     = mask_electron & muon_mask
mask_jet      = mask_muon & jet_mask
mask_dilepton = mask_jet & dilepton_mask

labels = ["Raw", "Electron", "Muon", "Jet", "Dilepton"]
masks  = [mask_raw, mask_electron, mask_muon, mask_jet, mask_dilepton]

# ------------------------------------------------
# Yield + error
# ------------------------------------------------
print("\nCutflow yields (weighted):\n")

for label, mask in zip(labels, masks):

    yield_val = float(ak.sum(w[mask]))
    err_val   = float(np.sqrt(ak.sum(w[mask] ** 2)))

    print(f"{label:10s} : {yield_val:14.6f} ± {err_val:10.6f}")


