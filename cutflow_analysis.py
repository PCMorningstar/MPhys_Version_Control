import uproot
import awkward as ak
import numpy as np

filename = "output_ntuples/ttll_601230_mc23a_fullsim.root"
treename = "reco"

# ------------------------------------------------
# MC normalisation
# ------------------------------------------------
luminosity = 29300.0
xsec = 85.482
filter_eff = 1.0
kfactor = 1.138433852
sum_of_weights = 4268786417.0

norm = luminosity * xsec * filter_eff * kfactor / sum_of_weights

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
# Total event weight
# ------------------------------------------------
w = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
    * norm
)

# ------------------------------------------------
# Masks
# ------------------------------------------------
mask_raw = (arr["selection_cuts_NOSYS"] == 1) | (arr["selection_cuts_NOSYS"] == 0)

mask_electron = mask_raw & (arr["electron_selections_paper_NOSYS"] == 1)
mask_muon     = mask_electron & (arr["muon_selections_paper_NOSYS"] == 1)
mask_jet      = mask_muon & (arr["jet_selections_paper_NOSYS"] == 1)
mask_dilepton = mask_jet & (arr["dilepton_selections_paper_NOSYS"] == 1)

labels = ["Raw", "Electron", "Muon", "Jet", "Dilepton"]
masks  = [mask_raw, mask_electron, mask_muon, mask_jet, mask_dilepton]

# ------------------------------------------------
# Yield + error
# ------------------------------------------------
print("\nCutflow yields (luminosity-normalised):\n")
print(f"Normalisation factor = {norm:.12e}\n")

for label, mask in zip(labels, masks):
    yield_val = float(ak.sum(w[mask]))
    err_val = float(np.sqrt(ak.sum(w[mask] ** 2)))

    print(f"{label:10s} : {yield_val:14.6f} ± {err_val:10.6f}")