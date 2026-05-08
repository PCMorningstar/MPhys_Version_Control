import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------
# Open ROOT file
# ------------------------------------------------
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

# ------------------------------------------------
# Branches
# ------------------------------------------------
branches = [
    "weight_mc_NOSYS",
    "weight_pileup_NOSYS",
    "weight_leptonSF_tight_NOSYS",
    "weight_jvt_effSF_NOSYS",
]

# ------------------------------------------------
# Read arrays
# ------------------------------------------------
with uproot.open(filename) as f:
    tree = f[treename]
    arr = tree.arrays(branches, library="ak")

# ------------------------------------------------
# Total event weight
# ------------------------------------------------
total_weights = (
    arr["weight_mc_NOSYS"]
    * arr["weight_pileup_NOSYS"]
    * arr["weight_leptonSF_tight_NOSYS"]
    * arr["weight_jvt_effSF_NOSYS"]
    * norm
)

# ------------------------------------------------
# Convert to numpy
# ------------------------------------------------
weights_np = ak.to_numpy(total_weights)

# ------------------------------------------------
# Print first 10 weights
# ------------------------------------------------
print("\nFirst 10 total event weights:\n")

for i, w in enumerate(weights_np[:10]):
    print(f"{i:2d} : {w:.12e}")

# ------------------------------------------------
# Plot distribution
# ------------------------------------------------
plt.figure(figsize=(8,6))

plt.hist(
    weights_np,
    bins=100,
    histtype="step",
)

plt.xlabel("Total Event Weight")
plt.ylabel("Number of Events")
plt.title("Total Event-Weight Distribution")

plt.yscale("log")
plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig("total_event_weight_distribution_debug.png")
plt.close()

print("\nSaved plot:")
print("total_event_weight_distribution_debug.png")