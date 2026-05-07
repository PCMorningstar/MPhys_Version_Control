import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# Open ROOT file
file = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")

# Access tree
tree = file["reco"]

# Read event weights
weights = tree["weight_mc_NOSYS"].array(library="ak")

# Convert to numpy
weights_np = ak.to_numpy(weights)

# Plot
plt.figure(figsize=(8,6))

plt.hist(
    weights_np,
    bins=100,
    histtype="step",
)

plt.xlabel("Event weight")
plt.ylabel("Events")
plt.title("Event-weight distribution")

plt.yscale("log")   # usually useful
plt.grid(True)

plt.savefig("event_weight_distribution_debug.png")
plt.close()