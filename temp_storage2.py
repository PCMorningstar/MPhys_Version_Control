import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Function to read a block from file
# -------------------------------------------------
def read_efficiency_block(filename, header_keyword):
    data = []
    reading = False

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line:
                reading = False
                continue

            if line.startswith("#"):
                reading = header_keyword.lower() in line.lower()
                continue

            if reading:
                parts = [float(x.strip()) for x in line.split(",")]
                data.append(parts)

    return np.array(data)   # shape (n_rows, 3)


# -------------------------------------------------
# Load data
# -------------------------------------------------
filename = "misms_chi2_hyper_v_njets.py"   # <-- make sure this is a DATA file, not .py

chi2_data  = read_efficiency_block(filename, "chi2")
misms_data = read_efficiency_block(filename, "MISMS")
hyper_data = read_efficiency_block(filename, "HyPER")

# Extract columns
nj_chi2,  eff_chi2,  err_chi2  = chi2_data.T
nj_misms, eff_misms, err_misms = misms_data.T
nj_hyper, eff_hyper, err_hyper = hyper_data.T


# -------------------------------------------------
# Plot
# -------------------------------------------------
plt.figure(figsize=(8, 6))

# χ2
plt.step(
    nj_chi2, eff_chi2,
    where="mid",
    color="blue",
    linewidth=1.5,
    label=r"$\chi^2$"
)

plt.errorbar(
    nj_chi2, eff_chi2,
    yerr=err_chi2,
    fmt="none",
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

# MISMS
plt.step(
    nj_misms, eff_misms,
    where="mid",
    color="red",
    linewidth=1.5,
    label=r"MISMS"
)

plt.errorbar(
    nj_misms, eff_misms,
    yerr=err_misms,
    fmt="none",
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

# HyPER
plt.step(
    nj_hyper, eff_hyper,
    where="mid",
    color="yellow",
    linewidth=1.5,
    label=r"HyPER"
)

plt.errorbar(
    nj_hyper, eff_hyper,
    yerr=err_hyper,
    fmt="none",
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

# -------------------------------------------------
# Styling
# -------------------------------------------------
plt.xlabel("Jet Multiplicity")
plt.ylabel("Full Reconstruction Efficiency")
plt.title("Full Event Reconstruction Efficiency: Technique Comparison")

plt.grid(True, linestyle=":", linewidth=0.7)
plt.legend(frameon=False)

plt.xlim(min(nj_chi2) - 0.5, max(nj_chi2) + 0.5)
plt.ylim(0.0, 1.1)

plt.tight_layout()
plt.savefig("full_recon_eff_v_njets.png", dpi=300)
plt.show()
plt.close()