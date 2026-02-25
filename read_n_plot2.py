import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Function to read a CSV-style data block from a file
# -------------------------------------------------
def read_efficiency_block(filename, header_keyword):
    """
    Reads a block of efficiency data from a text file.
    Expected format:
    # <header_keyword>
    nj, eff, err
    """
    data = []
    reading = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                reading = False
                continue

            # Check for header
            if line.startswith("#"):
                reading = header_keyword.lower() in line.lower()
                continue

            if reading:
                # Convert line to floats
                parts = [float(x.strip()) for x in line.split(",")]
                if len(parts) == 3:
                    data.append(parts)

    if not data:
        raise ValueError(f"No data found for header '{header_keyword}' in {filename}")

    return np.array(data)  # shape (n_rows, 3)


# -------------------------------------------------
# Load data
# -------------------------------------------------
filename = "misms_chi2_hyper_v_njets.py"  # Make sure this is a plain data file, not .py

chi2_data  = read_efficiency_block(filename, "Number of Jets, NEW Reconstruction efficiency - selected terms , Error - chi2")
mdrs_data  = read_efficiency_block(filename, "Number of Jets, NEW Reconstruction efficiency, Error - MDRS")
misms_data = read_efficiency_block(filename, "Number of Jets, NEW Reconstruction efficiency, Error - MISMS")
hyper_data = read_efficiency_block(filename, "Number of Jets, NEW Reconstruction efficiency, Error - HyPER")

# Extract columns
nj_chi2,  eff_chi2,  err_chi2  = chi2_data.T
nj_mdrs,  eff_mdrs,  err_mdrs  = mdrs_data.T
nj_misms, eff_misms, err_misms = misms_data.T
nj_hyper, eff_hyper, err_hyper = hyper_data.T


# -------------------------------------------------
# Plot
# -------------------------------------------------
plt.figure(figsize=(8, 6))

def plot_with_errorbars(nj, eff, err, color, label):
    plt.step(nj, eff, where="mid", color=color, linewidth=1.5, label=label)
    plt.errorbar(nj, eff, yerr=err, fmt="none", ecolor="black",
                 elinewidth=1, capsize=4, capthick=1)

plot_with_errorbars(nj_chi2, eff_chi2, err_chi2, "blue", r"$\chi^2$")
plot_with_errorbars(nj_mdrs, eff_mdrs, err_mdrs, "purple", r"MDRS")
plot_with_errorbars(nj_misms, eff_misms, err_misms, "red", r"MISMS")
plot_with_errorbars(nj_hyper, eff_hyper, err_hyper, "orange", r"HyPER")

# -------------------------------------------------
# Styling
# -------------------------------------------------
plt.xlabel("Jet Multiplicity")
plt.ylabel("Full Reconstruction Efficiency")
plt.title("Full Event Reconstruction Efficiency: Technique Comparison")

plt.grid(True, linestyle=":", linewidth=0.7)
plt.legend(frameon=False)

plt.xlim(0.0, max(nj_chi2) + 0.5)
plt.ylim(0.0, 1.0)

plt.tight_layout()
plt.savefig("full_recon_eff_v_njets.png", dpi=300)
plt.show()
plt.close()