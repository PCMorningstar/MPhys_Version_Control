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
filename = "purity_v_njets.py"  # Make sure this is a plain data file, not .py

bb_data  = read_efficiency_block(filename, "chi2: jets, purity, error - bb")
bc_data  = read_efficiency_block(filename, "chi2: jets, purity, error - b notb")
nbb_data = read_efficiency_block(filename, "chi2: jets, purity, error - notb b")
nbnb_data = read_efficiency_block(filename, "chi2: jets, purity, error - notb notb")

# Extract columns
nj_bb,  eff_bb,  err_bb  = bb_data.T
nj_bc,  eff_bc,  err_bc  = bc_data.T
nj_nbb, eff_nbb, err_nbb = nbb_data.T
nj_nbnb, eff_nbnb, err_nbnb = nbnb_data.T
# -------------------------------------------------
# Plot
# -------------------------------------------------
plt.figure(figsize=(8, 6))

def plot_with_errorbars(nj, eff, err, color, label):
    plt.step(nj, eff, where="mid", color=color, linewidth=1.5, label=label)
    plt.errorbar(nj, eff, yerr=err, fmt="none", ecolor="black",
                 elinewidth=1, capsize=4, capthick=1)

plot_with_errorbars(nj_bb, eff_bb, err_bb, "blue", r"{b, b}")
plot_with_errorbars(nj_bc, eff_bc, err_bc, "purple", r"{b, Not b}")
plot_with_errorbars(nj_nbb, eff_nbb, err_nbb, "red", r"{Not b, b}")
plot_with_errorbars(nj_nbnb, eff_nbnb, err_nbnb, "orange", r"{Not b, Not b}")

# -------------------------------------------------
# Styling
# -------------------------------------------------
plt.xlabel("Jet Multiplicity")
plt.ylabel("Probability Density")
plt.title(r"General Flavour Composition: Probability Density vs. Jet Multiplicity - Chi2")

plt.grid(True, linestyle=":", linewidth=0.7)
plt.legend(loc="upper left", frameon=False)
plt.xlim(0.0, 10.5)
plt.ylim(0.0, 1.0)

plt.tight_layout()
plt.savefig("general_flavour_composition_purity_v_njets_chi2.png", dpi=300)
plt.show()
plt.close()