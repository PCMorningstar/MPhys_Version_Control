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
filename = "chi2_algorithm_analysis.py"  # Make sure this is a plain data file, not .py

epsilon_both_ord_data  = read_efficiency_block(filename, "epsilon_both_ord(n)")
epsilon_1_data  = read_efficiency_block(filename, "epsilon_1(n)")
epsilon_0_data  = read_efficiency_block(filename, "epsilon_0(n)")
epsilon_unord_data  = read_efficiency_block(filename, "epsilon_unord(n)")

# Extract columns
nj_epsilon_both_ord,  eff_epsilon_both_ord,  err_epsilon_both_ord  = epsilon_both_ord_data.T
nj_epsilon_1,  eff_epsilon_1,  err_epsilon_1  = epsilon_1_data.T
nj_epsilon_0, eff_epsilon_0, err_epsilon_0 = epsilon_0_data.T
nj_epsilon_unord, eff_epsilon_unord, err_epsilon_unord = epsilon_unord_data.T
# -------------------------------------------------
# Plot
# -------------------------------------------------
plt.figure(figsize=(8, 6))

def plot_with_errorbars(nj, eff, err, color, label):
    plt.step(nj, eff, where="mid", color=color, linewidth=1.5, label=label)
    plt.errorbar(nj, eff, yerr=err, fmt="none", ecolor="black",
                 elinewidth=1, capsize=4, capthick=1)

plot_with_errorbars(nj_epsilon_both_ord, eff_epsilon_both_ord, err_epsilon_both_ord, "blue", r"Both Correct")
plot_with_errorbars(nj_epsilon_1, eff_epsilon_1, err_epsilon_1, "purple", r"One Correct")
plot_with_errorbars(nj_epsilon_0, eff_epsilon_0, err_epsilon_0, "red", r"Neither Correct")
plot_with_errorbars(nj_epsilon_unord, eff_epsilon_unord, err_epsilon_unord, "orange", r"Swapped Order")

# -------------------------------------------------
# Styling
# -------------------------------------------------
plt.xlabel("Jet Multiplicity")
plt.ylabel("Full Algorithm Efficiency")
plt.title(r"Full Algorithm Efficiency vs. Jet Multiplicity - Chi2")

plt.grid(True, linestyle=":", linewidth=0.7)
plt.legend(loc="upper right", frameon=False)
plt.xlim(0.0, 10.5)
plt.ylim(0.0, 1.0)

plt.tight_layout()
plt.savefig("full_algorithm_efficiency_v_njets_chi2.png", dpi=300)
plt.show()
plt.close()