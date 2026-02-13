import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Function to read a block from file
# ---------------------------
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
                reading = header_keyword in line
                continue
            if reading:
                # Split on comma, strip spaces, convert to float
                parts = [float(x.strip()) for x in line.split(",")]
                data.append(parts)
    return np.array(data)  # shape: (n_rows, 3)

# ---------------------------
# Load all three techniques
# ---------------------------
filename = "efficiency_v_bbjets.py"  # change to your actual file path

chi2_data  = read_efficiency_block(filename, "chi2")
misms_data = read_efficiency_block(filename, "MISMS")
hyper_data = read_efficiency_block(filename, "HyPER")

# Extract columns
nj_chi2, eff_chi2, err_chi2 = chi2_data[:,0], chi2_data[:,1], chi2_data[:,2]
nj_misms, eff_misms, err_misms = misms_data[:,0], misms_data[:,1], misms_data[:,2]
nj_hyper, eff_hyper, err_hyper = hyper_data[:,0], hyper_data[:,1], hyper_data[:,2]

# ---------------------------
# Plot grouped bar chart
# ---------------------------
bar_width = 0.25
x = np.arange(len(nj_chi2))  # positions for the bars

plt.figure(figsize=(10,6))

# Plot bars with error bars
plt.bar(x - bar_width, eff_chi2, width=bar_width, yerr=err_chi2, capsize=4,
        label=r"$\chi^2$", color="blue", alpha=0.7)
plt.bar(x, eff_misms, width=bar_width, yerr=err_misms, capsize=4,
        label=r"MISMS", color="red", alpha=0.7)
plt.bar(x + bar_width, eff_hyper, width=bar_width, yerr=err_hyper, capsize=4,
        label=r"HyPER", color="purple", alpha=0.7)

# Labels and title
plt.xticks(x, nj_chi2.astype(int))  # Jet multiplicity as ticks
plt.xlabel("Jet Multiplicity")
plt.ylabel("Reconstruction Efficiency")
plt.title(r"[$\bar{t}$] Event Reconstruction Efficiency: Technique Comparison")
plt.ylim(0,1)
plt.grid(True, which="both", linestyle=":", linewidth=0.7)
plt.legend(frameon=False)

plt.tight_layout()
plt.savefig("lmbb_recon_comparison.png", dpi=300)
plt.show()
plt.close()
