import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Read efficiency block from file
# -------------------------------------------------
def read_efficiency_block(filename, header):
    data = []
    reading = False

    header_lower = header.lower()

    with open(filename) as f:
        for line in f:
            line = line.strip()

            # Empty line stops reading
            if not line:
                reading = False
                continue

            # Header line: check if it contains header keyword
            if line.startswith("#"):
                reading = header_lower in line.lower()
                continue

            # Collect numeric lines if reading
            if reading:
                parts = line.split(",")
                if len(parts) != 3:
                    continue  # skip malformed lines
                try:
                    dR, eff, err = map(float, map(str.strip, parts))
                    data.append((dR, eff, err))
                except ValueError:
                    continue

    if not data:
        raise ValueError(f"No data found for header '{header}' in {filename}")

    data = np.array(data)
    return data[:, 0], data[:, 1], data[:, 2]

# -------------------------------------------------
# Load all efficiency blocks
# -------------------------------------------------
filename = "new_efficiency_v_jets_chi2.py"  # Make sure it's a .data file, not .py

nj_lpb, eff_lpb, err_lpb = read_efficiency_block(
    filename, "Number of Jets, efficiency, error - Default Terms (mlb_plus & mlb_minus)"
)

nj_pTdiff, eff_pTdiff, err_pTdiff = read_efficiency_block(
    filename, "Number of Terms, efficiency, error - Transverse Momentum Difference added (pTdiff)"
)

nj_sum_deltaR, eff_sum_deltaR, err_sum_deltaR = read_efficiency_block(
    filename, "Number of Terms, efficiency, error - Sum of Delta R added(sum_deltaR)"
)

nj_mllbb, eff_mllbb, err_mllbb = read_efficiency_block(
    filename, "Number of Terms, efficiency, error - Visible Mass of the ttbar system added (mllbb) (mllbb)"
)

nj_mT_ttbar, eff_mT_ttbar, err_mT_ttbar = read_efficiency_block(
    filename, "Number of Terms, efficiency, error - Transverse Mass of the ttbar system added (mT & MET) (mT_ttbar)"
)

# -------------------------------------------------
# Plot
# -------------------------------------------------
plt.figure(figsize=(7,5))

# Helper to plot step + errorbar
def plot_efficiency(nj, eff, err, color, label):
    plt.step(nj, eff, where="mid", color=color, linewidth=1.5, label=label)
    plt.errorbar(nj, eff, yerr=err, fmt="none", ecolor="black",
                 elinewidth=1, capsize=4, capthick=1)

plot_efficiency(nj_lpb, eff_lpb, err_lpb, "blue", "Default Terms")
plot_efficiency(nj_pTdiff, eff_pTdiff, err_pTdiff, "red", "pTdiff Added")
plot_efficiency(nj_sum_deltaR, eff_sum_deltaR, err_sum_deltaR, "brown", r"$\Delta R_{sum}$ Added")
plot_efficiency(nj_mllbb, eff_mllbb, err_mllbb, "purple", "mllbb Added")
plot_efficiency(nj_mT_ttbar, eff_mT_ttbar, err_mT_ttbar, "orange", "mT_ttbar Added")

plt.xlabel("Jet Multiplicity")
plt.ylabel("Full Reconstruction Efficiency")
plt.title("Full Reconstruction Efficiency vs. Jet Multiplicity [Cumulative Terms]")
plt.grid(True, which="both", linestyle=":", linewidth=0.7)
plt.legend(frameon=False)
plt.xlim(left=0)
plt.ylim(bottom=0.0, top=1.0)
plt.tight_layout()
plt.savefig("new_efficiency_v_jets_chi2.png", dpi=300)
plt.show()
plt.close()