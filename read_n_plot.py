import numpy as np
import matplotlib.pyplot as plt

def read_efficiency_block(filename, header):
    data = []
    reading = False

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line:
                reading = False
                continue

            if line.startswith("#"):
                reading = header in line
                continue

            if reading:
                dR, eff, err = map(float, line.split(","))
                data.append((dR, eff, err))

    data = np.array(data)
    return data[:,0], data[:,1], data[:,2]



nj_lpb, eff_lpb, err_lpb = read_efficiency_block(
    "new_efficiency_v_jets_chi2.py", "Number of Jets, NEW (ALL TERMS) Reconstruction efficiency, Error - chi2"
)

nj_lmbb, eff_lmbb, err_lmbb = read_efficiency_block(
    "new_efficiency_v_jets_chi2.py", "Number of Jets, NEW Reconstruction efficiency - all terms - (Constant), Error - chi2"
)

plt.figure(figsize=(7,5))

# ─────────────────────────────────────
# l+ b  (step + error bars)
# ─────────────────────────────────────
plt.step(
    nj_lpb, eff_lpb,
    where="mid",
    color="blue",
    linewidth=1.5,
    label=r"Non-Constant $\mu$ & $\sigma$"
)

plt.errorbar(
    nj_lpb, eff_lpb, yerr=err_lpb,
    fmt="none",              # <-- no markers, no line
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

# ─────────────────────────────────────
# l− b̄
# ─────────────────────────────────────
plt.step(
    nj_lmbb, eff_lmbb,
    where="mid",
    color="red",
    linewidth=1.5,
    label=r"Constant $\mu$ & $\sigma$"
)

plt.errorbar(
    nj_lmbb, eff_lmbb, yerr=err_lmbb,
    fmt="none",
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

plt.xlabel("Jet Multiplicity")
plt.ylabel("Full Reconstruction efficiency")
plt.title(r"$\chi^2$ Efficiency vs. Jet Multiplicity [Non-Constant vs. Constant $\mu$ & $\sigma$]")

plt.grid(True, which="both", linestyle=":", linewidth=0.7)
plt.legend(frameon=False)

plt.xlim(left=0)
plt.ylim(bottom=0.0, top=1.0)

plt.tight_layout()
plt.savefig("ncVc_newe_new_chi2_v_njets.png", dpi=300)
plt.show()
plt.close()
