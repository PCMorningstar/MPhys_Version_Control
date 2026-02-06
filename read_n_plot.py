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
    "efficiency_v_dR.py", "lpb"
)

nj_lmbb, eff_lmbb, err_lmbb = read_efficiency_block(
    "efficiency_v_dR.py", "lmbb"
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
    label=r"$b$-lepton pairing"
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
    label=r"$\bar{b}$-lepton pairing"
)

plt.errorbar(
    nj_lmbb, eff_lmbb, yerr=err_lmbb,
    fmt="none",
    ecolor="black",
    elinewidth=1,
    capsize=4,
    capthick=1
)

plt.xlabel("Truth-matched ΔR separation")
plt.ylabel("Reconstruction efficiency")
plt.title(r"$\chi^2$ Efficiency vs. Truth ΔR")

plt.grid(True, which="both", linestyle=":", linewidth=0.7)
plt.legend(frameon=False)

plt.xlim(left=0)
plt.ylim(bottom=0.8, top=0.86)

plt.tight_layout()
plt.savefig("efficiency_v_dR.png", dpi=300)
plt.show()
plt.close()
