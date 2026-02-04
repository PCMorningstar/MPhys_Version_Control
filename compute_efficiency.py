import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ROOT
# --------------------------------------------------
# Numerator and denominator
# --------------------------------------------------
                 # total valid events
N_fail = 8.04468
N_pass = 2.74022
N_tot  = N_fail + N_pass

err_fail = 0.6599
err_success = 0.3904

eff = N_pass / N_tot

# --------------------------------------------------
# Gaussian error propagation for efficiency
# --------------------------------------------------
# Efficiency = N_pass / N_tot
#
# σ_eff^2 = (∂ε/∂N_pass)^2 σ_pass^2 + (∂ε/∂N_tot)^2 σ_tot^2
#
# σ_pass = sqrt(N_pass), σ_tot = sqrt(N_tot)

ratio = 1/pow((N_pass + N_fail), 2)
first_part = pow((N_fail*err_success),2)
last_part = pow((N_pass*err_fail), 2)
together = first_part + last_part

sqrt = np.sqrt(together)    
sigma_eff = ratio*sqrt


# --------------------------------------------------
# Output
# --------------------------------------------------
print("Chi-squared Pairing Outcome (ΔR truth):")
print(f"  Efficiency = {eff:.4f}, {sigma_eff:.4f}")


labels = ["Correct", "Incorrect"]
counts = [N_pass, N_fail]
counts_err = [err_success, err_fail]
colour = ["tab:blue", "tab:red"]

plt.figure(figsize=(6, 4))
plt.bar(
    labels,
    counts,
    color=["tab:blue", "tab:red"],
    yerr = counts_err
)

plt.ylabel("Events")
plt.xlabel("Reconstruction Verdict")
plt.title(
    fr"χ² Reconstruction (Against Truth ΔR)"
)

plt.tight_layout()
plt.savefig("chi2_success_failure.png", dpi=300)
plt.show()
plt.close()

