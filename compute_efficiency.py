import uproot
import numpy as np

# -------------------------------
# Open ROOT file and get tree
# -------------------------------
f = uproot.open("/eos/user/p/pzeman/FFTutorial2/output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = f["reco"]

# -------------------------------
# Load arrays - Unbiased Truth
# -------------------------------
chi_lpb = tree["new_chi_indexed_NOSYS"].array()
#chi_lmbb = tree["new_misms_indexed_NOSYS"].array()
# -------------------------------
# Load arrays - Biased Truth
# -------------------------------
chi_lpb1 = tree["chi_lpb_NOSYS"].array()
#chi_lmbb1 = tree["misms_lpb_NOSYS"].array()

# Load all event weights
weight_mc = tree["weight_mc_NOSYS"].array()
weight_pu = tree["weight_pileup_NOSYS"].array()
weight_lep = tree["weight_leptonSF_tight_NOSYS"].array()
weight_jvt = tree["weight_jvt_effSF_NOSYS"].array()

# Compute total event weight
weights = weight_mc * weight_pu * weight_lep * weight_jvt

# -------------------------------
# Efficiency function
# -------------------------------
def efficiency_fail_success(values, weights):
    """
    Calculates weighted efficiency with:
    0 = invalid (ignored)
    1 = fail
    2 = success
    """
    values = np.asarray(values)
    weights = np.asarray(weights)

    # Only valid events (0 or 1)
    valid_mask = (values >= 0)
    values = values[valid_mask]
    weights = weights[valid_mask]

    pass_mask = (values == 1)
    fail_mask = (values == 0)

    N_pass = np.sum(weights[pass_mask])
    N_fail = np.sum(weights[fail_mask])
    N_tot  = N_pass + N_fail

    # Weighted errors
    err_pass = np.sqrt(np.sum(weights[pass_mask]**2))
    err_fail = np.sqrt(np.sum(weights[fail_mask]**2))
    err_tot  = np.sqrt(err_pass**2 + err_fail**2)

    # Efficiency
    eff = N_pass / N_tot

    # Gaussian error propagation
    dE_dNpass = 1.0 / N_tot
    dE_dNtot  = -N_pass / (N_tot**2)
    sigma_eff = np.sqrt((dE_dNpass * err_pass)**2 + (dE_dNtot * err_tot)**2)

    return eff, sigma_eff, N_pass, N_fail, N_tot

# -------------------------------
# Compute efficiencies
# -------------------------------
eff_lpb, sigma_lpb, Np_lpb, Nf_lpb, Nt_lpb = efficiency_fail_success(chi_lpb, weights)
#eff_lmbb, sigma_lmbb, Np_lmbb, Nf_lmbb, Nt_lmbb = efficiency_fail_success(chi_lmbb, weights)
eff_lpb1, sigma_lpb1, Np_lpb1, Nf_lpb1, Nt_lpb1 = efficiency_fail_success(chi_lpb1, weights)
#eff_lmbb1, sigma_lmbb1, Np_lmbb1, Nf_lmbb1, Nt_lmbb1 = efficiency_fail_success(chi_lmbb1, weights)

# -------------------------------
# Print results
# -------------------------------
#print("Efficiency for lpb_NOSYS:")
print("Efficiency for chi2 (Old Truth):")
print(f"{eff_lpb1:.4f}, {sigma_lpb1:.4f}")
print("Efficiency for chi2 (Better Truth):")
print(f"{eff_lpb:.4f}, {sigma_lpb:.4f}\n")

#print("Efficiency for lmbb_NOSYS:")
#print("Efficiency for misms (Old Truth):")
#print(f"{eff_lmbb1:.4f}, {sigma_lmbb1:.4f}")
#print("Efficiency for misms (Better Truth):")
#print(f"{eff_lmbb:.4f}, {sigma_lmbb:.4f}")
