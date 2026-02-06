import uproot
import numpy as np

# -------------------------------
# Open ROOT file and get tree
# -------------------------------
f = uproot.open("/eos/user/p/pzeman/FFTutorial2/output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = f["reco"]

# -------------------------------
# Load arrays
# -------------------------------
chi_lpb1 = tree["chi_lpb_one_NOSYS"].array()
chi_lpb2 = tree["chi_lpb_two_NOSYS"].array()
chi_lpb3 = tree["chi_lpb_three_NOSYS"].array()
chi_lpb4 = tree["chi_lpb_four_NOSYS"].array()
chi_lpb5 = tree["chi_lpb_five_NOSYS"].array()

chi_lmbb1 = tree["chi_lmbb_one_NOSYS"].array()
chi_lmbb2 = tree["chi_lmbb_two_NOSYS"].array()
chi_lmbb3 = tree["chi_lmbb_three_NOSYS"].array()
chi_lmbb4 = tree["chi_lmbb_four_NOSYS"].array()
chi_lmbb5 = tree["chi_lmbb_five_NOSYS"].array()

# Load all event weights
weight_mc = tree["weight_mc_NOSYS"].array()
weight_pu = tree["weight_pileup_NOSYS"].array()
weight_lep = tree["weight_leptonSF_tight_NOSYS"].array()
weight_jvt = tree["weight_jvt_effSF_NOSYS"].array()

# Compute total event weight
weights = weight_mc * weight_pu * weight_lep * weight_jvt

# Optional: apply selection mask if you have one
if "selection_cuts_NOSYS" in tree.keys():
    selection_mask = tree["selection_cuts_NOSYS"].array() == 1
    chi_lpb = chi_lpb[selection_mask]
    chi_lmbb = chi_lmbb[selection_mask]
    weights = weights[selection_mask]

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

    # Only valid events (1 or 2)
    valid_mask = (values > 0)
    values = values[valid_mask]
    weights = weights[valid_mask]

    pass_mask = (values == 2)
    fail_mask = (values == 1)

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
eff_lpb1, sigma_lpb1, Np_lpb1, Nf_lpb1, Nt_lpb1 = efficiency_fail_success(chi_lpb1, weights)
eff_lpb2, sigma_lpb2, Np_lpb2, Nf_lpb2, Nt_lpb2 = efficiency_fail_success(chi_lpb2, weights)
eff_lpb3, sigma_lpb3, Np_lpb3, Nf_lpb3, Nt_lpb3 = efficiency_fail_success(chi_lpb3, weights)
eff_lpb4, sigma_lpb4, Np_lpb4, Nf_lpb4, Nt_lpb4 = efficiency_fail_success(chi_lpb4, weights)
eff_lpb5, sigma_lpb5, Np_lpb5, Nf_lpb5, Nt_lpb5 = efficiency_fail_success(chi_lpb5, weights)

eff_lmbb1, sigma_lmbb1, Np_lmbb1, Nf_lmbb1, Nt_lmbb1 = efficiency_fail_success(chi_lmbb1, weights)
eff_lmbb2, sigma_lmbb2, Np_lmbb2, Nf_lmbb2, Nt_lmbb2 = efficiency_fail_success(chi_lmbb2, weights)
eff_lmbb3, sigma_lmbb3, Np_lmbb3, Nf_lmbb3, Nt_lmbb3 = efficiency_fail_success(chi_lmbb3, weights)
eff_lmbb4, sigma_lmbb4, Np_lmbb4, Nf_lmbb4, Nt_lmbb4 = efficiency_fail_success(chi_lmbb4, weights)
eff_lmbb5, sigma_lmbb5, Np_lmbb5, Nf_lmbb5, Nt_lmbb5 = efficiency_fail_success(chi_lmbb5, weights)


# -------------------------------
# Print results
# -------------------------------
print("Efficiency for chi_lpb_NOSYS:")
print(f"{eff_lpb1:.4f}, {sigma_lpb1:.4f}")
print(f"{eff_lpb2:.4f}, {sigma_lpb2:.4f}")
print(f"{eff_lpb3:.4f}, {sigma_lpb3:.4f}")
print(f"{eff_lpb4:.4f}, {sigma_lpb4:.4f}")
print(f"{eff_lpb5:.4f}, {sigma_lpb5:.4f}\n")

print("Efficiency for chi_lmbb_NOSYS:")
print(f"{eff_lmbb1:.4f}, {sigma_lmbb1:.4f}")
print(f"{eff_lmbb2:.4f}, {sigma_lmbb2:.4f}")
print(f"{eff_lmbb3:.4f}, {sigma_lmbb3:.4f}")
print(f"{eff_lmbb4:.4f}, {sigma_lmbb4:.4f}")
print(f"{eff_lmbb5:.4f}, {sigma_lmbb5:.4f}")