import uproot
import numpy as np

# --------------------------------------------------
# Open ROOT ntuple
# --------------------------------------------------
ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
f = uproot.open(ntuple_path)
tree = f["reco"]

# --------------------------------------------------
# Selection variables
# --------------------------------------------------
selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
jet_size       = tree["jet_size_NOSYS"].array(library="np")

# Weighted event (FastFrames already multiplies all factors)
weights_all = tree["weight_mc_NOSYS"].array(library="np") * \
              tree["weight_pileup_NOSYS"].array(library="np") * \
              tree["weight_leptonSF_tight_NOSYS"].array(library="np") * \
              tree["weight_jvt_effSF_NOSYS"].array(library="np")

# Indexed branches
chi_branch  = "new_chi_indexed_jets_NOSYS"
mdrs_branch = "new_mdrs_indexed_jets_NOSYS"
misms_branch= "new_misms_indexed_jets_NOSYS"

# --------------------------------------------------
# Weighted efficiency function
# --------------------------------------------------
def weighted_efficiency(values, weights):
    """
    Compute weighted efficiency:
    eff = sum(weights of success events) / sum(weights of valid events)
    Only considers values 0 (fail) or 1 (success)
    """
    values  = np.asarray(values)
    weights = np.asarray(weights)
    valid_mask = (values == 0) | (values == 1)
    if not np.any(valid_mask):
        return 0.0, 0.0, 0.0, 0.0, 0.0

    values  = values[valid_mask]
    weights = weights[valid_mask]

    pass_mask = values == 1
    fail_mask = values == 0

    w_pass = np.sum(weights[pass_mask])
    w_fail = np.sum(weights[fail_mask])
    w_tot  = w_pass + w_fail

    eff = w_pass / w_tot if w_tot > 0 else 0.0

    # Weighted binomial uncertainty
    sigma_eff = np.sqrt( np.sum((weights[pass_mask])**2) ) / w_tot if w_tot > 0 else 0.0

    return eff, sigma_eff, w_pass, w_fail, w_tot

# --------------------------------------------------
# Regions as in YAML
# --------------------------------------------------
regions = {
    "2jets_region"   : (selection_cuts == 1) & (jet_size == 2),
    "3jets_region"   : (selection_cuts == 1) & (jet_size == 3),
    "4jets_region"   : (selection_cuts == 1) & (jet_size == 4),
    "5jets_region"   : (selection_cuts == 1) & (jet_size == 5),
    "6jets_region"   : (selection_cuts == 1) & (jet_size == 6),
    "7jets_region"   : (selection_cuts == 1) & (jet_size == 7),
    "8jets_region"   : (selection_cuts == 1) & (jet_size == 8),
    "9jets_region"   : (selection_cuts == 1) & (jet_size == 9),
    "10jets_region"  : (selection_cuts == 1) & (jet_size == 10)
}

# --------------------------------------------------
# Compute weighted efficiencies
# --------------------------------------------------
eff_results = {}
for region_name, mask in regions.items():
    if not np.any(mask):
        eff_results[region_name] = {
            "chi2": (0,0), "mdrs": (0,0), "misms": (0,0),
            "w_pass": 0, "w_fail": 0, "w_tot": 0
        }
        continue

    region_weights = weights_all[mask]
    chi_vals  = tree[chi_branch].array(library="np")[mask]
    mdrs_vals = tree[mdrs_branch].array(library="np")[mask]
    misms_vals= tree[misms_branch].array(library="np")[mask]

    eff_chi, sigma_chi, wpass_chi, wfail_chi, wtot_chi = weighted_efficiency(chi_vals, region_weights)
    eff_mdrs, sigma_mdrs, wpass_mdrs, wfail_mdrs, wtot_mdrs = weighted_efficiency(mdrs_vals, region_weights)
    eff_misms, sigma_misms, wpass_misms, wfail_misms, wtot_misms = weighted_efficiency(misms_vals, region_weights)

    eff_results[region_name] = {
        "chi2": (eff_chi, sigma_chi),
        "mdrs": (eff_mdrs, sigma_mdrs),
        "misms": (eff_misms, sigma_misms),
        "w_pass": wpass_chi, "w_fail": wfail_chi, "w_tot": wtot_chi
    }

# --------------------------------------------------
# Print weighted summary
# --------------------------------------------------
print(f"{'Region':<15} | {'chi2':>10} ± {'sigma':>6} | {'mdrs':>10} ± {'sigma':>6} | {'misms':>10} ± {'sigma':>6} | {'w_pass':>10} {'w_fail':>10} {'w_tot':>10}")
print("-"*110)
for region_name, res in eff_results.items():
    print(f"{region_name:<15} | "
          f"{res['chi2'][0]:10.4f} ± {res['chi2'][1]:6.4f} | ")
print("-"*110)