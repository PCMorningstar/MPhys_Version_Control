import uproot
import numpy as np

# --------------------------------------------------
# Open ROOT ntuple
# --------------------------------------------------
ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
with uproot.open(ntuple_path) as f:
    tree = f["reco"]

    # --------------------------------------------------
    # Core branches
    # --------------------------------------------------
    selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
    jet_size = tree["jet_size_NOSYS"].array(library="np")

    # Weighted event
    weights_all = (
        tree["weight_mc_NOSYS"].array(library="np")
        * tree["weight_pileup_NOSYS"].array(library="np")
        * tree["weight_leptonSF_tight_NOSYS"].array(library="np")
        * tree["weight_jvt_effSF_NOSYS"].array(library="np")
    )

    # Indexed branches
    chi_vals_all = tree["new_chi_indexed_jets_NOSYS"].array(library="np")
    mdrs_vals_all = tree["new_mdrs_indexed_jets_NOSYS"].array(library="np")
    misms_vals_all = tree["new_misms_indexed_jets_NOSYS"].array(library="np")


# --------------------------------------------------
# Weighted efficiency function
# --------------------------------------------------
def weighted_efficiency(values, weights):
    """
    Weighted efficiency using only values equal to 0 or 1.
    Returns: eff, sigma_eff, w_pass, w_fail, w_tot
    """
    values = np.asarray(values)
    weights = np.asarray(weights)

    valid_mask = (values == 0) | (values == 1)
    if not np.any(valid_mask):
        return 0.0, 0.0, 0.0, 0.0, 0.0

    values = values[valid_mask]
    weights = weights[valid_mask]

    pass_mask = values == 1
    fail_mask = values == 0

    w_pass = np.sum(weights[pass_mask])
    w_fail = np.sum(weights[fail_mask])
    w_tot = w_pass + w_fail

    if w_tot <= 0:
        return 0.0, 0.0, w_pass, w_fail, w_tot

    eff = w_pass / w_tot

    sumw2_pass = np.sum(weights[pass_mask] ** 2)
    sumw2_fail = np.sum(weights[fail_mask] ** 2)

    var_eff = (
        ((1.0 - eff) ** 2) * sumw2_pass
        + (eff ** 2) * sumw2_fail
    ) / (w_tot ** 2)

    sigma_eff = np.sqrt(max(var_eff, 0.0))
    return eff, sigma_eff, w_pass, w_fail, w_tot


# --------------------------------------------------
# Regions as a function of jet multiplicity
# Matches the YAML:
# selection_cuts_NOSYS == 1 && jet_size_NOSYS == N
# --------------------------------------------------
regions = {
    "2jets_region": (selection_cuts == 1) & (jet_size == 2),
    "3jets_region": (selection_cuts == 1) & (jet_size == 3),
    "4jets_region": (selection_cuts == 1) & (jet_size == 4),
    "5jets_region": (selection_cuts == 1) & (jet_size == 5),
    "6jets_region": (selection_cuts == 1) & (jet_size == 6),
    "7jets_region": (selection_cuts == 1) & (jet_size == 7),
    "8jets_region": (selection_cuts == 1) & (jet_size == 8),
    "9jets_region": (selection_cuts == 1) & (jet_size == 9),
    "10jets_region": (selection_cuts == 1) & (jet_size == 10),
    "geq2jets_region": (selection_cuts == 1) & (jet_size >= 2),
}


# --------------------------------------------------
# Compute weighted efficiencies
# --------------------------------------------------
eff_results = {}

for region_name, mask in regions.items():
    if not np.any(mask):
        eff_results[region_name] = {
            "chi2": (0.0, 0.0),
            "mdrs": (0.0, 0.0),
            "misms": (0.0, 0.0),
            "w_pass": 0.0,
            "w_fail": 0.0,
            "w_tot": 0.0,
        }
        continue

    region_weights = weights_all[mask]
    chi_vals = chi_vals_all[mask]
    mdrs_vals = mdrs_vals_all[mask]
    misms_vals = misms_vals_all[mask]

    eff_chi, sigma_chi, wpass_chi, wfail_chi, wtot_chi = weighted_efficiency(chi_vals, region_weights)
    eff_mdrs, sigma_mdrs, _, _, _ = weighted_efficiency(mdrs_vals, region_weights)
    eff_misms, sigma_misms, _, _, _ = weighted_efficiency(misms_vals, region_weights)

    eff_results[region_name] = {
        "chi2": (eff_chi, sigma_chi),
        "mdrs": (eff_mdrs, sigma_mdrs),
        "misms": (eff_misms, sigma_misms),
        "w_pass": wpass_chi,
        "w_fail": wfail_chi,
        "w_tot": wtot_chi,
    }


# --------------------------------------------------
# Print weighted summary
# --------------------------------------------------
print(
    f"{'Region':<15} | "
    f"{'chi2':>10} ± {'sigma':>8} | "
    f"{'mdrs':>10} ± {'sigma':>8} | "
    f"{'misms':>10} ± {'sigma':>8}"
)
print("-" * 80)

for region_name, res in eff_results.items():
    print(
        f"{region_name:<15} | "
        f"{res['chi2'][0]:10.4f} ± {res['chi2'][1]:8.4f} | "
        f"{res['mdrs'][0]:10.4f} ± {res['mdrs'][1]:8.4f} | "
        f"{res['misms'][0]:10.4f} ± {res['misms'][1]:8.4f}"
    )

print("-" * 80)