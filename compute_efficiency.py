import uproot
import numpy as np

# --------------------------------------------------
# Open the ROOT ntuple produced by FastFrames
# --------------------------------------------------
ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
f = uproot.open(ntuple_path)
tree = f["reco"]

# --------------------------------------------------
# Load needed base arrays
# --------------------------------------------------
# Event weights
weight_mc  = tree["weight_mc_NOSYS"].array(library="np")
weight_pu  = tree["weight_pileup_NOSYS"].array(library="np")
weight_lep = tree["weight_leptonSF_tight_NOSYS"].array(library="np")
weight_jvt = tree["weight_jvt_effSF_NOSYS"].array(library="np")
weights_all = weight_mc * weight_pu * weight_lep * weight_jvt

# Selection variables
selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
jet_size       = tree["jet_size_NOSYS"].array(library="np")

# Chi2 indexed branch (present in your output ntuple)
chi_branch = "new_chi_indexed_jets_NOSYS"

# --------------------------------------------------
# Weighted efficiency function
# --------------------------------------------------
def efficiency_fail_success(values, weights):
    """
    Calculates weighted efficiency:
      - values >= 0 are considered valid
      - values == 1 => success
      - values == 0 => fail
    """
    values  = np.asarray(values)
    weights = np.asarray(weights)
    valid_mask = (values == 0) | (values == 1)
    values  = values[valid_mask]
    weights = weights[valid_mask]

    pass_mask = values == 1
    fail_mask = values == 0

    N_pass = np.sum(weights[pass_mask])
    N_fail = np.sum(weights[fail_mask])
    N_tot  = N_pass + N_fail
    if N_tot == 0:
        return 0.0, 0.0, 0, 0, 0

    err_pass = np.sqrt(np.sum(weights[pass_mask]**2))
    err_fail = np.sqrt(np.sum(weights[fail_mask]**2))
    err_tot  = np.sqrt(err_pass**2 + err_fail**2)

    eff = N_pass / N_tot
    dE_dNpass = 1.0 / N_tot
    dE_dNtot  = -N_pass / (N_tot**2)
    sigma_eff = np.sqrt((dE_dNpass*err_pass)**2 + (dE_dNtot*err_tot)**2)

    return eff, sigma_eff, N_pass, N_fail, N_tot

# --------------------------------------------------
# Define regions exactly as in your YAML
# --------------------------------------------------
regions = {
    "sec_region"     : (selection_cuts == 1),
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
# Loop and compute efficiencies
# --------------------------------------------------
eff_results = {}

for region_name, mask in regions.items():
    # Skip if mask is empty
    if not np.any(mask):
        eff_results[region_name] = {
            "eff": 0.0, "sigma": 0.0, "N_pass": 0, "N_fail": 0, "N_tot": 0
        }
        continue

    # Extract chi-indexed values for valid events in this region
    chi_vals = tree[chi_branch].array(library="np")[mask]
    region_weights = weights_all[mask]

    eff, sigma, Np, Nf, Nt = efficiency_fail_success(chi_vals, region_weights)
    eff_results[region_name] = {
        "eff": eff, "sigma": sigma,
        "N_pass": Np, "N_fail": Nf, "N_tot": Nt
    }

# --------------------------------------------------
# Print summary
# --------------------------------------------------
print(f"{'Region':<15} | {'efficiency':>10} ± {'sigma':>6} | {'N_pass':>8} {'N_fail':>8} {'N_tot':>8}")
print("-"*60)
for region_name, res in eff_results.items():
    print(f"{region_name:<15} | {res['eff']:10.4f} ± {res['sigma']:6.4f} | "
          f"{int(res['N_pass']):8d} {int(res['N_fail']):8d} {int(res['N_tot']):8d}")