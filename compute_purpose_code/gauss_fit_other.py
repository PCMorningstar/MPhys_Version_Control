import uproot
import numpy as np
from scipy.optimize import curve_fit

# --------------------------------------------------
# Config: regions and observables
# --------------------------------------------------
REGIONS = {
    "2jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 2),
    "3jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 3),
    "4jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 4),
    "5jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 5),
    "6jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 6),
    "7jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 7),
    "8jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 8),
    "9jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 9),
    "10jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 10),
}

OBSERVABLES = [
    "new_truth_mlpb_NOSYS",
    "new_truth_mlmbb_NOSYS"
]

# --------------------------------------------------
# Gaussian model
# --------------------------------------------------
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

# --------------------------------------------------
# Fit Gaussian robustly
# --------------------------------------------------
def fit_gaussian(data, bins=100):
    data = np.asarray(data)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]
    if len(data) < 10:
        return 0.0, 0.0

    counts, edges = np.histogram(data, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    A0, mu0, sigma0 = float(counts.max()), float(np.mean(data)), float(np.std(data))
    sigma0 = max(sigma0, 1e-2)

    try:
        popt, _ = curve_fit(
            gaussian, centers, counts,
            p0=[A0, mu0, sigma0],
            bounds=([0, -np.inf, 1e-3], [np.inf, np.inf, np.inf]),
            maxfev=50000
        )
        mu, sigma = popt[1], popt[2]
    except Exception:
        mu, sigma = mu0, sigma0

    return mu, sigma

# --------------------------------------------------
# Load ROOT file and arrays
# --------------------------------------------------
file = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = file["reco"]

arrays = {branch: tree[branch].array(library="np") for branch in ["jet_size_NOSYS","selection_cuts_NOSYS"] + OBSERVABLES}

# --------------------------------------------------
# Fit per region
# --------------------------------------------------
results = {}
for region_name, region_mask_func in REGIONS.items():
    mask = region_mask_func(arrays)
    results[region_name] = {}
    for obs in OBSERVABLES:
        data = arrays[obs][mask]
        mu, sigma = fit_gaussian(data)
        results[region_name][obs] = (mu, sigma)

# --------------------------------------------------
# Print results
# --------------------------------------------------
for region, obs_dict in results.items():
    print(f"Region: {region}")
    for obs, (mu, sigma) in obs_dict.items():
        print(f"  {obs:<20} : mu = {mu:.2f}, sigma = {sigma:.2f}")
    print()