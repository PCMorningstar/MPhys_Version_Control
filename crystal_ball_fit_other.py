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
    "new_truth_mlmbb_NOSYS",
    "new_truth_pTdiff_NOSYS",
    "new_truth_sum_deltaR_NOSYS",
    "new_truth_mllbb_NOSYS",
    "new_truth_mT_ttbar_NOSYS"
]

# --------------------------------------------------
# Crystal Ball function (left-tail version)
# A = amplitude
# mu = mean
# sigma = width
# alpha = tail transition
# n = tail exponent
# --------------------------------------------------
def crystal_ball(x, A, mu, sigma, alpha, n):

    sigma = np.abs(sigma)
    alpha = np.abs(alpha)
    n = np.abs(n)

    t = (x - mu) / sigma

    # Gaussian region
    gaussian = np.exp(-0.5 * t**2)

    # Power-law tail
    A_tail = (n / alpha)**n * np.exp(-0.5 * alpha**2)
    B_tail = n / alpha - alpha

    tail = A_tail * (B_tail - t)**(-n)

    return A * np.where(t > -alpha, gaussian, tail)

# --------------------------------------------------
# Fit Crystal Ball
# --------------------------------------------------
def fit_crystal_ball(data, bins=100):

    data = np.asarray(data)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]  # keep if observable must be positive

    if len(data) < 10:
        return 0.0, 0.0

    counts, edges = np.histogram(data, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    A0     = float(counts.max())
    mu0    = float(np.mean(data))
    sigma0 = max(float(np.std(data)), 1e-2)

    alpha0 = 1.5
    n0     = 5.0

    try:
        popt, _ = curve_fit(
            crystal_ball,
            centers,
            counts,
            p0=[A0, mu0, sigma0, alpha0, n0],
            bounds=(
                [0, -np.inf, 1e-3, 0.1, 0.1],
                [np.inf, np.inf, np.inf, 10.0, 100.0]
            ),
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
# Fit per region (Crystal Ball)
# --------------------------------------------------
results = {}
for region_name, region_mask_func in REGIONS.items():
    mask = region_mask_func(arrays)
    results[region_name] = {}

    for obs in OBSERVABLES:
        data = arrays[obs][mask]

        # Fit Crystal Ball: returns (mu, sigma)
        mu, sigma = fit_crystal_ball(data)
        results[region_name][obs] = (mu, sigma)

# --------------------------------------------------
# Print results
# --------------------------------------------------
for region, obs_dict in results.items():
    print(f"Region: {region}")
    for obs, (mu, sigma) in obs_dict.items():
        print(f"  {obs:<20} : mu = {mu:.2f}, sigma = {sigma:.2f}")
    print()