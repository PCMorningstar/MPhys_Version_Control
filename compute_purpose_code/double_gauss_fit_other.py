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
# Double-Gaussian model (shared mean)
# --------------------------------------------------
def double_gaussian(x, A1, A2, mu, sigma1, sigma2):
    sigma1 = np.abs(sigma1)
    sigma2 = np.abs(sigma2)
    return (
        A1 * np.exp(-((x - mu) ** 2) / (2.0 * sigma1**2)) +
        A2 * np.exp(-((x - mu) ** 2) / (2.0 * sigma2**2))
    )

# --------------------------------------------------
# Fit Double Gaussian (replaces fit_gaussian)
# Returns (mu, sigma_eff) where sigma_eff is mixture RMS around mu
# --------------------------------------------------
def fit_double_gaussian(data, bins=100):
    data = np.asarray(data)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]   # keep your original behaviour

    if len(data) < 10:
        return 0.0, 0.0

    counts, edges = np.histogram(data, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # fit only non-empty bins
    m = counts > 0
    x = centers[m]
    y = counts[m]
    if len(x) < 5:
        mu0, sigma0 = float(np.mean(data)), max(float(np.std(data)), 1e-2)
        return mu0, sigma0

    mu0 = float(np.mean(data))
    sigma0 = max(float(np.std(data)), 1e-2)
    A0 = float(y.max())

    # initial guesses: narrow+wide
    p0 = [
        0.7 * A0,          # A1
        0.3 * A0,          # A2
        mu0,               # mu
        0.6 * sigma0,      # sigma1
        1.8 * sigma0,      # sigma2
    ]

    bounds = (
        [0.0, 0.0, -np.inf, 1e-3, 1e-3],
        [np.inf, np.inf, np.inf, np.inf, np.inf],
    )

    try:
        popt, _ = curve_fit(
            double_gaussian, x, y,
            p0=p0,
            bounds=bounds,
            maxfev=100000
        )
        A1, A2, mu, s1, s2 = popt
        A1, A2 = float(A1), float(A2)
        mu = float(mu)
        s1, s2 = float(abs(s1)), float(abs(s2))

        # enforce ordering (optional)
        if s2 < s1:
            s1, s2 = s2, s1
            A1, A2 = A2, A1

        # effective sigma = mixture RMS
        denom = A1 + A2
        if denom <= 0:
            return mu0, sigma0
        w1 = A1 / denom
        w2 = A2 / denom
        sigma_eff = np.sqrt(w1 * s1**2 + w2 * s2**2)

        return mu, float(sigma_eff)

    except Exception:
        return mu0, sigma0

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
        mu, sigma = fit_double_gaussian(data)
        results[region_name][obs] = (mu, sigma)

# --------------------------------------------------
# Print results
# --------------------------------------------------
for region, obs_dict in results.items():
    print(f"Region: {region}")
    for obs, (mu, sigma) in obs_dict.items():
        print(f"  {obs:<20} : mu = {mu:.2f}, sigma = {sigma:.2f}")
    print()