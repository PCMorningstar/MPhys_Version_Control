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
# Crystal Ball function (left-tail version)
# A = amplitude
# mu = mean
# sigma = width
# alpha = tail transition
# n = tail exponent
# --------------------------------------------------
def crystal_ball(x, A, mu, sigma, alpha, n):
    x = np.asarray(x, dtype=float)

    sigma = max(abs(float(sigma)), 1e-12)
    alpha = max(abs(float(alpha)), 1e-12)
    n = max(abs(float(n)), 1e-12)

    t = (x - mu) / sigma

    # Precompute constants
    A_tail = (n / alpha) ** n * np.exp(-0.5 * alpha**2)
    B_tail = n / alpha - alpha

    # Allocate output
    y = np.empty_like(t, dtype=float)

    # Gaussian core: t > -alpha
    core_mask = t > -alpha
    y[core_mask] = np.exp(-0.5 * t[core_mask] ** 2)

    # Left tail: t <= -alpha
    tail_mask = ~core_mask
    if np.any(tail_mask):
        den = B_tail - t[tail_mask]

        # Numerical protection against zero/negative denominator
        den = np.maximum(den, 1e-300)

        y[tail_mask] = A_tail * den ** (-n)

    return A * y


# --------------------------------------------------
# Fit Crystal Ball
# --------------------------------------------------
def fit_crystal_ball(data, bins=100):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]  # keep if observable must be positive

    if data.size < 10:
        return 0.0, 0.0

    counts, edges = np.histogram(data, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Keep only non-empty bins for the fit
    nonzero = counts > 0
    counts_fit = counts[nonzero]
    centers_fit = centers[nonzero]

    if counts_fit.size < 5:
        mu0 = float(np.mean(data))
        sigma0 = max(float(np.std(data)), 1e-2)
        return mu0, sigma0

    A0 = float(counts_fit.max())
    mu0 = float(np.mean(data))
    sigma0 = max(float(np.std(data)), 1e-2)
    alpha0 = 1.5
    n0 = 5.0

    try:
        popt, _ = curve_fit(
            crystal_ball,
            centers_fit,
            counts_fit,
            p0=[A0, mu0, sigma0, alpha0, n0],
            bounds=(
                [0.0, -np.inf, 1e-3, 0.1, 0.5],
                [np.inf, np.inf, np.inf, 10.0, 100.0]
            ),
            maxfev=50000
        )

        mu = float(popt[1])
        sigma = abs(float(popt[2]))

    except Exception:
        mu = mu0
        sigma = sigma0

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