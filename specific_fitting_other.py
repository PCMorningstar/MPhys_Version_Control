import uproot
import numpy as np
import matplotlib.pyplot as plt
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

OBSERVABLES_Gaussian = [
    "new_truth_mlpb_NOSYS",
    "new_truth_mlmbb_NOSYS",
]

OBSERVABLES_CrystalBall = [
    "new_truth_sum_deltaR_NOSYS",
    "new_truth_mllbb_NOSYS",
    "new_truth_mT_ttbar_NOSYS"
]

OBSERVABLES_DoubleGaussian = [
    "new_truth_pTdiff_NOSYS"
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
    data = data[data > 0.0]  # keep if observable must be positive

    if len(data) < 10:
        return 0.0, 0.0, None, None, None

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
        A, mu, sigma = popt
    except Exception:
        A, mu, sigma = A0, mu0, sigma0

    return mu, sigma, counts, edges, A

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
    out = np.empty_like(t, dtype=float)

    # Gaussian region
    m_gaus = t > -alpha
    out[m_gaus] = np.exp(-0.5 * t[m_gaus]**2)

    # Tail region
    A_tail = (n / alpha)**n * np.exp(-0.5 * alpha**2)
    B_tail = n / alpha - alpha

    m_tail = ~m_gaus
    u = B_tail - t[m_tail]           # should be > 0
    u = np.maximum(u, 1e-12)         # guard domain
    out[m_tail] = A_tail * u**(-n)

    return A * out

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

arrays_Gaussian = {branch: tree[branch].array(library="np")
          for branch in ["jet_size_NOSYS", "selection_cuts_NOSYS"] + OBSERVABLES_Gaussian}
arrays_CrystalBall = {branch: tree[branch].array(library="np")
          for branch in ["jet_size_NOSYS", "selection_cuts_NOSYS"] + OBSERVABLES_CrystalBall}
arrays_DoubleGaussian = {branch: tree[branch].array(library="np")
          for branch in ["jet_size_NOSYS", "selection_cuts_NOSYS"] + OBSERVABLES_DoubleGaussian}

# --------------------------------------------------
# Fit per region - NO PLOTTING
# --------------------------------------------------
for region_name, region_mask_func in REGIONS.items():
    mask = region_mask_func(arrays_Gaussian)

    for obs_Gaussian in OBSERVABLES_Gaussian:
        data_Gaussian = arrays_Gaussian[obs_Gaussian][mask]
        mu_Gaussian, sigma_Gaussian, counts, edges, A = fit_gaussian(data_Gaussian)
        print(f"Region: {region_name}, Observable: {obs_Gaussian} - Gaussian")
        print(f"  mu = {mu_Gaussian:.2f}, sigma = {sigma_Gaussian:.2f}")

    mask = region_mask_func(arrays_CrystalBall)
    for obs_CrystalBall in OBSERVABLES_CrystalBall:
        data_CrystalBall = arrays_CrystalBall[obs_CrystalBall][mask]
        mu_CrystalBall, sigma_CrystalBall = fit_crystal_ball(data_CrystalBall)
        print(f"Region: {region_name}, Observable: {obs_CrystalBall} - Crystal Ball")
        print(f"  mu = {mu_CrystalBall:.2f}, sigma = {sigma_CrystalBall:.2f}")

    mask = region_mask_func(arrays_DoubleGaussian)
    for obs_DoubleGaussian in OBSERVABLES_DoubleGaussian:
        data_DoubleGaussian = arrays_DoubleGaussian[obs_DoubleGaussian][mask]
        mu_DoubleGaussian, sigma_DoubleGaussian = fit_double_gaussian(data_DoubleGaussian)
        print(f"Region: {region_name}, Observable: {obs_DoubleGaussian} - Double Gaussian")
        print(f"  mu = {mu_DoubleGaussian:.2f}, sigma = {sigma_DoubleGaussian:.2f}")