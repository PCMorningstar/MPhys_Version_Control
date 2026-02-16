import uproot
import numpy as np

# --- REQUIRED ON lxplus ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# --------------------------------------------------
# Load ROOT file
# --------------------------------------------------
file = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = file["reco"]

m_eb = tree["truth_pTdiff_NOSYS"].array(library="np")
m_eb1 = tree["truth_sum_deltaR_NOSYS"].array(library="np")
m_eb2 = tree["truth_mllbb_NOSYS"].array(library="np")
m_eb3 = tree["truth_mT_ttbar_NOSYS"].array(library="np")

# --------------------------------------------------
# Clean arrays
# --------------------------------------------------
def clean(arr):
    arr = arr[np.isfinite(arr)]
    return arr[arr > 0.0]

m_ebc = clean(m_eb)
m_ebc1 = clean(m_eb1)
m_ebc2 = clean(m_eb2)
m_ebc3 = clean(m_eb3)

# --------------------------------------------------
# Gaussian model
# --------------------------------------------------
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

# --------------------------------------------------
# Gaussian fit
# --------------------------------------------------
def fit_gaussian(data, bins=100, xmin=None, xmax=None):
    data = np.array(data)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]

    # Use quantiles to define initial guess if xmin/xmax not given
    q_low, q_high = np.percentile(data, [5, 95])
    if xmin is None:
        xmin = 0
    if xmax is None:
        xmax = max(data.max(), q_high + 0.5*(q_high-q_low))

    counts, edges = np.histogram(data, bins=bins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Robust initial guesses
    A0 = counts.max()
    mu0 = centers[np.argmax(counts)]
    sigma0 = 0.5*(q_high - q_low)  # 1σ approx from 16-84 percentile

    try:
        popt, pcov = curve_fit(
            gaussian,
            centers,
            counts,
            p0=[A0, mu0, sigma0],
            bounds=([0, 0, 1e-3], [np.inf, np.inf, np.inf]),
            maxfev=50000
        )
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        # fallback
        popt = [A0, mu0, sigma0]
        perr = [0, 0, 0]

    return popt, perr

# --------------------------------------------------
# Fits
# --------------------------------------------------
popt_ej1, perr_ej1 = fit_gaussian(m_ebc)
A_ej1, mu_ej1, sigma_ej1 = popt_ej1
popt_ej2, perr_ej2 = fit_gaussian(m_ebc1)
A_ej2, mu_ej2, sigma_ej2 = popt_ej2
popt_ej3, perr_ej3 = fit_gaussian(m_ebc2)
A_ej3, mu_ej3, sigma_ej3 = popt_ej3
popt_ej4, perr_ej4 = fit_gaussian(m_ebc3)
A_ej4, mu_ej4, sigma_ej4 = popt_ej4

# --------------------------------------------------
# Print results
# --------------------------------------------------
print("truth_pTdiff_NOSYS:")
print(f"mu, sigma = {mu_ej1:.2f}, {sigma_ej1:.2f}")
print("")
print("truth_sum_deltaR_NOSYS:")
print(f"mu, sigma = {mu_ej2:.2f}, {sigma_ej2:.2f}")
print("")
print("truth_mllbb_NOSYS:")
print(f"mu, sigma = {mu_ej3:.2f}, {sigma_ej3:.2f}")
print("")
print("truth_mT_ttbar_NOSYS:")
print(f"mu, sigma = {mu_ej4:.2f}, {sigma_ej4:.2f}")