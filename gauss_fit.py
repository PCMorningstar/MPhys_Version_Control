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

m_eb = tree["truth_m_lpb_NOSYS"].array(library="np")
m_mbb = tree["truth_m_lmbb_NOSYS"].array(library="np")


# --------------------------------------------------
# Clean arrays
# --------------------------------------------------
def clean(arr):
    arr = arr[np.isfinite(arr)]
    return arr[arr > 0.0]

m_ebc = clean(m_eb)
m_mbbc = clean(m_mbb)

# --------------------------------------------------
# Gaussian model
# --------------------------------------------------
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

# --------------------------------------------------
# Gaussian fit
# --------------------------------------------------
def fit_gaussian(data, bins=100, xmin=0, xmax=400):
    counts, edges = np.histogram(data, bins=bins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    A0 = counts.max()
    mu0 = centers[np.argmax(counts)]
    sigma0 = np.std(data)

    popt, pcov = curve_fit(
        gaussian,
        centers,
        counts,
        p0=[A0, mu0, sigma0],
        bounds=([0, 0, 1e-3], [np.inf, xmax, np.inf]),
        maxfev=10000
    )

    perr = np.sqrt(np.diag(pcov))
    return popt, perr

# --------------------------------------------------
# Fits
# --------------------------------------------------
popt_ej1, perr_ej1 = fit_gaussian(m_ebc)
A_ej1, mu_ej1, sigma_ej1 = popt_ej1

popt_mj1, perr_mj1 = fit_gaussian(m_mbbc)
A_mj1, mu_mj1, sigma_mj1 = popt_mj1

# --------------------------------------------------
# Print results
# --------------------------------------------------
print("m_lpb_truth:")
print(f"  mu, sigma    = {mu_ej1:.2f}, {sigma_ej1:.2f} GeV")

print("m_lmbbar_truth:")
print(f"  mu, sigma    = {mu_mj1:.2f}, {sigma_mj1:.2f} GeV")
