import uproot
import numpy as np
from scipy.optimize import curve_fit

# -------------------------------
# Gaussian function
# -------------------------------
def gauss(x, mu, sigma, A):
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)

# -------------------------------
# Open ROOT file
# -------------------------------
f = uproot.open("/eos/user/p/pzeman/FFTutorial2/output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = f["reco"]

# -------------------------------
# Collect arrays
# -------------------------------
truth_m_lpb = [tree[f"truth_m_lpb{i}_NOSYS"].array() for i in range(1,6)]
truth_m_lmbb = [tree[f"truth_m_lmbb{i}_NOSYS"].array() for i in range(1,6)]

# -------------------------------
# Helper to flatten jagged arrays
# -------------------------------
def safe_flatten(arr):
    arr = np.asarray(arr)
    if arr.ndim == 1:
        return arr
    else:
        return np.concatenate(arr)

# -------------------------------
# Function to fit and get mean/std
# -------------------------------
def fit_gaussian(arr, branch_name):
    data = safe_flatten(arr)
    data = data[data > 0]  # remove invalid entries (0 = invalid)
    
    if len(data) == 0:
        print(f"{branch_name}: No valid entries")
        return

    hist, bins = np.histogram(data, bins=50)
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    
    # Initial guesses: mean, std, amplitude
    p0 = [np.mean(data), np.std(data), np.max(hist)]
    
    try:
        popt, _ = curve_fit(gauss, bin_centers, hist, p0=p0)
        mu, sigma, _ = popt
        print(f"{branch_name}: mean = {mu:.3f}, std = {sigma:.3f}")
        
    except RuntimeError:
        print(f"{branch_name}: Gaussian fit failed")

# -------------------------------
# Loop over branches
# -------------------------------
for i, arr in enumerate(truth_m_lpb, start=1):
    fit_gaussian(arr, f"truth_m_lpb{i}")

for i, arr in enumerate(truth_m_lmbb, start=1):
    fit_gaussian(arr, f"truth_m_lmbb{i}")
