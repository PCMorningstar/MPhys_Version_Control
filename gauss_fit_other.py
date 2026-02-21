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

variables = {
    "Truth-Matched m(lp,b) [GeV]": tree["new_truth_mlpb_NOSYS"].array(library="np"),
    "Truth-Matched m(lm,bbar) [GeV]": tree["new_truth_mlmbb_NOSYS"].array(library="np"),
    "Truth-Matched pTdiff [GeV]": tree["new_truth_pTdiff_NOSYS"].array(library="np"),
    "Truth-Matched sum_deltaR [rad]": tree["new_truth_sum_deltaR_NOSYS"].array(library="np"),
    "Truth-Matched m(lp,lm,b,bbar) [GeV]": tree["new_truth_mllbb_NOSYS"].array(library="np"),
    "Truth-Matched mT(t,tbar) [GeV]": tree["new_truth_mT_ttbar_NOSYS"].array(library="np"),
}

# --------------------------------------------------
# Clean arrays
# --------------------------------------------------
def clean(arr):
    arr = arr[np.isfinite(arr)]
    return arr[arr > 0.0]

# --------------------------------------------------
# Gaussian model
# --------------------------------------------------
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

# --------------------------------------------------
# Gaussian fit function
# --------------------------------------------------
def fit_gaussian(data, bins=100, xmin=None, xmax=None):
    data = np.array(data)
    data = data[np.isfinite(data)]
    data = data[data > 0.0]

    q_low, q_high = np.percentile(data, [5, 95]) # Best against outliers
    if xmin is None:
        xmin = 0
    if xmax is None:
        xmax = max(data.max(), q_high + 0.5*(q_high-q_low))

    counts, edges = np.histogram(data, bins=bins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    A0 = counts.max()
    mu0 = centers[np.argmax(counts)]
    sigma0 = 0.5*(q_high - q_low)

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
        popt = [A0, mu0, sigma0]
        perr = [0, 0, 0]

    return popt, perr, xmin, xmax

# --------------------------------------------------
# Loop over variables and plot
# --------------------------------------------------
for label, array in variables.items():
    data = clean(array)
    popt, perr, xmin, xmax = fit_gaussian(data)
    A, mu, sigma = popt

    print(f"{label}: mu = {mu:.2f}, sigma = {sigma:.2f}")

    x_fit = np.linspace(xmin, xmax, 1000)

    plt.figure(figsize=(7,5))
    plt.hist(
        data, bins=100, range=(xmin, xmax),
        histtype="step", linewidth=1.5, color="blue", label="Distribution"
    )
    plt.plot(
        x_fit, gaussian(x_fit, *popt),
        color="red", linewidth=1.5,
        label=rf"Gaussian Fit ($\mu={mu:.2f}$, $\sigma={sigma:.2f}$)"
    )

    plt.xlabel(f"{label}")
    plt.ylabel("Events")
    plt.title(f"Gaussian Fit of {label} [2 Jets]")
    plt.grid(True, which="both", linestyle=":", linewidth=0.7)
    plt.legend(frameon=False)
    plt.xlim(xmin, xmax)
    plt.tight_layout()
    plt.savefig(f"gauss_fit_{label}.png", dpi=300)
    plt.close()