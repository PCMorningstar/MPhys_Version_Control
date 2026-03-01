import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --------------------------------------------------
# Config: regions and observables
# --------------------------------------------------
REGIONS = {
    "2jets_region": lambda arrays: (arrays["selection_cuts_NOSYS"] == 1) & (arrays["jet_size_NOSYS"] == 2),
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
# Gaussian model
# --------------------------------------------------
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

# --------------------------------------------------
# Fit Gaussian robustly
# --------------------------------------------------
def fit_gaussian(data, bins=100, allow_negative=False):
    data = np.asarray(data)
    data = data[np.isfinite(data)]

    # Keep the old behaviour for most observables,
    # but allow negatives for pTdiff (or anything you flag).
    if not allow_negative:
        data = data[data > 0.0]

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
# Load ROOT file and arrays
# --------------------------------------------------
file = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")
tree = file["reco"]

arrays = {branch: tree[branch].array(library="np")
          for branch in ["jet_size_NOSYS", "selection_cuts_NOSYS"] + OBSERVABLES}

# --------------------------------------------------
# Fit per region and plot for 2jets
# --------------------------------------------------
for region_name, region_mask_func in REGIONS.items():
    mask = region_mask_func(arrays)

    for obs in OBSERVABLES:
        data = arrays[obs][mask]

        allow_negative = (obs == "new_truth_pTdiff_NOSYS")
        mu, sigma, counts, edges, A = fit_gaussian(data, allow_negative=allow_negative)

        # Print results
        print(f"Region: {region_name}, Observable: {obs}")
        print(f"  mu = {mu:.2f}, sigma = {sigma:.2f}")

        # Plot histogram + Gaussian fit
        centers = 0.5 * (edges[:-1] + edges[1:])
        plt.figure(figsize=(7, 5))
        plt.bar(centers, counts, width=(edges[1] - edges[0]), alpha=0.6, label="Data")

        x_fit = np.linspace(edges[0], edges[-1], 500)
        plt.plot(
            x_fit,
            gaussian(x_fit, A, mu, sigma),
            "r-",
            label=f"Gaussian fit\nmu={mu:.2f}, sigma={sigma:.2f}"
        )

        # Only force negative visibility for pTdiff
        if allow_negative:
            xmin, xmax = edges[0], edges[-1]
            if xmin >= 0:
                # If histogram range is still non-negative (rare), force symmetric around 0
                xmax = max(abs(xmax), abs(xmin))
                xmin = -xmax
            plt.xlim(xmin, xmax)

        plt.xlabel(obs)
        plt.ylabel("Counts")
        plt.title(f"{region_name} - {obs}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"fit_{region_name}_{obs}.png")
        plt.close()