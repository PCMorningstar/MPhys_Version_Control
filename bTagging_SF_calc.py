import numpy as np
from scipy.optimize import minimize_scalar

# Columns: SV centre, yield, error
# Real data
data_yield = np.array([
    [-0.25,   32.000000,   5.656854],
    [ 0.25,   29.000000,   5.385165],
    [ 0.75,  243.000000,  15.588457],
    [ 1.25,  936.000000,  30.594117],
    [ 1.75, 2059.000000,  45.376205],
    [ 2.25, 3046.000000,  55.190579],
    [ 2.75, 3379.000000,  58.129167],
    [ 3.25, 2884.000000,  53.702886],
    [ 3.75, 2064.000000,  45.431267],
    [ 4.25, 1137.000000,  33.719431],
    [ 4.75,  556.000000,  23.579652],
    [ 5.25,  237.000000,  15.394804],
], dtype=float)

# MC data
b_yield = np.array([
    [-0.25,   41.007889,   1.560346],
    [ 0.25,   26.688442,   1.252819],
    [ 0.75,  268.798218,   4.010045],
    [ 1.25, 1010.887878,   7.747545],
    [ 1.75, 2300.191162,  11.700373],
    [ 2.25, 3494.957031,  14.429359],
    [ 2.75, 3931.815430,  15.309745],
    [ 3.25, 3409.349609,  14.267761],
    [ 3.75, 2327.449219,  11.786798],
    [ 4.25, 1322.819702,   8.882799],
    [ 4.75,  625.590210,   6.121982],
    [ 5.25,  256.905914,   3.928616],
], dtype=float)

nonb_yield = np.array([
    [-0.25, 0.077817, 0.059005],
    [ 0.25, 0.000000, 0.000000],
    [ 0.75, 0.072826, 0.072826],
    [ 1.25, 0.570535, 0.181255],
    [ 1.75, 0.574948, 0.181241],
    [ 2.25, 1.089280, 0.253113],
    [ 2.75, 0.712352, 0.207429],
    [ 3.25, 1.403130, 0.293136],
    [ 3.75, 0.450786, 0.172465],
    [ 4.25, 0.099179, 0.070988],
    [ 4.75, 0.000000, 0.000000],
    [ 5.25, 0.053373, 0.053373],
], dtype=float)

N_data = data_yield[:, 1]
s_data = data_yield[:, 2]

N_b = b_yield[:, 1]
s_b = b_yield[:, 2]

N_nonb = nonb_yield[:, 1]
s_nonb = nonb_yield[:, 2]


def chi2(sf):
    model = sf * N_b + N_nonb
    sigma2 = s_data**2 + (sf * s_b)**2 + s_nonb**2

    return np.sum((N_data - model)**2 / sigma2)


result = minimize_scalar(
    chi2,
    bounds=(0.0, 2.0),
    method="bounded",
)

sf_best = result.x
chi2_min = result.fun
ndof = len(N_data) - 1

print(f"SF = {sf_best:.6f}")
print(f"chi2_min = {chi2_min:.6f}")
print(f"ndof = {ndof}")
print(f"chi2/ndof = {chi2_min / ndof:.6f}")

# =================================================
# Error analysis - statistical since SF not obtained from direct algebraic analysis
# We want a high confidence 68.3 "standard"
# https://scikit-hep.org/iminuit/reference.html?utm_source=chatgpt.com
# =================================================

from scipy.optimize import brentq

target = chi2_min + 1.0

def delta_chi2(sf):
    return chi2(sf) - target

sf_low = brentq(delta_chi2, 0.0, sf_best)
sf_high = brentq(delta_chi2, sf_best, 2.0)

sf_err_minus = sf_best - sf_low
sf_err_plus = sf_high - sf_best

print(f"SF = {sf_best:.6f} -{sf_err_minus:.6f} +{sf_err_plus:.6f}")
