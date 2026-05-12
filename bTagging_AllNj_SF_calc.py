import numpy as np
from scipy.optimize import minimize_scalar

# Columns: SV centre, yield, error
# =================================================
# 2022 real data
# =================================================
data_yield = np.array([
    [-0.25, 1172.000000, 34.234486],
    [ 0.25,  165.000000, 12.845233],
    [ 0.75, 1021.000000, 31.953091],
    [ 1.25, 2333.000000, 48.301139],
    [ 1.75, 4057.000000, 63.694584],
    [ 2.25, 5356.000000, 73.184698],
    [ 2.75, 5684.000000, 75.392307],
    [ 3.25, 4705.000000, 68.593003],
    [ 3.75, 3379.000000, 58.129167],
    [ 4.25, 1834.000000, 42.825226],
    [ 4.75,  897.000000, 29.949958],
    [ 5.25,  386.000000, 19.646883],
], dtype=float)

# =================================================
# MC data
# =================================================
b_yield = np.array([
    [-0.25, 1226.807739,  8.499174],
    [ 0.25,  167.533417,  3.145997],
    [ 0.75, 1012.231750,  7.735671],
    [ 1.25, 2440.517578, 12.020829],
    [ 1.75, 4492.460449, 16.329641],
    [ 2.25, 6115.364258, 19.059971],
    [ 2.75, 6570.312012, 19.759647],
    [ 3.25, 5565.834473, 18.206757],
    [ 3.75, 3781.346191, 15.006824],
    [ 4.25, 2141.918457, 11.291224],
    [ 4.75, 1014.655823,  7.788656],
    [ 5.25,  422.212189,  5.028871],
], dtype=float)

nonb_yield = np.array([
    [-0.25, 2.363751, 0.362003],
    [ 0.25, 0.293880, 0.132216],
    [ 0.75, 1.601413, 0.304359],
    [ 1.25, 3.826025, 0.468395],
    [ 1.75, 3.991644, 0.478313],
    [ 2.25, 4.380769, 0.504499],
    [ 2.75, 3.733285, 0.472312],
    [ 3.25, 3.953874, 0.488992],
    [ 3.75, 1.571640, 0.309339],
    [ 4.25, 0.459528, 0.158132],
    [ 4.75, 0.028092, 0.028092],
    [ 5.25, 0.310903, 0.130178],
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

print("Nj >= 2")
print(f"SF = {sf_best:.6f}")
print(f"chi2_min = {chi2_min:.6f}")
print(f"ndof = {ndof}")
print(f"chi2/ndof = {chi2_min / ndof:.6f}")

# =================================================
# Error analysis - statistical since SF not obtained from direct algebraic analysis
# 68.3% CL / 1 sigma statistical uncertainty# https://scikit-hep.org/iminuit/reference.html?utm_source=chatgpt.com
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
