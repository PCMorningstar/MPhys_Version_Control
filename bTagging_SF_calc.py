import numpy as np
from scipy.optimize import minimize_scalar

# Columns: SV centre, yield, error
data_yield = np.array([
    [-0.25,   32.000000,   5.656854],
    [ 0.25,   26.000000,   5.099020],
    [ 0.75,  243.000000,  15.588457],
    [ 1.25,  936.000000,  30.594117],
    [ 1.75, 2064.000000,  45.431267],
    [ 2.25, 3056.000000,  55.281100],
    [ 2.75, 3388.000000,  58.206529],
    [ 3.25, 2898.000000,  53.833075],
    [ 3.75, 2070.000000,  45.497253],
    [ 4.25, 1139.000000,  33.749074],
    [ 4.75,  560.000000,  23.664319],
    [ 5.25,  236.000000,  15.362291],
], dtype=float)

b_yield = np.array([
    [-0.25,   40.829327,   1.556309],
    [ 0.25,   26.749590,   1.253481],
    [ 0.75,  268.019135,   4.004860],
    [ 1.25, 1010.403015,   7.744958],
    [ 1.75, 2304.259277,  11.712639],
    [ 2.25, 3504.570312,  14.448512],
    [ 2.75, 3945.906494,  15.336808],
    [ 3.25, 3418.974854,  14.287642],
    [ 3.75, 2336.287598,  11.808280],
    [ 4.25, 1326.139038,   8.893599],
    [ 4.75,  626.579956,   6.126587],
    [ 5.25,  257.413513,   3.932517],
], dtype=float)

nonb_yield = np.array([
    [-0.25, 0.130102, 0.078837],
    [ 0.25, 0.035520, 0.035520],
    [ 0.75, 0.187201, 0.108833],
    [ 1.25, 0.594767, 0.182868],
    [ 1.75, 0.778982, 0.208000],
    [ 2.25, 1.379906, 0.286233],
    [ 2.75, 0.837302, 0.226112],
    [ 3.25, 1.499193, 0.301050],
    [ 3.75, 0.583836, 0.196469],
    [ 4.25, 0.163026, 0.095476],
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
