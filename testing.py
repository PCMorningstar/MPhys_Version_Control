import numpy as np
from scipy.optimize import minimize_scalar

# ================================================= 
# b - Events
# ================================================= 
b_yield = np.array([
    [-0.25,   59613.984375,   2300.355957],
    [ 0.25,   38874.960938,   1847.597656],
    [ 0.75,  389292.562500,   5906.194824],
    [ 1.25, 1463274.500000,  11407.054688],
    [ 1.75, 3339577.500000,  17253.572266],
    [ 2.25, 5080530.000000,  21285.535156],
    [ 2.75, 5715139.000000,  22584.531250],
    [ 3.25, 4958822.000000,  21053.396484],
    [ 3.75, 3384584.250000,  17387.451172],
    [ 4.25, 1927614.250000,  13121.601562],
    [ 4.75,  909799.500000,   9031.028320],
    [ 5.25,  374631.187500,   5805.943848],
], dtype=float)


# ================================================= 
# non-b - Events
# ================================================= 
nonb_yield = np.array([
    [-0.25, 194.777802, 118.028427],
    [ 0.25,  53.177876,  53.177876],
    [ 0.75, 280.261353, 162.934464],
    [ 1.25, 890.433167, 273.773407],
    [ 1.75,1166.223755, 311.399109],
    [ 2.25,2065.873535, 428.523071],
    [ 2.75,1253.535034, 338.515320],
    [ 3.25,2162.624756, 443.213593],
    [ 3.75, 874.067261, 294.135834],
    [ 4.25, 244.068115, 142.938339],
    [ 4.75,   0.000000,   0.000000],   # undefined ratio
    [ 5.25,  79.905853,  79.905853],
], dtype=float)

# ================================================= 
# 2022 real data (scaled to MC normalisation)
# Columns: SV centre, N_data_scaled, err_scaled
# ================================================= 
data_yield_scaled = np.array([
    [-0.25,    54659.784816,   9662.576125],
    [ 0.25,    42702.956888,   8540.591378],
    [ 0.75,   406532.149569,  26351.565309],
    [ 1.25,  1549263.275880,  51442.442740],
    [ 1.75,  3440150.206860,  76656.268099],
    [ 2.25,  5071403.159964,  93072.855441],
    [ 2.75,  5614584.771574,  97930.459295],
    [ 3.25,  4789563.644507,  90449.661099],
    [ 3.75,  3441858.325136,  76675.296588],
    [ 4.25,  1907968.113735,  57087.960238],
    [ 4.75,   932632.578424,  39912.989759],
    [ 5.25,   399699.676467,  26129.185255],
], dtype=float)



N_data = data_yield_scaled[:, 1]
err_data = data_yield_scaled[:, 2]

N_b_MC = b_yield[:, 1]
err_b_MC = b_yield[:, 2]

N_nonb_MC = nonb_yield[:, 1]
err_nonb_MC = nonb_yield[:, 2]

def chi2(sf):
    model = sf * N_b_MC + N_nonb_MC

    sigma2 = (
        err_data**2
        + (sf * err_b_MC)**2
        + err_nonb_MC**2
    )

    return np.sum((N_data - model)**2 / sigma2)

result = minimize_scalar(
    chi2,
    bounds=(0.5, 1.5),
    method="bounded"
)

SF_best = result.x
chi2_min = result.fun
ndof = len(N_data) - 1

print(f"Best-fit SF = {SF_best:.6f}")
print(f"chi2_min    = {chi2_min:.6f}")
print(f"ndof        = {ndof}")
print(f"chi2/ndof   = {chi2_min / ndof:.6f}")

sf_grid = np.linspace(0.5, 1.5, 10000)
chi2_grid = np.array([chi2(sf) for sf in sf_grid])

chi2_min = np.min(chi2_grid)
sf_best = sf_grid[np.argmin(chi2_grid)]

mask_1sigma = chi2_grid <= chi2_min + 3.84

sf_low = sf_grid[mask_1sigma][0]
sf_high = sf_grid[mask_1sigma][-1]

err_low = sf_best - sf_low
err_high = sf_high - sf_best

print(f"SF = {sf_best:.6f} +{err_high:.6f} -{err_low:.6f}")