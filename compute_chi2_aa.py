import ROOT
import numpy as np

FILES = ["output_ntuples/ttll_601230_mc23a_fullsim.root"]
TREE  = "reco"

CHI2_COL = "raw_chi2_minval_truthall_NOSYS"
NJET_COL = "jet_size_NOSYS"
SEL_COL  = "selection_cuts_NOSYS"

JET_BINS = list(range(2, 11))

SELECTION_EXPR = f"({SEL_COL} == 1)"

WEIGHT_EXPR = (
    "weight_mc_NOSYS * "
    "weight_pileup_NOSYS * "
    "weight_leptonSF_tight_NOSYS * "
    "weight_jvt_effSF_NOSYS"
)

# ---------------------------------
# Helper functions
# ---------------------------------

def safe_div(num, den):
    return float(num) / float(den) if den else float("nan")

def njet_to_bin_idx_expr(jet_bins):
    pieces = []
    for ib, n in enumerate(jet_bins):
        pieces.append(f"(njet == {n}) ? {ib} : ")
    return "".join(pieces) + "-1"

def weighted_eff_and_unc(sumw_num, sumw2_num, sumw_den, sumw2_den):

    if sumw_den <= 0:
        return float("nan"), float("nan")

    p = sumw_num / sumw_den

    var = (sumw2_num*(1-2*p) + (p*p)*sumw2_den) / (sumw_den*sumw_den)

    if var < 0 and var > -1e-15:
        var = 0

    return float(p), float(np.sqrt(var))

# ---------------------------------
# Dataframe
# ---------------------------------

df = ROOT.RDataFrame(TREE, FILES)

df = df.Filter(SELECTION_EXPR)

df = df.Define("w", WEIGHT_EXPR)
df = df.Define("w2", "w*w")

df = df.Define("njet", f"static_cast<int>({NJET_COL})")
df = df.Define("njet_bin", njet_to_bin_idx_expr(JET_BINS))
df = df.Filter("njet_bin >= 0")

df = df.Define("best_i",     f"{CHI2_COL}[0]")
df = df.Define("best_j",     f"{CHI2_COL}[1]")
df = df.Define("truth_b",    f"{CHI2_COL}[9]")
df = df.Define("truth_bbar", f"{CHI2_COL}[10]")

df = df.Define("hasTruth", "(truth_b >= 0.0f) && (truth_bbar >= 0.0f)")

# ---------------------------------
# Categories
# ---------------------------------

df = df.Define(
    "is_both_ord",
    "hasTruth && (best_i == truth_b) && (best_j == truth_bbar)"
)

df = df.Define(
    "is_flipped",
    "hasTruth && (best_i == truth_bbar) && (best_j == truth_b)"
)

df = df.Define(
    "picked_truth_b",
    "hasTruth && ((best_i == truth_b) || (best_j == truth_b))"
)

df = df.Define(
    "picked_truth_bbar",
    "hasTruth && ((best_i == truth_bbar) || (best_j == truth_bbar))"
)

df = df.Define(
    "is_f1",
    "hasTruth && (picked_truth_b != picked_truth_bbar)"
)

df = df.Define(
    "is_f0",
    "hasTruth && !(picked_truth_b) && !(picked_truth_bbar)"
)

# ---------------------------------
# Histograms
# ---------------------------------

nbins = len(JET_BINS)

def hist_sum(expr, weight):
    return df.Filter(expr).Histo1D(("h","",nbins,0,nbins),"njet_bin",weight)

h_den_w  = hist_sum("hasTruth", "w")
h_den_w2 = hist_sum("hasTruth", "w2")

h_both_w  = hist_sum("is_both_ord", "w")
h_both_w2 = hist_sum("is_both_ord", "w2")

h_flip_w  = hist_sum("is_flipped", "w")
h_flip_w2 = hist_sum("is_flipped", "w2")

h_f1_w  = hist_sum("is_f1", "w")
h_f1_w2 = hist_sum("is_f1", "w2")

h_f0_w  = hist_sum("is_f0", "w")
h_f0_w2 = hist_sum("is_f0", "w2")

# ---------------------------------
# Arrays
# ---------------------------------

def arr(h):
    return np.array([float(h.GetBinContent(i+1)) for i in range(nbins)])

W_den  = arr(h_den_w)
W2_den = arr(h_den_w2)

W_both  = arr(h_both_w)
W2_both = arr(h_both_w2)

W_flip  = arr(h_flip_w)
W2_flip = arr(h_flip_w2)

W_f1  = arr(h_f1_w)
W2_f1 = arr(h_f1_w2)

W_f0  = arr(h_f0_w)
W2_f0 = arr(h_f0_w2)

# ---------------------------------
# Efficiencies
# ---------------------------------

eps_both, sig_both = [], []
eps_flip, sig_flip = [], []
eps_1, sig_1 = [], []
eps_0, sig_0 = [], []

for i in range(nbins):

    p,s = weighted_eff_and_unc(W_both[i],W2_both[i],W_den[i],W2_den[i])
    eps_both.append(p); sig_both.append(s)

    p,s = weighted_eff_and_unc(W_flip[i],W2_flip[i],W_den[i],W2_den[i])
    eps_flip.append(p); sig_flip.append(s)

    p,s = weighted_eff_and_unc(W_f1[i],W2_f1[i],W_den[i],W2_den[i])
    eps_1.append(p); sig_1.append(s)

    p,s = weighted_eff_and_unc(W_f0[i],W2_f0[i],W_den[i],W2_den[i])
    eps_0.append(p); sig_0.append(s)

# ---------------------------------
# Printing blocks for storage
# ---------------------------------

def print_block(name, vals, unc):

    print(f"\n{name}")
    print("njet value unc")

    for i,n in enumerate(JET_BINS):
        print(f"{n:>4d} {vals[i]:.8f} {unc[i]:.8f}")

print_block("epsilon_both_ord", eps_both, sig_both)
print_block("epsilon_flipped",  eps_flip, sig_flip)
print_block("epsilon_1",        eps_1, sig_1)
print_block("epsilon_0",        eps_0, sig_0)