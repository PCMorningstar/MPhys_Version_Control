import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import os

# -------------------------------------------------
# INPUT
# -------------------------------------------------
ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree = uproot.open(ntuple_path)["reco"]

selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
jet_size       = tree["jet_size_NOSYS"].array(library="np")

weights_all = (
    tree["weight_mc_NOSYS"].array(library="np") *
    tree["weight_pileup_NOSYS"].array(library="np") *
    tree["weight_leptonSF_tight_NOSYS"].array(library="np") *
    tree["weight_jvt_effSF_NOSYS"].array(library="np")
)

regions = {
    "2jets_region":  (selection_cuts == 1) & (jet_size == 2),
    "4jets_region":  (selection_cuts == 1) & (jet_size == 4),
    "6jets_region":  (selection_cuts == 1) & (jet_size == 6)
}

# -------------------------------------------------
# CHI2 OUTPUT VEC
#   [0]=best_i, [1]=best_j, [4]=chi2_min, [5]=is_correct (ordered, -1/0/1)
# -------------------------------------------------
chi2vec = tree["raw_chi2_minval_truthall_NOSYS"].array(library="ak")

best_i     = ak.to_numpy(chi2vec[:, 0]).astype(np.int32)
best_j     = ak.to_numpy(chi2vec[:, 1]).astype(np.int32)
chi2_min   = ak.to_numpy(chi2vec[:, 4])
is_correct = ak.to_numpy(chi2vec[:, 5])  # -1 no truth, 0 wrong, 1 correct (ordered)

# -------------------------------------------------
# TRUTH (needed to split wrong into subclasses)
# same convention as your C++: truth_b = idx[0], truth_bbar = idx[3] with candidates==1
# -------------------------------------------------
truth_idx = tree["event_jet_truth_idx"].array(library="ak")
truth_cand = tree["event_jet_truth_candidates"].array(library="ak")

truth_b    = ak.to_numpy(truth_idx[:, 0]).astype(np.int32)
truth_bbar = ak.to_numpy(truth_idx[:, 3]).astype(np.int32)
cand_b     = ak.to_numpy(truth_cand[:, 0]).astype(np.int32)
cand_bbar  = ak.to_numpy(truth_cand[:, 3]).astype(np.int32)

# truth availability mask (must match what C++ assumes)
have_truth = (
    (cand_b == 1) & (cand_bbar == 1) &
    (truth_b >= 0) & (truth_bbar >= 0) &
    (truth_b != truth_bbar)
)

# -------------------------------------------------
# GLOBAL VALIDITY
# -------------------------------------------------
valid = (
    np.isfinite(chi2_min) &
    (chi2_min >= 0) &
    np.isfinite(weights_all) &
    have_truth &
    (best_i >= 0) & (best_j >= 0)
)

# -------------------------------------------------
# Plot helper: step + errorbars, density-normalised
# -------------------------------------------------
def plot_hist_with_errors(x, w, bins, color, label):
    if x.size == 0:
        return

    counts, edges = np.histogram(x, bins=bins, weights=w)
    widths = np.diff(edges)
    centers = edges[:-1] + widths / 2

    total_w = np.sum(counts)
    if total_w <= 0:
        return

    density = counts / (total_w * widths)

    # weighted variance per bin: sum(w^2) gives Poisson-like uncertainty for weighted hist
    var = np.histogram(x, bins=bins, weights=w**2)[0]
    err = np.sqrt(var) / (total_w * widths)

    plt.step(centers, density, where="mid", color=color, linewidth=1.5, label=label)

# -------------------------------------------------
# PLOTTING
# -------------------------------------------------
output_dir = "chi2_plots_matching_v_mismatching"
os.makedirs(output_dir, exist_ok=True)

NBINS = 40
XMIN, XMAX = 0.0, 50.0
bins = np.linspace(XMIN, XMAX, NBINS + 1)

for region_name, region_mask in regions.items():
    m = valid & region_mask

    x = chi2_min[m]
    w = weights_all[m]

    bi = best_i[m]
    bj = best_j[m]
    tb = truth_b[m]
    tbb = truth_bbar[m]

    # Diagnostic categories:
    # A) ordered correct
    mask_correct = (bi == tb) & (bj == tbb)

    # B) swapped truth (both truth jets, reversed)
    mask_swapped = (bi == tbb) & (bj == tb)

    # C) one-of-two correct (exactly one matches its ordered target)
    mask_one_correct = ((bi == tb) ^ (bj == tbb))

    # D) none correct (neither matches its ordered target) AND not swapped
    mask_none_correct = ~(mask_correct | mask_swapped | mask_one_correct)

    if x.size == 0:
        print(f"[skip] {region_name}: no events after cuts")
        continue

    plt.figure(figsize=(8, 6))

    # Colors: requested diagnostic emphasis with purple + orange
    plot_hist_with_errors(x[mask_correct],      w[mask_correct],      bins, "blue",   "Correct")
    plot_hist_with_errors(x[mask_one_correct],  w[mask_one_correct],  bins, "purple", "One Correct")
    plot_hist_with_errors(x[mask_swapped],      w[mask_swapped],      bins, "red", "Swapped Truth")
    plot_hist_with_errors(x[mask_none_correct], w[mask_none_correct], bins, "orange",  "None Correct")

    plt.xlabel(r"$\chi^2_{\min}$")
    plt.ylabel("Probability Density")
    plt.title(f"Chi2 Reconstruction Diagnostics: Probability Density - {region_name}")

    plt.grid(True, linestyle=":", linewidth=0.7)
    plt.legend(loc="upper right", frameon=False)
    plt.xlim(XMIN, XMAX)

    plt.tight_layout()
    outfile = os.path.join(output_dir, f"chi2_pdf_diagnostic_{region_name}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"[saved] {outfile}")