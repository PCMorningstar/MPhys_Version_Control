import uproot
import numpy as np

ntuple_path = "output_ntuples/ttll_601230_mc23a_fullsim.root"
tree = uproot.open(ntuple_path)["reco"]

# --------------------------------------------------
# Event selection + weights
# --------------------------------------------------
selection_cuts = tree["selection_cuts_NOSYS"].array(library="np")
jet_size       = tree["jet_size_NOSYS"].array(library="np")

weights_all = (
    tree["weight_mc_NOSYS"].array(library="np") *
    tree["weight_pileup_NOSYS"].array(library="np") *
    tree["weight_leptonSF_tight_NOSYS"].array(library="np") *
    tree["weight_jvt_effSF_NOSYS"].array(library="np")
)

regions = {nj: (selection_cuts == 1) & (jet_size == nj) for nj in range(2, 11)}

# --------------------------------------------------
# Inputs (chi2 only)
# --------------------------------------------------
ordered_truth_flav = tree["ordered_jet_truth_flavour_NOSYS"].array(library="np")
pair_idx_chi2      = tree["raw_chi2_pairing_NOSYS"].array(library="np")

# --------------------------------------------------
# non-b, non-b (oriented but symmetric here):
#   (fi != b) AND (fj != b)
#
# Denominator = all valid chi2 pairs (after selection & njets region)
# --------------------------------------------------
def weighted_nonb_nonb(pair_indices, truth_flavs, weights, b_code=5):
    pair_indices = np.asarray(pair_indices, dtype=object)
    truth_flavs  = np.asarray(truth_flavs,  dtype=object)
    weights      = np.asarray(weights,      dtype=float)

    w_den  = 0.0
    w_num  = 0.0
    w2_num = 0.0

    n_valid = 0

    for ev_p, ev_f, w in zip(pair_indices, truth_flavs, weights):
        if ev_p is None or ev_f is None:
            continue
        if np.isscalar(ev_p) or len(ev_p) < 2:
            continue

        i, j = int(ev_p[0]), int(ev_p[1])
        if i < 0 or j < 0:
            continue

        n = len(ev_f)
        if i >= n or j >= n or i == j:
            continue

        fi = int(ev_f[i])
        fj = int(ev_f[j])

        # denominator: any valid reco pair
        w_den += w
        n_valid += 1

        # numerator: both are non-b
        if (fi != b_code) and (fj != b_code):
            w_num  += w
            w2_num += w * w

    if w_den <= 0.0:
        return 0.0, 0.0, 0.0, 0

    frac = w_num / w_den
    err  = np.sqrt(w2_num) / w_den
    return float(frac), float(err), float(w_den), int(n_valid)


def print_block(title, rows):
    print(title)
    print("njets, frac_nonb_nonb, err_nonb_nonb")
    for nj, f, e in rows:
        print(f"{nj}, {f:.4f}, {e:.4f}")
    print()


B_CODE = 5

rows = []
for nj in range(2, 11):
    m = regions[nj]
    frac, err, wden, n_valid = weighted_nonb_nonb(
        pair_idx_chi2[m],
        ordered_truth_flav[m],
        weights_all[m],
        b_code=B_CODE,
    )
    rows.append((nj, frac, err))

print_block("chi2: non-b, non-b over valid chi2 pairs", rows)