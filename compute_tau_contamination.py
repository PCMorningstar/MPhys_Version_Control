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
# Inputs:
#   - pairing indices per method: {best_i, best_j}
#   - ordered truth flavour per reco jet
# --------------------------------------------------
ordered_truth_flav = tree["ordered_jet_truth_flavour_NOSYS"].array(library="np")

pairing_idx_branches = {
    "chi2":  "raw_chi2_pairing_NOSYS",
    "mdrs":  "raw_mdrs_pairing_NOSYS",
    "misms": "raw_misms_pairing_NOSYS",
}

# sanity: ensure branches exist
keys = set(tree.keys())
for m, b in pairing_idx_branches.items():
    if b not in keys:
        raise RuntimeError(f"Missing branch '{b}' for method '{m}'. Check your ntuple branch names.")

pair_idx = {m: tree[b].array(library="np") for m, b in pairing_idx_branches.items()}

# --------------------------------------------------
# tau-contamination:
#   contam = N(valid_pair & (fi==tau OR fj==tau)) / N(valid_pair)
#
# Weighted:
#   denom = sum_w(valid_pair)
#   num   = sum_w(valid_pair & (fi==tau OR fj==tau))
#   err   = sqrt(sum(w_pass^2)) / denom
# --------------------------------------------------
def weighted_pair_tau_contamination_oriented(pair_indices, truth_flavs, weights, tau_code=15):
    pair_indices = np.asarray(pair_indices, dtype=object)
    truth_flavs  = np.asarray(truth_flavs,  dtype=object)
    weights      = np.asarray(weights,      dtype=float)

    w_den  = 0.0
    w_num  = 0.0
    w2_num = 0.0

    n_valid_events = 0
    example_pairs = []
    example_flavs = []

    for ev_p, ev_f, w in zip(pair_indices, truth_flavs, weights):
        if ev_p is None or ev_f is None:
            continue
        if np.isscalar(ev_p):
            continue
        if len(ev_p) < 2:
            continue

        i, j = int(ev_p[0]), int(ev_p[1])
        if i < 0 or j < 0:
            continue

        n = len(ev_f)
        if i >= n or j >= n:
            continue
        if i == j:
            continue

        # valid reco pair -> denominator
        w_den += w
        n_valid_events += 1

        fi = int(ev_f[i])
        fj = int(ev_f[j])

        if len(example_pairs) < 5:
            example_pairs.append((i, j))
            example_flavs.append((fi, fj))

        # numerator: at least one is tau (code 15)
        if (fi == tau_code) or (fj == tau_code):
            w_num  += w
            w2_num += w * w

    if w_den <= 0.0:
        return 0.0, 0.0, 0.0, 0.0, n_valid_events, example_pairs, example_flavs

    contam = w_num / w_den
    err    = np.sqrt(w2_num) / w_den
    return float(contam), float(err), float(w_num), float(w_den), n_valid_events, example_pairs, example_flavs


def print_block(title, rows):
    print(title)
    print("njets, tau_contam, error")
    for nj, p, e in rows:
        print(f"{nj}, {p:.4f}, {e:.4f}")
    print()


# --------------------------------------------------
# Run tau-contamination per method, per region
# --------------------------------------------------
TAU_CODE = 15  # tau PDG ID

for method in ["chi2", "mdrs", "misms"]:
    rows = []

    for nj in range(2, 11):
        m = regions[nj]

        tcont, e, wnum, wden, n_valid, ex_pairs, ex_flavs = weighted_pair_tau_contamination_oriented(
            pair_idx[method][m],
            ordered_truth_flav[m],
            weights_all[m],
            tau_code=TAU_CODE,
        )

        rows.append((nj, tcont, e))

    print_block(
        f"{method}: tau-contamination = N(valid_pair & (fi=={TAU_CODE} OR fj=={TAU_CODE})) / N(valid_pair)",
        rows
    )