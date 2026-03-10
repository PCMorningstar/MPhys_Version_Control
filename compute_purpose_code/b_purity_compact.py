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
# Helper
# --------------------------------------------------
def flavour_label(f, b_code=5, c_code=4, l_code=0):
    if f == b_code:
        return "b"
    elif f == c_code:
        return "c"
    elif f == l_code:
        return "l"
    else:
        return None  # ignore everything outside b/c/light


def weighted_pair_composition(pair_indices, truth_flavs, weights,
                              b_code=5, c_code=4, l_code=0):
    pair_indices = np.asarray(pair_indices, dtype=object)
    truth_flavs  = np.asarray(truth_flavs, dtype=object)
    weights      = np.asarray(weights, dtype=float)

    cats = ["bb", "bc", "bl", "cc", "cl", "ll"]

    w_num  = {k: 0.0 for k in cats}
    w2_num = {k: 0.0 for k in cats}

    w_den = 0.0
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

        li = flavour_label(fi, b_code=b_code, c_code=c_code, l_code=l_code)
        lj = flavour_label(fj, b_code=b_code, c_code=c_code, l_code=l_code)

        # keep only events where BOTH selected jets are in {b,c,l}
        if li is None or lj is None:
            continue

        pair = "".join(sorted([li, lj]))   # order-independent
        if pair not in w_num:
            continue

        w_den += w
        n_valid += 1

        w_num[pair]  += w
        w2_num[pair] += w * w

    out = {"w_den": float(w_den), "n_valid": int(n_valid)}

    for k in cats:
        if w_den > 0.0:
            out[f"frac_{k}"] = float(w_num[k] / w_den)
            out[f"err_{k}"]  = float(np.sqrt(w2_num[k]) / w_den)
        else:
            out[f"frac_{k}"] = 0.0
            out[f"err_{k}"]  = 0.0

    return out


def print_block(title, rows, frac_key, err_key, header):
    print(title)
    print(header)
    for nj, res in rows:
        print(f"{nj}, {res[frac_key]:.6f}, {res[err_key]:.6f}")
    print()


B_CODE = 5
C_CODE = 4
L_CODE = 0

rows = []
for nj in range(2, 11):
    m = regions[nj]
    res = weighted_pair_composition(
        pair_idx_chi2[m],
        ordered_truth_flav[m],
        weights_all[m],
        b_code=B_CODE,
        c_code=C_CODE,
        l_code=L_CODE,
    )
    rows.append((nj, res))

print_block("chi2: bb over valid chi2 pairs", rows, "frac_bb", "err_bb", "njets, frac_bb, err_bb")
print_block("chi2: bc over valid chi2 pairs", rows, "frac_bc", "err_bc", "njets, frac_bc, err_bc")
print_block("chi2: bl over valid chi2 pairs", rows, "frac_bl", "err_bl", "njets, frac_bl, err_bl")
print_block("chi2: cc over valid chi2 pairs", rows, "frac_cc", "err_cc", "njets, frac_cc, err_cc")
print_block("chi2: cl over valid chi2 pairs", rows, "frac_cl", "err_cl", "njets, frac_cl, err_cl")
print_block("chi2: ll over valid chi2 pairs", rows, "frac_ll", "err_ll", "njets, frac_ll, err_ll")