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
# Flavour composition of valid chi2-selected jet pairs
#
# These are three separate flavour studies:
#
# b-study:
#   both_b, one_b_only, nonb_nonb
#
# c-study:
#   both_c, one_c_only, nonc_nonc
#
# light-study:
#   both_light, one_light_only, nonlight_nonlight
#
# Denominator = all valid chi2 pairs (after selection & njets region)
# --------------------------------------------------
def weighted_pair_flavour_fractions(pair_indices, truth_flavs, weights,
                                    b_code=5, c_code=4, l_code=0):
    pair_indices = np.asarray(pair_indices, dtype=object)
    truth_flavs  = np.asarray(truth_flavs, dtype=object)
    weights      = np.asarray(weights, dtype=float)

    w_den = 0.0
    n_valid = 0

    w_both_b = 0.0
    w2_both_b = 0.0
    w_one_b = 0.0
    w2_one_b = 0.0
    w_nonb_nonb = 0.0
    w2_nonb_nonb = 0.0

    w_both_c = 0.0
    w2_both_c = 0.0
    w_one_c = 0.0
    w2_one_c = 0.0
    w_nonc_nonc = 0.0
    w2_nonc_nonc = 0.0

    w_both_light = 0.0
    w2_both_light = 0.0
    w_one_light = 0.0
    w2_one_light = 0.0
    w_nonlight_nonlight = 0.0
    w2_nonlight_nonlight = 0.0

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

        # ---------- b study ----------
        is_b_i = (fi == b_code)
        is_b_j = (fj == b_code)

        if is_b_i and is_b_j:
            w_both_b += w
            w2_both_b += w * w
        elif is_b_i ^ is_b_j:
            w_one_b += w
            w2_one_b += w * w
        else:
            w_nonb_nonb += w
            w2_nonb_nonb += w * w

        # ---------- c study ----------
        is_c_i = (fi == c_code)
        is_c_j = (fj == c_code)

        if is_c_i and is_c_j:
            w_both_c += w
            w2_both_c += w * w
        elif is_c_i ^ is_c_j:
            w_one_c += w
            w2_one_c += w * w
        else:
            w_nonc_nonc += w
            w2_nonc_nonc += w * w

        # ---------- light study ----------
        is_l_i = (fi == l_code)
        is_l_j = (fj == l_code)

        if is_l_i and is_l_j:
            w_both_light += w
            w2_both_light += w * w
        elif is_l_i ^ is_l_j:
            w_one_light += w
            w2_one_light += w * w
        else:
            w_nonlight_nonlight += w
            w2_nonlight_nonlight += w * w

    if w_den <= 0.0:
        return {
            "frac_both_b": 0.0,
            "err_both_b": 0.0,
            "frac_one_b_only": 0.0,
            "err_one_b_only": 0.0,
            "frac_nonb_nonb": 0.0,
            "err_nonb_nonb": 0.0,

            "frac_both_c": 0.0,
            "err_both_c": 0.0,
            "frac_one_c_only": 0.0,
            "err_one_c_only": 0.0,
            "frac_nonc_nonc": 0.0,
            "err_nonc_nonc": 0.0,

            "frac_both_light": 0.0,
            "err_both_light": 0.0,
            "frac_one_light_only": 0.0,
            "err_one_light_only": 0.0,
            "frac_nonlight_nonlight": 0.0,
            "err_nonlight_nonlight": 0.0,

            "w_den": 0.0,
            "n_valid": 0,
        }

    return {
        "frac_both_b": float(w_both_b / w_den),
        "err_both_b": float(np.sqrt(w2_both_b) / w_den),
        "frac_one_b_only": float(w_one_b / w_den),
        "err_one_b_only": float(np.sqrt(w2_one_b) / w_den),
        "frac_nonb_nonb": float(w_nonb_nonb / w_den),
        "err_nonb_nonb": float(np.sqrt(w2_nonb_nonb) / w_den),

        "frac_both_c": float(w_both_c / w_den),
        "err_both_c": float(np.sqrt(w2_both_c) / w_den),
        "frac_one_c_only": float(w_one_c / w_den),
        "err_one_c_only": float(np.sqrt(w2_one_c) / w_den),
        "frac_nonc_nonc": float(w_nonc_nonc / w_den),
        "err_nonc_nonc": float(np.sqrt(w2_nonc_nonc) / w_den),

        "frac_both_light": float(w_both_light / w_den),
        "err_both_light": float(np.sqrt(w2_both_light) / w_den),
        "frac_one_light_only": float(w_one_light / w_den),
        "err_one_light_only": float(np.sqrt(w2_one_light) / w_den),
        "frac_nonlight_nonlight": float(w_nonlight_nonlight / w_den),
        "err_nonlight_nonlight": float(np.sqrt(w2_nonlight_nonlight) / w_den),

        "w_den": float(w_den),
        "n_valid": int(n_valid),
    }


def print_block(title, rows, frac_key, err_key, header):
    print(title)
    print(header)
    for nj, res in rows:
        print(f"{nj}, {res[frac_key]:.4f}, {res[err_key]:.4f}")
    print()


C_CODE = 4
L_CODE = 0
B_CODE = 5

rows = []
for nj in range(2, 11):
    m = regions[nj]
    res = weighted_pair_flavour_fractions(
        pair_idx_chi2[m],
        ordered_truth_flav[m],
        weights_all[m],
        b_code=B_CODE,
        c_code=C_CODE,
        l_code=L_CODE,
    )
    rows.append((nj, res))

# ---------- b study ----------
print_block(
    "chi2: both b over valid chi2 pairs",
    rows,
    "frac_both_b",
    "err_both_b",
    "njets, frac_both_b, err_both_b",
)

print_block(
    "chi2: one b only over valid chi2 pairs",
    rows,
    "frac_one_b_only",
    "err_one_b_only",
    "njets, frac_one_b_only, err_one_b_only",
)

print_block(
    "chi2: non-b, non-b over valid chi2 pairs",
    rows,
    "frac_nonb_nonb",
    "err_nonb_nonb",
    "njets, frac_nonb_nonb, err_nonb_nonb",
)

# ---------- c study ----------
print_block(
    "chi2: both c over valid chi2 pairs",
    rows,
    "frac_both_c",
    "err_both_c",
    "njets, frac_both_c, err_both_c",
)

print_block(
    "chi2: one c only over valid chi2 pairs",
    rows,
    "frac_one_c_only",
    "err_one_c_only",
    "njets, frac_one_c_only, err_one_c_only",
)

print_block(
    "chi2: non-c, non-c over valid chi2 pairs",
    rows,
    "frac_nonc_nonc",
    "err_nonc_nonc",
    "njets, frac_nonc_nonc, err_nonc_nonc",
)

# ---------- light study ----------
print_block(
    "chi2: both light over valid chi2 pairs",
    rows,
    "frac_both_light",
    "err_both_light",
    "njets, frac_both_light, err_both_light",
)

print_block(
    "chi2: one light only over valid chi2 pairs",
    rows,
    "frac_one_light_only",
    "err_one_light_only",
    "njets, frac_one_light_only, err_one_light_only",
)

print_block(
    "chi2: non-light, non-light over valid chi2 pairs",
    rows,
    "frac_nonlight_nonlight",
    "err_nonlight_nonlight",
    "njets, frac_nonlight_nonlight, err_nonlight_nonlight",
)