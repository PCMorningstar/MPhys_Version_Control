ROOT::VecOps::RVec<float> PairLeptonsToJets_Chi2TruthAll(
    const ROOT::VecOps::RVec<TLV>& el_TLV,
    const ROOT::VecOps::RVec<TLV>& mu_TLV,
    const ROOT::VecOps::RVec<TLV>& jet_TLV,
    float met_met,
    float met_phi,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    unsigned long n_electrons,
    unsigned long n_muons
) {
    
    float mean_mlb  = 92000.0f;  // Default
    float sigma_mlb = 21000.0f;

    if (jet_TLV.size() == 2) {
        mean_mlb  = 96374.0f;
        sigma_mlb = 27445.0f;
    }
    else if (jet_TLV.size() == 3) {
        mean_mlb  = 93600.0f;
        sigma_mlb = 21131.0f;
    }
    else if (jet_TLV.size() == 4) {
        mean_mlb  = 93000.0f;
        sigma_mlb = 26653.0f;
    }
    else if (jet_TLV.size() >= 5) {
        mean_mlb  = 90600.0f;
        sigma_mlb = 24341.0f;
    }

    const float mean_mllbb  = 357000.0f;  // MeV
    const float sigma_mllbb = 140000.0f;  // MeV

    const float mean_pTdiff = -363.5f;       // MeV
    const float sigma_pTdiff = 65000.0f;  // MeV

    const float mean_mT     = 389760.0f;  // MeV
    const float sigma_mT    = 88627.1f;  // MeV

    const float mean_deltaR = 3.567f;
    const float sigma_deltaR = 1.211f;

    const float w_mlb   = 1.0f;
    const float w_mllbb = 0.0f;
    const float w_pT    = 1.0f;
    const float w_mT    = 1.0f;
    const float w_deltaR = 1.0f;

    // Output layout:
    // 0: jet index for lepton 0 (positive lepton)
    // 1: jet index for lepton 1 (negative lepton)
    // 2: m(l0, j0)
    // 3: m(l1, j1)
    // 4: chi2 of chosen assignment
    // 5: correctness flag  (-1 no truth, 0 wrong, 1 correct)
    // 6: chi2 of truth (b, bbar) assignment (if available)
    // 7: delta chi2 (chosen - truth)  (if truth available)
    // 8: chi2 of swapped truth (bbar, b) assignment (if available)
    const float SENTINEL_NEG = -1.0f;
    ROOT::VecOps::RVec<float> out(9, SENTINEL_NEG);

    // Require exactly one e and one mu.
    if (!(n_electrons == 1 && n_muons == 1)) return out;
    if (el_TLV.empty() || mu_TLV.empty()) return out;
    if (jet_TLV.size() < 2) return out;

    // Charge vectors must be present; require opposite sign.
    if (el_charge.size() < 1 || mu_charge.size() < 1) return out;
    if (el_charge[0] * mu_charge[0] >= 0.0f) {
        // Same-sign or zero charge; leave sentinels.
        return out;
    }

    // Order leptons so leps[0] is positive (defines our convention).
    ROOT::VecOps::RVec<TLV> leps(2);
    if (el_charge[0] > 0.0f) {
        leps[0] = el_TLV[0];
        leps[1] = mu_TLV[0];
    } else {
        leps[0] = mu_TLV[0];
        leps[1] = el_TLV[0];
    }

    // Precompute MET components (assumed MeV).
    const float met_px = met_met * std::cos(met_phi);
    const float met_py = met_met * std::sin(met_phi);

    struct Chi2Terms {
        float chi2;
        float mlb0;
        float mlb1;
        float mllbb;
        float pTdiff;
        float mT_ttbar;
        float sum_deltaR;
    };

    // Helper to compute chi2 (and keep terms) for assignment (ii, jj).
    auto evaluate_pair = [&](size_t ii, size_t jj) -> Chi2Terms {
        Chi2Terms info{};
        info.mlb0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], jet_TLV[ii]);
        info.mlb1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], jet_TLV[jj]);

        const TLV vis0 = leps[0] + jet_TLV[ii];
        const TLV vis1 = leps[1] + jet_TLV[jj];
        info.pTdiff = vis0.Pt() - vis1.Pt();
        info.sum_deltaR = deltaR(leps[0], jet_TLV[ii]) + deltaR(leps[1], jet_TLV[jj]);

        const TLV vis = leps[0] + leps[1] + jet_TLV[ii] + jet_TLV[jj];
        info.mllbb = vis.M();
        const float m_vis = vis.M();
        const float pT_vis = vis.Pt();
        const float Et_vis = std::sqrt(std::max(0.0f, m_vis * m_vis + pT_vis * pT_vis));
        const float px_sum = vis.Px() + met_px;
        const float py_sum = vis.Py() + met_py;
        const float arg = (Et_vis + met_met) * (Et_vis + met_met) - (px_sum * px_sum + py_sum * py_sum);
        info.mT_ttbar = (arg > 0.0f) ? std::sqrt(arg) : 0.0f;

        const float term_mlb0 = (info.mlb0 - mean_mlb) / sigma_mlb;
        const float term_mlb1 = (info.mlb1 - mean_mlb) / sigma_mlb;
        const float term_mllbb = (info.mllbb - mean_mllbb) / sigma_mllbb;
        const float term_pT = (info.pTdiff - mean_pTdiff) / sigma_pTdiff;
        const float term_mT = (info.mT_ttbar - mean_mT) / sigma_mT;
        const float term_deltaR = (info.sum_deltaR - mean_deltaR) / sigma_deltaR;

        info.chi2 = w_mlb   * (term_mlb0 * term_mlb0 + term_mlb1 * term_mlb1)
                  + w_mllbb * (term_mllbb * term_mllbb)
                  + w_pT    * (term_pT * term_pT)
                  + w_mT    * (term_mT * term_mT)
                  + w_deltaR * (term_deltaR * term_deltaR);
        return info;
    };

    // Scan all ordered jet pairs.
    float best_chi2 = std::numeric_limits<float>::infinity();
    int best_i = -1;
    int best_j = -1;
    Chi2Terms best_terms{};

    const size_t nj = jet_TLV.size();
    for (size_t i = 0; i < nj; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            if (i == j) continue;
            Chi2Terms info = evaluate_pair(i, j);
            if (info.chi2 < best_chi2) {
                best_chi2 = info.chi2;
                best_i = static_cast<int>(i);
                best_j = static_cast<int>(j);
                best_terms = info;
            }
        }
    }

    if (best_i < 0 || best_j < 0) {
        // No valid pair found.
        return out;
    }

    out[0] = static_cast<float>(best_i);
    out[1] = static_cast<float>(best_j);
    out[2] = best_terms.mlb0;
    out[3] = best_terms.mlb1;
    out[4] = best_terms.chi2;

    // --- Truth handling (assumed: event_jet_truth_idx already in the same ordering as jet_TLV) ---
    int truth_b = -1;
    int truth_bbar = -1;
    bool haveTruth = false;

    if (event_jet_truth_idx.size() >= 4 && event_jet_truth_candidates.size() >= 4) {
        const bool validB    = (event_jet_truth_idx[0] != -1 && event_jet_truth_candidates[0] == 1);
        const bool validBbar = (event_jet_truth_idx[3] != -1 && event_jet_truth_candidates[3] == 1);
        if (validB && validBbar) {
            truth_b    = event_jet_truth_idx[0];
            truth_bbar = event_jet_truth_idx[3];
            // NOTE: this assumes truth indices already correspond to positions in jet_TLV.
            if (truth_b >= 0 && truth_bbar >= 0 &&
                truth_b < static_cast<int>(nj) &&
                truth_bbar < static_cast<int>(nj) &&
                truth_b != truth_bbar) {
                haveTruth = true;
            }
        }
    }

    if (!haveTruth) {
        // out[5], out[6], out[7], out[8] remain sentinel -1.
        return out;
    }

    // Compute chi2 for truth orientations.
    const Chi2Terms truth_terms = evaluate_pair(static_cast<size_t>(truth_b),
                                                static_cast<size_t>(truth_bbar));
    const Chi2Terms truth_terms_swapped = evaluate_pair(static_cast<size_t>(truth_bbar),
                                                        static_cast<size_t>(truth_b));

    out[6] = truth_terms.chi2;
    out[7] = best_terms.chi2 - truth_terms.chi2;
    out[8] = truth_terms_swapped.chi2;

    const bool correct =
        (best_i == truth_b && best_j == truth_bbar);
        //(best_i == truth_bbar && best_j == truth_b);
    out[5] = correct ? 1.0f : 0.0f;

    return out;
}
