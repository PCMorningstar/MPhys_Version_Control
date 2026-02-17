// ============================================================
// NEW Minimum Invariant Squared Mass Sum (MISMS)
// pairing in the (l+, l-) basis
// Minimises:  M(l+ b)^2 + M(l- b)^2
// ============================================================

RVec<int> new_misms_pairing(
    const RVec<float>& jet_pt,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi,
    const RVec<float>& jet_e,
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi,
    const RVec<float>& el_e,
    const RVec<float>& el_charge,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& mu_e,
    const RVec<float>& mu_charge
  ){
  
    constexpr float GeV = 1.f/1000.f;
  
    // -------------------------
    // Collect leptons
    // -------------------------
    RVec<Lepton> leptons;
  
    for (size_t i=0; i<el_pt.size(); ++i)
      leptons.push_back(
        Lepton{ V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV),
                el_charge[i] });
  
    for (size_t i=0; i<mu_pt.size(); ++i)
      leptons.push_back(
        Lepton{ V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV),
                mu_charge[i] });
  
    if (leptons.size() < 2) return idx;
  
    // -------------------------
    // Pick one l+ and one l−
    // -------------------------
    const V4* lplus  = nullptr;
    const V4* lminus = nullptr;
  
    for (const auto& lep : leptons){
      if (lep.charge > 0.f && !lplus)  lplus  = &lep.p4;
      if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
  
    if (!lplus || !lminus) return idx;
  
    // -------------------------
    // Build jets
    // -------------------------
    std::vector<V4> jets;
  
    for (size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(
        jet_pt[i]*GeV,
        jet_eta[i],
        jet_phi[i],
        jet_e[i]*GeV
      );
  
    if (jets.size() < 2) return idx;
  
    // -------------------------
    // MISMS minimisation
    // -------------------------
    float min_sum = 1e12f;
    int best_i = -1;
    int best_j = -1;
  
    for (size_t i=0; i<jets.size(); ++i){
      for (size_t j=i+1; j<jets.size(); ++j){
  
        // Assignment A: l+ → i, l− → j
        float sumA =
        std::pow((*lplus  + jets[i]).M2(), 2) +
        std::pow((*lminus + jets[j]).M2(), 2);
  
        // Assignment B: l+ → j, l− → i
        float sumB =
          std::pow((*lplus  + jets[j]).M2(), 2) +
          std::pow((*lminus + jets[i]).M2(), 2);
  
        if (sumA < min_sum){
          min_sum = sumA;
          best_i  = i;
          best_j  = j;
        }
  
        if (sumB < min_sum){
          min_sum = sumB;
          best_i  = j;
          best_j  = i;
        }
      }
    }
    // =========================
    // Truth extraction
    // =========================
    int truth_b = -1;
    int truth_bbar = -1;

    if (event_jet_truth_idx_b.size() > 0 &&
        event_jet_truth_candidates_b.size() > 0 &&
        event_jet_truth_idx_bbar.size() > 0 &&
        event_jet_truth_candidates_bbar.size() > 0) {

        bool validB =
            (event_jet_truth_idx_b != -1 &&
                event_jet_truth_candidates_b == 1);

        bool validBbar =
            (event_jet_truth_idx_bbar != -1 &&
                event_jet_truth_candidates_bbar == 1);

        if (validB && validBbar) {
            truth_b    = event_jet_truth_idx_b;
            truth_bbar = event_jet_truth_idx_bbar;
        }
    }
    if (truth_b < 0 || truth_bbar < 0) return -1;

    // =========================
    // Correctness flag (out[5])
    // =========================
    bool correct =
    (best_i == truth_b && best_j == truth_bbar);

    return correct ? 1 : 0;
}