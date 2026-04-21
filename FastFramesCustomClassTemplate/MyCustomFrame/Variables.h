# pragma once

#include <string>
#include <vector>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <tuple>


using V4 = ROOT::Math::PtEtaPhiEVector;
using ROOT::VecOps::RVec;

namespace ttZ{ 

  struct Lepton {
    ROOT::Math::PtEtaPhiEVector p4;
    float charge;
    Lepton(const ROOT::Math::PtEtaPhiEVector& vec, float q)
      : p4(vec), charge(q) {}
  };
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Jet pT selection (changed in yaml) /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool jet_pt_region_0to30_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_30to60_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_60to90_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_90to120_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_120to150_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_150to180_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_180to210_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_210to240_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_240to270_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_270to300_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_300to360_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);
  bool jet_pt_region_360to900_GeV(const RVec<float>& jet_pt,
    const RVec<float>& chi_pair);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// SV Invariant Mass selection (changed in yaml) /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool sv_invariant_mass_region_neg0point5upto0_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_0to0point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_0point5to1_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_1to1point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_1point5to2_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_2to2point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_2point5to3_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_3to3point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_3point5to4_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_4to4point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_4point5to5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_5to5point5_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);
  bool sv_invariant_mass_region_5point5to6_GeV(const RVec<float>& sv_mass, const RVec<float>& chi_pair);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Order by pT /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  RVec<float> pt_order(const RVec<float>& fv_comp, const RVec<float>& pt);

  RVec<char> pt_order_nonfloat(const RVec<char>& comp, const RVec<float>& pt);

  RVec<int> pt_order_int(const RVec<int>& comp, const RVec<float>& pt);

  int b_selector(const RVec<int>& comp);

  int bbar_selector(const RVec<int>& comp);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int jet_size(const RVec<float>& j_pt);
  int jet_size_two(const RVec<float>& j_pt);
  int jet_size_three(const RVec<float>& j_pt);
  int jet_size_four(const RVec<float>& j_pt);
  int jet_size_five(const RVec<float>& j_pt);
  int jet_size_six(const RVec<float>& j_pt);
  int jet_size_seven(const RVec<float>& j_pt);
  int jet_size_eight(const RVec<float>& j_pt);
  int jet_size_nine(const RVec<float>& j_pt);
  int jet_size_ten(const RVec<float>& j_pt);
  
  // -----------------------------------------------------------------------------
  // A1 — Electron transverse momentum requirement
  // Requirement: pT(e) ≥ 28 GeV
  // -----------------------------------------------------------------------------
  bool cutA1_el_pt(const RVec<float>& el_pt);
  // -----------------------------------------------------------------------------
  // A2 — Electron pseudorapidity acceptance
  // Requirement: |η(e)| < 2.47
  // -----------------------------------------------------------------------------
  bool cutA2_el_eta(const RVec<float>& el_eta);
  // -----------------------------------------------------------------------------
  // A3 — Electron calorimeter crack veto
  // Requirement: Reject 1.37 ≤ |η| < 1.52
  // -----------------------------------------------------------------------------
  bool cutA3_el_crack(const RVec<float>& el_eta);
  // -----------------------------------------------------------------------------
  // A4 — Electron identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA4_el_tight(const RVec<char>& el_tight);
  // -----------------------------------------------------------------------------
  // A5 — Muon transverse momentum requirement
  // Requirement: pT(μ) ≥ 28 GeV
  // -----------------------------------------------------------------------------
  bool cutA5_mu_pt(const RVec<float>& mu_pt);
  // -----------------------------------------------------------------------------
  // A6 — Muon pseudorapidity acceptance
  // Requirement: |η(μ)| < 2.5
  // -----------------------------------------------------------------------------
  bool cutA6_mu_eta(const RVec<float>& mu_eta);
  // -----------------------------------------------------------------------------
  // A7 — Muon identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA7_mu_tight(const RVec<char>& mu_tight);
  // -----------------------------------------------------------------------------
  // A7.1 — Jet identification quality
  // Requirement: select jvt flag == 1 (like the lepton tight_ID)
  // -----------------------------------------------------------------------------
  bool cutA7_1_jet_jvt(const RVec<char>& jet_jvt);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.2 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------
  // C2.1 — Exactly two prompt leptons
  // Requirement: N(leptons) == 2
  // -----------------------------------------------------------------------------
  bool cutC21_exactly2leptons(const RVec<float>& el_pt,
    const RVec<float>& mu_pt);
  // -----------------------------------------------------------------------------
  // C2.2 — Opposite-charge requirement for dileptons
  // Requirement: q(e) * q(μ) == -1
  // -----------------------------------------------------------------------------
  bool cutC22_opposite_charge(const RVec<float>& el_charge,
    const RVec<float>& mu_charge);
  // -----------------------------------------------------------------------------
  // C2.3 — eμ final-state requirement
  // Requirement: Exactly 1 electron and exactly 1 muon
  // -----------------------------------------------------------------------------
  bool cutC23_one_el_one_mu(const RVec<float>& el_pt,
    const RVec<float>& mu_pt);

  // -----------------------------------------------------------------------------
  // C2.4 — Invariant mass threshold m(eμ) ≥ 50 GeV
  // Requirement: m(eμ) ≥ 50 GeV
  // -----------------------------------------------------------------------------
  bool cutC24_meemu_gt50(const RVec<float>& el_pt,  const RVec<float>& mu_pt,
    const RVec<float>& el_eta, const RVec<float>& mu_eta,
    const RVec<float>& el_phi, const RVec<float>& mu_phi);
  // -----------------------------------------------------------------------------
  // C2.5 — Exactly two reconstructed jets
  // Requirement: N(jets) == 2
  // -----------------------------------------------------------------------------
  bool cutC25_exactly2jets(const RVec<float>& j_pt);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.3 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct PairingResult {
    bool valid = false;
    float mj1 = -1.f;
    float mj2 = -1.f;
  };
  // -----------------------------------------------------------------------------
  // Compute Section 6.3 minimal pairing
  //
  // Purpose:
  //   Given 1 electron, 1 muon, and 2 jets (j1 = leading, j2 = subleading),
  //   compute the two possible top-decay hypotheses:
  //
  //      Pairing A: (e ↔ j1) and (μ ↔ j2)
  //      Pairing B: (e ↔ j2) and (μ ↔ j1)
  //
  //   Select the pairing that MINIMISES:
  //        m(j1,l_i)^2 + m(j2,l_j)^2
  //
  // Returned values:
  //   - mj1, mj2 : invariant masses of the chosen pairing
  //   - valid    : false if jets/leptons missing
  // -----------------------------------------------------------------------------
  PairingResult compute_pairing_S6_3(
    const RVec<float>& el_pt,  const RVec<float>& el_eta,
    const RVec<float>& el_phi, const RVec<float>& el_e,
    const RVec<float>& mu_pt,  const RVec<float>& mu_eta,
    const RVec<float>& mu_phi, const RVec<float>& mu_e,
    const RVec<float>& j_pt,   const RVec<float>& j_eta,
    const RVec<float>& j_phi,  const RVec<float>& j_e);
  // -----------------------------------------------------------------------------
  // D3.1 — Valid pairing computed
  // Requirement: compute_pairing_S6_3() returned valid = true
  // -----------------------------------------------------------------------------
  bool cutD31_pairing_valid(const PairingResult& P);
  // -----------------------------------------------------------------------------
  // D3.2 — Leading pairing mass threshold
  // Requirement: mj1 ≥ 20 GeV
  // -----------------------------------------------------------------------------
  bool cutD32_mj1(const PairingResult& P);
  // -----------------------------------------------------------------------------
  // D3.3 — Subleading pairing mass threshold
  // Requirement: mj2 ≥ 20 GeV
  // -----------------------------------------------------------------------------
  bool cutD33_mj2(const PairingResult& P);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  RVec<int> jets_clean_from_e(
      const RVec<float>& jet_eta,
      const RVec<float>& jet_phi,
      const RVec<float>& el_eta,
      const RVec<float>& el_phi
  );

  int electron_clean_from_jets(
      const RVec<float>& el_eta,
      const RVec<float>& el_phi,
      const RVec<float>& jet_eta,
      const RVec<float>& jet_phi
  );

  int muon_clean_from_jets(
      const RVec<float>& mu_eta,
      const RVec<float>& mu_phi,
      const RVec<float>& jet_eta,
      const RVec<float>& jet_phi
  );

  bool cutA8_jet_clean(const RVec<int>& jet_keep_flags);
  bool cutA9_el_clean(const int& el_keep_flag);
  bool cutA10_mu_clean(const int& mu_keep_flag);

  
  bool section_6_1(
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<char>&  el_tight,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<char>&  mu_tight,
    const RVec<char>& jet_jvt
  );
  bool section_6_2(
    const RVec<float>& el_pt,
    const RVec<float>& mu_pt,
    const RVec<float>& el_eta,
    const RVec<float>& mu_eta,
    const RVec<float>& el_phi,
    const RVec<float>& mu_phi,
    const RVec<float>& el_charge,
    const RVec<float>& mu_charge,
    const RVec<float>& j_pt
  );
  bool section_6_3(
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi,
    const RVec<float>& el_e,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& mu_e,
    const RVec<float>& j_pt,
    const RVec<float>& j_eta,
    const RVec<float>& j_phi,
    const RVec<float>& j_e
  );

  bool selection_cuts(
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi,
    const RVec<float>& el_e,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& mu_e,
    const RVec<float>& j_pt,
    const RVec<float>& j_eta,
    const RVec<float>& j_phi,
    const RVec<float>& j_e,
    const RVec<float>& el_charge,
    const RVec<float>& mu_charge,
    const RVec<int>& jet_keep_flags,
    int& el_keep_flag,
    int& mu_keep_flag,
    const RVec<char>& el_tight,
    const RVec<char>& mu_tight,
    const RVec<char>& jet_jvt
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Selections for cut-flow /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool electron_selections_paper(
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<char>&  el_tight,
    const int& el_keep_flag
  );
  bool muon_selections_paper(
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<char>&  mu_tight,
    const int& mu_keep_flag
  );
  bool jet_selections_paper(
    const RVec<char>& jet_jvt,
    const RVec<int>& jet_keep_flags,
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi,
    const RVec<float>& el_e,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& mu_e,
    const RVec<float>& j_pt,
    const RVec<float>& j_eta,
    const RVec<float>& j_phi,
    const RVec<float>& j_e
  );
  bool dilepton_selections_paper(
    const RVec<float>& el_pt,
    const RVec<float>& mu_pt,
    const RVec<float>& el_eta,
    const RVec<float>& mu_eta,
    const RVec<float>& el_phi,
    const RVec<float>& mu_phi,
    const RVec<float>& el_charge,
    const RVec<float>& mu_charge
  );

  // New attempt at detailed truths - hence new mus and sigmas needed for chi2!
  // ============================================================
  // Detailed truth observables using fixed truth jet indices
  // ============================================================
  RVec<float> new_detailed_truth(
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
    const RVec<float>& mu_charge,
    const float& met_met,
    const float& met_phi,
    const int& event_jet_truth_idx_b,
    const int& event_jet_truth_idx_bbar);

  // Safe scalar extractors - additional for chi2
  float new_truth_mlpb (const RVec<float>& v);
  float new_truth_mlmbb(const RVec<float>& v);
  float new_truth_pTdiff (const RVec<float>& v);
  float new_truth_sum_deltaR(const RVec<float>& v);
  float new_truth_mllbb (const RVec<float>& v);
  float new_truth_mT_ttbar(const RVec<float>& v);

  // ============================================================
  // new chi2 pairing in the (l+, l-) basis
  // returns:
  //  1  -> correct pairing
  //  0  -> wrong pairing
  // -1  -> no valid truth
  // ============================================================
  int new_chi_indexed(
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
    const RVec<float>& mu_charge,
    const float& met_met,
    const float& met_phi,
    // --- truth inputs (pT ordered)
    const int& event_jet_truth_idx_b,
    const int& event_jet_truth_idx_bbar,
    const int& event_jet_truth_candidates_b,
    const int& event_jet_truth_candidates_bbar);

  // ============================================================
  // NEW Minimum Invariant Squared Mass Sum (MISMS)
  // pairing in the (l+, l-) basis
  // Minimises:  M(l+ b)^2 + M(l- b)^2
  // ============================================================

  int new_misms_pairing(
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
    const RVec<float>& mu_charge,
    // --- truth inputs (pT ordered)
    const int& event_jet_truth_idx_b,
    const int& event_jet_truth_idx_bbar,
    const int& event_jet_truth_candidates_b,
    const int& event_jet_truth_candidates_bbar);

  // ============================================================
  // NEW Minimum Delta R Sum (MDRS)
  // pairing in the (l+, l-) basis
  // Minimises:  dR(l+, j1) + dR(l-, j2)
  // ============================================================
  int new_MDRS_pairing(
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
    const RVec<float>& mu_charge,
    const int& event_jet_truth_idx_b,
    const int& event_jet_truth_idx_bbar,
    const int& event_jet_truth_candidates_b,
    const int& event_jet_truth_candidates_bbar
  );

  // ============================================================
  // Raw pairing algorithms in the (l+, l-) basis
  // Returns:
  //  best_i  -> index of the best jet for l+ (for b)
  //  best_j  -> index of the best jet for l- (for bbar)
  // -1      -> no valid pairing
  // ============================================================
  RVec<int> raw_chi2_pairing(
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
  );

  // ============================================================
  // Raw MDRS pairing in the (l+, l-) basis
  // Returns:
  //  best_i  -> index of the best jet for l+ (for b)
  //  best_j  -> index of the best jet for l- (for bbar)
  // -1      -> no valid pairing
  // ============================================================
  RVec<int> raw_MDRS_pairing(
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
  );

  // ============================================================
  // Raw MISMS pairing in the (l+, l-) basis
  // Returns:
  //  best_i  -> index of the best jet for l+ (for b)
  //  best_j  -> index of the best jet for l- (for bbar)
  // -1      -> no valid pairing
  // ============================================================
  RVec<int> raw_MISMS_pairing(
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
  );

  // ============================================================
  // Chi2 - return minimum chi2 only
  // ============================================================
  RVec<float> raw_chi2_minval_truthall(
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
    const RVec<float>& mu_charge,
    const RVec<int>& event_jet_truth_idx,
    const RVec<int>& event_jet_truth_candidates
  );

}