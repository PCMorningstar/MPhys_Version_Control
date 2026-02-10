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

    // creating multiply dRs for efficiency v. dR
    float dr_one();
    float dr_two();
    float dr_three();
    float dr_four();
    float dr_five();
    float dr_truth();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Order by pT /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  RVec<float> pt_order(const RVec<float>& fv_comp, const RVec<float>& pt);

  RVec<char> pt_order_nonfloat(const RVec<char>& comp, const RVec<float>& pt);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int jet_size(const RVec<float>& j_pt);
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
  // ============================================================
  // dR_matched: returns truth-clean masses in the (l+, l-) basis
  // out[0] = m(l+, b)   , out[1] = m(l-, bbar)
  // ============================================================
  struct Lepton {
    V4 p4;
    float charge;
  
    Lepton(const V4& _p4, float _charge)
      : p4(_p4), charge(_charge) {}
  };
  

  RVec<int> dR_matched(
    const RVec<float>& b_pt,
    const RVec<float>& b_eta,
    const RVec<float>& b_phi,
    const RVec<float>& b_e,
    const RVec<float>& bbar_pt,
    const RVec<float>& bbar_eta,
    const RVec<float>& bbar_phi,
    const RVec<float>& bbar_e,
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
    const float& dR_cut
  );

  // Safe scalar extractors (now for l+/l- masses)
  float truth_m_lpb (const RVec<int>& v);
  float truth_m_lmbb(const RVec<int>& v);

  // ============================================================
  // dR truth pairing indices in the (l+, l-) basis
  // returns [jet_for_lplus, jet_for_lminus] among (0,1), else [-1,-1]
  // (Since b ↔ l+ and bbar ↔ l-, this is [jet_for_b, jet_for_bbar])
  // ============================================================
  RVec<int> dR_truth_pairing_idx_lp_lm(
    const RVec<float>& b_pt,
    const RVec<float>& b_eta,
    const RVec<float>& b_phi,
    const RVec<float>& b_e,
    const RVec<float>& bbar_pt,
    const RVec<float>& bbar_eta,
    const RVec<float>& bbar_phi,
    const RVec<float>& bbar_e,
    const RVec<float>& jet_pt,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi,
    const RVec<float>& jet_e,
    const float& dR_cut
  );

  // ============================================================
  // Enum comparison for each dR condition: truth (l+,l-) vs chi2 (l+,l-)
  // returns 0=invalid, 1=wrong, 2=correct
  // ============================================================
  // dR ONE
  RVec<int> chi2_pairing_min_mlb_by_charge1(
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
  // dR TWO
  RVec<int> chi2_pairing_min_mlb_by_charge2(
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
  // dR THREE
  RVec<int> chi2_pairing_min_mlb_by_charge3(
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
  // dR FOUR
  RVec<int> chi2_pairing_min_mlb_by_charge4(
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
  // dR FIVE
  RVec<int> chi2_pairing_min_mlb_by_charge5(
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
  // chi2 pairing in the (l+, l-) basis
  // returns [jet_for_lplus, jet_for_lminus] among (0,1), else [-1,-1]
  // NOTE: el_charge and mu_charge are floats in your ntuples.
  // ============================================================
  RVec<int> chi2_pairing_min_mlb_by_charge(
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

  int chi2_vs_dR_enum_lpb(
    const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  );

  int chi2_vs_dR_enum_lmbb(
    const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  );
}