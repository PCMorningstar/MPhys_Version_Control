#include "MyCustomFrame/Variables.h"
#include "FastFrames/DefineHelpers.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include "ROOT/RVec.hxx"
#include <Math/VectorUtil.h>
#include <numeric>   // std::iota
#include <algorithm> // std::sort (if used nearby)
#include <type_traits> // std::is_same

#include "Math/Vector4D.h"
#include "TLorentzVector.h"

#include "TH1F.h"
#include "TF1.h"
#include <tuple>

using V4 = ROOT::Math::PtEtaPhiEVector;
using ROOT::VecOps::RVec;

namespace ttZ{ //GPT aid

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Order by pT /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // object ordering, by pT for each 4 vec component
  RVec<float> pt_order(const RVec<float>& fv_comp, const RVec<float>& pt)
  {
    auto obj_pt_sorted = ROOT::VecOps::Take(fv_comp,
      ROOT::VecOps::Argsort(pt, [](float a, float b){ return a > b;}));

    return obj_pt_sorted;
  }

  RVec<char> pt_order_nonfloat(const RVec<char>& comp, const RVec<float>& pt)
  {
    auto obj_pt_sorted = ROOT::VecOps::Take(comp,
      ROOT::VecOps::Argsort(pt, [](float a, float b){ return a > b;}));

    return obj_pt_sorted;
  }

  // Usefull formulae for the analysis
  // Invariant mass calculator

  // Settings changer - for HyPER comparison ////////////////////////
  // Number of jets /////////////////////////////////////////////////
  // -----------------------------------------------------------------------------
  // C2.5 — Exactly two reconstructed jets -changed from 1 to 13
  // -----------------------------------------------------------------------------
  bool cutC25_exactly2jets(const RVec<float>& j_pt)
  {
      return (j_pt.size() >= 0); // CHANGED ACCORDING TO THE NEW REQUIREMENT
  }

  int jet_size(const RVec<float>& j_pt)
  {
    return (j_pt.size());
  }

  // creating multiply dRs for efficiency v. dR
  float dr_one()
  {
    return (0.5f);
  }
  float dr_two()
  {
    return (2.0f);
  }
  float dr_three()
  {
    return (3.0f);
  }
  float dr_four()
  {
    return (4.0f);
  }
  float dr_five()
  {
    return (5.0f);
  }


  float dr_truth()
  {
    return (0.4f);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Compute ΔR between two objects
  float deltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2) {
    return ROOT::Math::VectorUtil::DeltaR(
        ROOT::Math::PtEtaPhiEVector(1, eta1, phi1, 1),
        ROOT::Math::PtEtaPhiEVector(1, eta2, phi2, 1)
    );
  }

  // Reject jets if ΔR(jet, electron) < 0.2
  RVec<int> jets_clean_from_e(
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi)
  {
    size_t nJets = jet_eta.size();
    RVec<int> keep(nJets, 1);

    if (el_eta.empty()) return keep;

    for (size_t m = 0; m < el_eta.size(); ++m) {
      float eta_val = el_eta[m];
      float phi_val = el_phi[m];

      for (size_t j = 0; j < nJets; ++j) {
          float dR = deltaR(jet_eta[j], jet_phi[j], eta_val, phi_val);
          if (dR < 0.2) keep[j] = 0;  // reject jet
      }
    }

    return keep;
  }

  // Reject electrons if ΔR(e, jet) < 0.4
  int electron_clean_from_jets(
    const RVec<float>& el_eta,
    const RVec<float>& el_phi,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi)
  {
    // if none, consider "fail" or "pass" depending on your convention
    if (el_eta.empty()) return 0;

    for (size_t m = 0; m < el_eta.size(); ++m) {
        float eta_val = el_eta[m];
        float phi_val = el_phi[m];

        for (size_t j = 0; j < jet_eta.size(); ++j) {
            float dR_val = deltaR(eta_val, phi_val, jet_eta[j], jet_phi[j]);
            if (dR_val < 0.4) {
                return 0;  // reject: this muon is too close to a jet
            }
        }
    }
    return 1;  // all muons are clean
  }

  // Reject muons if ΔR(μ, jet) < 0.4
  int muon_clean_from_jets(
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi)
  {
    // if no muons, consider "fail" or "pass" depending on your convention
    if (mu_eta.empty()) return 0;

    for (size_t m = 0; m < mu_eta.size(); ++m) {
        float m_eta = mu_eta[m];
        float m_phi = mu_phi[m];

        for (size_t j = 0; j < jet_eta.size(); ++j) {
            float dR_val = deltaR(m_eta, m_phi, jet_eta[j], jet_phi[j]);
            if (dR_val < 0.4) {
                return 0;  // reject: this muon is too close to a jet
            }
        }
    }
    return 1;  // all muons are clean
  }

  bool cutA8_jet_clean(const RVec<int>& jet_keep_flags) {
    // True if at least 1 jet remains after cleaning
    return ROOT::VecOps::Sum(jet_keep_flags) >= 1;
  }

  bool cutA9_el_clean(const int& el_keep_flag) {
      return el_keep_flag == 1;
  }

  bool cutA10_mu_clean(const int& mu_keep_flag) {
      return mu_keep_flag == 1;
  }


  // -----------------------------------------------------------------------------
  // A1 — Electron transverse momentum requirement
  // Requirement: pT(e) ≥ 28 GeV
  // -----------------------------------------------------------------------------
  bool cutA1_el_pt(const RVec<float>& el_pt) {
    if (el_pt.empty()) return false;  // fail if no electrons
    for (float pt : el_pt) {
        if (pt < 28000) return false;  // fail if any electron below threshold
    }
    return true;  // all electrons above threshold
  }
  // -----------------------------------------------------------------------------
  // A2 — Electron pseudorapidity acceptance
  // Requirement: |η(e)| < 2.47
  // -----------------------------------------------------------------------------
  bool cutA2_el_eta(const RVec<float>& el_eta) {

    if (el_eta.empty()) return false;  // fail if none
    for (float eta : el_eta) {
      if (std::abs(eta) >= 2.47) return false;  // fail if any electron below threshold
    }
    return true;  // all electrons above threshold
  }
  // -----------------------------------------------------------------------------
  // A3 — Electron calorimeter crack veto
  // Requirement: Reject 1.37 ≤ |η| < 1.52
  // -----------------------------------------------------------------------------
  bool cutA3_el_crack(const RVec<float>& el_eta) {
    if (el_eta.empty()) return false;  // fail if none
    for (float eta : el_eta) {
      if(std::abs(eta) >= 1.37 && std::abs(eta) < 1.52) return false;
    }
    return true;
  }
  // -----------------------------------------------------------------------------
  // A4 — Electron identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA4_el_tight(const RVec<char>& el_tight) {
    if (el_tight.empty()) return false;  // fail if none
    for (float id : el_tight) {
      if(id == 0) return false;
    }
    return true;
  }
  // -----------------------------------------------------------------------------
  // A5 — Muon transverse momentum requirement
  // Requirement: pT(μ) ≥ 28 GeV
  // -----------------------------------------------------------------------------
  bool cutA5_mu_pt(const RVec<float>& mu_pt) {
    if (mu_pt.empty()) return false;  // fail if none
    for (float pt : mu_pt) {
        if (pt < 28000) return false;  // fail if any below threshold
    }
    return true;  // all above threshold
  }
  // -----------------------------------------------------------------------------
  // A6 — Muon pseudorapidity acceptance
  // Requirement: |η(μ)| < 2.5
  // -----------------------------------------------------------------------------
  bool cutA6_mu_eta(const RVec<float>& mu_eta) {
    if (mu_eta.empty()) return false;  // fail if none
    for (float eta : mu_eta) {
      if (std::abs(eta) >= 2.5) return false;
    }
    return true;
  }
  // -----------------------------------------------------------------------------
  // A7 — Muon identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA7_mu_tight(const RVec<char>& mu_tight) {
    if (mu_tight.empty()) return false;  // fail if none
    for (float id : mu_tight) {
      if(id == 0) return false;
    }
    return true;
  }
  // -----------------------------------------------------------------------------
  // A7.1 — Jet identification quality
  // Requirement: select jvt flag == 1 (like the lepton tight_ID)
  // -----------------------------------------------------------------------------
  bool cutA7_1_jet_jvt(const RVec<char>& jet_jvt)
  {
    if (jet_jvt.empty()) return false;  // fail if none
    for (float id : jet_jvt) {
      if(id == 0) return false;
    }
    return true;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.2 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------
  // C2.1 — Exactly two prompt leptons
  // Requirement: N(leptons) == 2
  // -----------------------------------------------------------------------------
  bool cutC21_exactly2leptons(const RVec<float>& el_pt,
    const RVec<float>& mu_pt)
  {
    int nlep = el_pt.size() + mu_pt.size();
    return (nlep == 2);
  }
  // -----------------------------------------------------------------------------
  // C2.2 — Opposite-charge requirement for dileptons
  // Requirement: q(e) * q(μ) == -1
  // -----------------------------------------------------------------------------
  bool cutC22_opposite_charge(const RVec<float>& el_charge,
    const RVec<float>& mu_charge)
  {
    if (el_charge.empty() || mu_charge.empty()) return false;

    // Loop over all electron-muon pairs
    for (size_t i = 0; i < el_charge.size(); ++i) {
      for (size_t j = 0; j < mu_charge.size(); ++j) {
          if (el_charge[i] * mu_charge[j] == -1) return true;  // any opposite-charge pair
      }
    }
    return false;  // no opposite-charge pair found
  }
  // -----------------------------------------------------------------------------
  // C2.3 — eμ final-state requirement
  // Requirement: Exactly 1 electron and exactly 1 muon
  // -----------------------------------------------------------------------------
  bool cutC23_one_el_one_mu(const RVec<float>& el_pt,
    const RVec<float>& mu_pt)
  {
    return (el_pt.size() == 1 && mu_pt.size() == 1);
  }
  //-----------------------------------------------------------------------------
  // C2.4 — Invariant mass threshold m(eμ) ≥ 50 GeV
  // Requirement: m(eμ) ≥ 50 GeV
  // -----------------------------------------------------------------------------
  bool cutC24_meemu_gt50(const RVec<float>& el_pt,  const RVec<float>& mu_pt,
    const RVec<float>& el_eta, const RVec<float>& mu_eta,
    const RVec<float>& el_phi, const RVec<float>& mu_phi)
  {
    if (el_pt.empty() || mu_pt.empty()) return false;

    // Loop over all electron-muon pairs
    for (size_t i = 0; i < el_pt.size(); ++i) {
      for (size_t j = 0; j < mu_pt.size(); ++j) {
        float pt1 = el_pt[i] / 1000.f;  // MeV -> GeV
        float pt2 = mu_pt[j] / 1000.f;
        float dEta = el_eta[i] - mu_eta[j];
        float dPhi = el_phi[i] - mu_phi[j];
        float mass = std::sqrt(2 * pt1 * pt2 * (std::cosh(dEta) - std::cos(dPhi)));

        if (mass >= 50.0f) return true;  // passes if any pair satisfies
      }
    }
    return false;  // no pair passed
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.3 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------
  // Compute Section 6.3 minimal pairing
  //
  // Purpose:
  //   Given all electron, all muon, and 2 jets (j1 = leading, j2 = subleading),
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
    const RVec<float>& j_phi,  const RVec<float>& j_e)
  {
    PairingResult R;
    constexpr float g = 1.f / 1000.f; // MeV → GeV

    // -------------------------------
    // Require at least 1 electron, 1 muon, 2 jets
    // -------------------------------
    if (el_pt.empty() || mu_pt.empty() || j_pt.size() < 2) return R;

    using V4 = ROOT::Math::PtEtaPhiEVector;

    // -------------------------------
    // Identify the two leading jets
    // -------------------------------
    size_t i_j1 = 0, i_j2 = 1;
    if (j_pt[1] > j_pt[0]) std::swap(i_j1, i_j2);

    V4 j1(j_pt[i_j1]*g, j_eta[i_j1], j_phi[i_j1], j_e[i_j1]*g);
    V4 j2(j_pt[i_j2]*g, j_eta[i_j2], j_phi[i_j2], j_e[i_j2]*g);

    float min_sum = std::numeric_limits<float>::max();

    // -------------------------------
    // Loop over all electron-muon pairs
    // -------------------------------
    for (size_t i_e = 0; i_e < el_pt.size(); ++i_e) {
      V4 e(el_pt[i_e]*g, el_eta[i_e], el_phi[i_e], el_e[i_e]*g);

    for (size_t i_mu = 0; i_mu < mu_pt.size(); ++i_mu) {
      V4 mu(mu_pt[i_mu]*g, mu_eta[i_mu], mu_phi[i_mu], mu_e[i_mu]*g);

      // Pairing A: e→j1, μ→j2
      float m_ej1 = (e + j1).M();
      float m_mj2 = (mu + j2).M();
      float sumA  = m_ej1*m_ej1 + m_mj2*m_mj2;

      // Pairing B: e→j2, μ→j1
      float m_ej2 = (e + j2).M();
      float m_mj1 = (mu + j1).M();
      float sumB  = m_ej2*m_ej2 + m_mj1*m_mj1;

        // -------------------------------
        // Select the pairing with minimal sum
        // -------------------------------
        if (sumA < min_sum) {
          min_sum = sumA;
          R.mj1 = m_ej1;
          R.mj2 = m_mj2;
          R.valid = true;
        }
        if (sumB < min_sum) {
          min_sum = sumB;
          R.mj1 = m_mj1;
          R.mj2 = m_ej2;
          R.valid = true;
        }
      } // mu loop
    } // e loop

    // -------------------------------
    // Apply Section 6.3 thresholds
    // -------------------------------
    if (R.valid) {
        if (R.mj1 < 20.f || R.mj2 < 20.f) R.valid = false;
    }

    return R;
  }
  // -----------------------------------------------------------------------------
  // D3.1 — Valid pairing computed
  // Requirement: compute_pairing_S6_3() returned valid = true
  // -----------------------------------------------------------------------------
  bool cutD31_pairing_valid(const PairingResult& P) {
    return P.valid;
  }
  // -----------------------------------------------------------------------------
  // D3.2 — Leading pairing mass threshold
  // Requirement: mj1 ≥ 20 GeV
  // -----------------------------------------------------------------------------
  bool cutD32_mj1(const PairingResult& P) {
    return (P.mj1 >= 20.f);
  }
  // -----------------------------------------------------------------------------
  // D3.3 — Subleading pairing mass threshold
  // Requirement: mj2 ≥ 20 GeV
  // -----------------------------------------------------------------------------
  bool cutD33_mj2(const PairingResult& P) {
    return (P.mj2 >= 20.f);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool section_6_1(
    const RVec<float>& el_pt,
    const RVec<float>& el_eta,
    const RVec<char>&  el_tight,
    const RVec<float>& mu_pt,
    const RVec<float>& mu_eta,
    const RVec<char>&  mu_tight,
    const RVec<char>& jet_jvt
  )
  {
      return (
        cutA1_el_pt(el_pt)     &&
        cutA2_el_eta(el_eta)   &&
        cutA3_el_crack(el_eta) &&
        cutA4_el_tight(el_tight) &&
        cutA5_mu_pt(mu_pt)     &&
        cutA6_mu_eta(mu_eta)   &&
        cutA7_mu_tight(mu_tight) &&
        cutA7_1_jet_jvt(jet_jvt)
      );
  }
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
  )
  {
    return (
      cutC21_exactly2leptons(el_pt, mu_pt) &&
      cutC22_opposite_charge(el_charge, mu_charge) &&
      cutC23_one_el_one_mu(el_pt, mu_pt) &&
      cutC24_meemu_gt50(el_pt, mu_pt, el_eta, mu_eta, el_phi, mu_phi) &&
      cutC25_exactly2jets(j_pt)
    );
  }
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
  )
  {
    // Compute minimal pairing
    PairingResult P = compute_pairing_S6_3(
        el_pt, el_eta, el_phi, el_e,
        mu_pt, mu_eta, mu_phi, mu_e,
        j_pt, j_eta, j_phi, j_e
    );
    // Apply classification cuts
    return (
        cutD31_pairing_valid(P) &&
        cutD32_mj1(P) &&
        cutD33_mj2(P)
    );
  }

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
    const RVec<char>&  el_tight,
    const RVec<char>&  mu_tight,
    const RVec<char>& jet_jvt
  )
  {
    return (
      section_6_1(el_pt, el_eta, el_tight, mu_pt, mu_eta, mu_tight, jet_jvt) &&
      section_6_2(el_pt, mu_pt, el_eta, mu_eta, el_phi, mu_phi, el_charge, mu_charge, j_pt) &&
      section_6_3(el_pt, el_eta, el_phi, el_e, mu_pt, mu_eta, mu_phi, mu_e, j_pt, j_eta, j_phi, j_e) &&

      cutA8_jet_clean(jet_keep_flags) &&
      cutA9_el_clean(el_keep_flag) &&
      cutA10_mu_clean(mu_keep_flag)
    );
  }


  // ============================================================
  // dR_matched: returns truth-clean masses in the (l+, l-) basis
  // out[0] = m(l+, b)   , out[1] = m(l-, bbar)
  // ============================================================
  RVec<int> dR_matched(
    const RVec<float>& b_pti,
    const RVec<float>& b_etai,
    const RVec<float>& b_phii,
    const RVec<float>& b_ei,
    const RVec<float>& bbar_pti,
    const RVec<float>& bbar_etai,
    const RVec<float>& bbar_phii,
    const RVec<float>& bbar_ei,
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
  ){
    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float GeV = 1.f / 1000.f;

    RVec<float> out;

    // ----------------------------------------------------------
    // Require at least two leptons and two jets
    // ----------------------------------------------------------
    if ((el_pt.size() + mu_pt.size()) < 2) return out;
    if (jet_pt.size() < 2) return out;
    if (b_pti.empty() || bbar_pti.empty()) return out;

    // ----------------------------------------------------------
    // Sort truth b and bbar by pT
    // ----------------------------------------------------------
    auto idx_b    = ROOT::VecOps::Argsort(b_pti, [](float a,float b){return a>b;});
    auto idx_bbar = ROOT::VecOps::Argsort(bbar_pti, [](float a,float b){return a>b;});

    V4 b    (b_pti[idx_b[0]]    *GeV, b_etai[idx_b[0]],    b_phii[idx_b[0]],    b_ei[idx_b[0]]    *GeV);
    V4 bbar (bbar_pti[idx_bbar[0]]*GeV, bbar_etai[idx_bbar[0]], bbar_phii[idx_bbar[0]], bbar_ei[idx_bbar[0]]*GeV);

    // ----------------------------------------------------------
    // Build lepton collection
    // ----------------------------------------------------------
    RVec<Lepton> leptons;

    for(size_t i=0;i<el_pt.size();++i)
      leptons.emplace_back(
        V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV),
        el_charge[i]
      );

    for(size_t i=0;i<mu_pt.size();++i)
      leptons.emplace_back(
        V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV),
        mu_charge[i]
      );

    if(leptons.size() != 2) return out;
    if(leptons[0].charge * leptons[1].charge >= 0) return out;

    const V4& lplus  = (leptons[0].charge > 0.f) ? leptons[0].p4 : leptons[1].p4;
    const V4& lminus = (leptons[0].charge < 0.f) ? leptons[0].p4 : leptons[1].p4;

    // ----------------------------------------------------------
    // Build reco jets
    // ----------------------------------------------------------
    std::vector<V4> jets;
    for(size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);

    auto dR = [](const V4& a,const V4& b){
      return ROOT::Math::VectorUtil::DeltaR(a,b);
    };

    // ----------------------------------------------------------
    // GLOBAL 2D ΔR MINIMISATION
    // ----------------------------------------------------------
    float best_metric = 1e9;
    int best_ib = -1;
    int best_ibbar = -1;

    for(size_t i=0;i<jets.size();++i){
      for(size_t j=0;j<jets.size();++j){

        if(i==j) continue;

        float dr_b    = dR(jets[i], b);
        float dr_bbar = dR(jets[j], bbar);

        if(dr_b > dR_cut || dr_bbar > dR_cut) continue;

        float metric = dr_b*dr_b + dr_bbar*dr_bbar;

        if(metric < best_metric){
          best_metric = metric;
          best_ib = i;
          best_ibbar = j;
        }
      }
    }

    if(best_ib < 0 || best_ibbar < 0) return out;

    // ----------------------------------------------------------
    // Return invariant masses
    // ----------------------------------------------------------
    out.push_back( (lplus  + jets[best_ib]).M() );
    out.push_back( (lminus + jets[best_ibbar]).M() );

    return out;
  }

  // Safe scalar extractors
  float truth_m_lpb (const RVec<int>& v){ return (v.size()>0 ? v[0] : -1.f); }
  float truth_m_lmbb(const RVec<int>& v){ return (v.size()>1 ? v[1] : -1.f); }

  // ============================================================
  // dR truth pairing indices in the (l+, l-) basis
  // returns [jet_for_lplus, jet_for_lminus] else [-1,-1]
  // ============================================================
  RVec<int> dR_truth_pairing_idx_lp_lm(
    const RVec<float>& b_pti,
    const RVec<float>& b_etai,
    const RVec<float>& b_phii,
    const RVec<float>& b_ei,
    const RVec<float>& bbar_pti,
    const RVec<float>& bbar_etai,
    const RVec<float>& bbar_phii,
    const RVec<float>& bbar_ei,
    const RVec<float>& jet_pt,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi,
    const RVec<float>& jet_e,
    const float& dR_cut
  ){
    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float GeV = 1.f / 1000.f;

    RVec<int> out;

    // Require at least two jets and truth b/bbar
    if (jet_pt.size() < 2) return out;
    if (b_pti.empty() || bbar_pti.empty()) return out;

    // ----------------------------------------------------------
    // Sort truth b and bbar by pT (highest pT first)
    // ----------------------------------------------------------
    auto idx_b    = ROOT::VecOps::Argsort(b_pti, [](float a,float b){return a>b;});
    auto idx_bbar = ROOT::VecOps::Argsort(bbar_pti, [](float a,float b){return a>b;});

    V4 b    (b_pti[idx_b[0]]*GeV, b_etai[idx_b[0]],    b_phii[idx_b[0]],    b_ei[idx_b[0]]*GeV);
    V4 bbar (bbar_pti[idx_bbar[0]]*GeV, bbar_etai[idx_bbar[0]], bbar_phii[idx_bbar[0]], bbar_ei[idx_bbar[0]]*GeV);

    // ----------------------------------------------------------
    // Build reco jets
    // ----------------------------------------------------------
    std::vector<V4> jets;
    jets.reserve(jet_pt.size());

    for(size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);

    auto dR = [](const V4& a,const V4& b){
      return ROOT::Math::VectorUtil::DeltaR(a,b);
    };

    // ----------------------------------------------------------
    // GLOBAL 2D MINIMISATION OVER JET PAIRS
    // ----------------------------------------------------------
    float best_metric = 1e9f;
    int best_ib = -1;
    int best_ibbar = -1;

    for(size_t i=0;i<jets.size();++i){
      for(size_t j=0;j<jets.size();++j){

        if(i == j) continue;

        float dr_b    = dR(jets[i], b);
        float dr_bbar = dR(jets[j], bbar);

        if(dr_b > dR_cut || dr_bbar > dR_cut) continue;

        float metric = dr_b*dr_b + dr_bbar*dr_bbar;

        if(metric < best_metric){
          best_metric = metric;
          best_ib = i;
          best_ibbar = j;
        }
      }
    }

    if(best_ib < 0 || best_ibbar < 0) return out;

    out.push_back(best_ib);
    out.push_back(best_ibbar);

    return out;
  }

  // ============================================================
  // chi2 pairing in the (l+, l-) basis
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2, -1);  // [jet for l+, jet for l-]

    // chi2 constants
    constexpr float mu_lpb = 89.27f, mu_lmbb = 95.45f;
    constexpr float s_lp = 28.25f, s_lm = 29.33f;
    constexpr float s2_lp = s_lp*s_lp, s2_lm = s_lm*s_lm;

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i = 0; i < el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i = 0; i < mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() != 2) return idx; // Enforces dilepton (precaution)

    if (leptons[0].charge * leptons[1].charge >= 0) // Enforces opposite charge (precaution)
    return idx;
    const V4& lplus  = (leptons[0].charge > 0) ? leptons[0].p4 : leptons[1].p4;
    const V4& lminus = (leptons[0].charge < 0) ? leptons[0].p4 : leptons[1].p4;

    // Build jets
    std::vector<V4> jets;
    for (size_t i = 0; i < jet_pt.size(); ++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 
    auto chi2 = [&](float a, float b) {
        return (a - mu_lpb)*(a - mu_lpb)/s2_lp + (b - mu_lmbb)*(b - mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs for 2D minimization
    float min_chi2 = 1e12f;
    int best_i = -1, best_j = -1;

    for (size_t i = 0; i < jets.size(); ++i) {
      for (size_t j = 0; j < jets.size(); ++j) {
        if (i == j) continue;  // avoid assigning same jet to both leptons

        float chi_val = chi2((lplus + jets[i]).M(), (lminus + jets[j]).M());

        if (chi_val < min_chi2) {
          min_chi2 = chi_val;
          best_i = i;
          best_j = j;
        }
      }
    }

    if (best_i >= 0 && best_j >= 0) {
      idx[0] = best_i;   // jet assigned to l+
      idx[1] = best_j;   // jet assigned to l-
    }
    return idx;
  }
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2,-1);

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i=0; i<el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i=0; i<mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() < 2) return idx;

    // Pick one l+ and one l−
    const V4* lplus = nullptr;
    const V4* lminus = nullptr;
    for (const auto& lep : leptons){
        if (lep.charge > 0.f && !lplus)  lplus = &lep.p4;
        if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
    if (!lplus || !lminus) return idx;

    // Build jets
    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 constants
    constexpr float mu_lpb = 97.79f, mu_lmbb = 97.95f;
    constexpr float s2_lp = 29.95f*29.95f, s2_lm = 29.99f*29.99f;
    auto chi2 = [&](float a, float b){
        return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs
    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (*lplus + jets[i]).M(), (*lminus + jets[j]).M() );
            float chiB = chi2( (*lplus + jets[j]).M(), (*lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2,-1);

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i=0; i<el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i=0; i<mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() < 2) return idx;

    // Pick one l+ and one l−
    const V4* lplus = nullptr;
    const V4* lminus = nullptr;
    for (const auto& lep : leptons){
        if (lep.charge > 0.f && !lplus)  lplus = &lep.p4;
        if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
    if (!lplus || !lminus) return idx;

    // Build jets
    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 constants
    constexpr float mu_lpb = 97.79f, mu_lmbb = 97.95f;
    constexpr float s2_lp = 29.95f*29.95f, s2_lm = 29.99f*29.99f;
    auto chi2 = [&](float a, float b){
        return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs
    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (*lplus + jets[i]).M(), (*lminus + jets[j]).M() );
            float chiB = chi2( (*lplus + jets[j]).M(), (*lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2,-1);

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i=0; i<el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i=0; i<mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() < 2) return idx;

    // Pick one l+ and one l−
    const V4* lplus = nullptr;
    const V4* lminus = nullptr;
    for (const auto& lep : leptons){
        if (lep.charge > 0.f && !lplus)  lplus = &lep.p4;
        if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
    if (!lplus || !lminus) return idx;

    // Build jets
    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 constants
    constexpr float mu_lpb = 97.79f, mu_lmbb = 97.95f;
    constexpr float s2_lp = 29.95f*29.95f, s2_lm = 29.99f*29.99f;
    auto chi2 = [&](float a, float b){
        return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs
    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (*lplus + jets[i]).M(), (*lminus + jets[j]).M() );
            float chiB = chi2( (*lplus + jets[j]).M(), (*lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2,-1);

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i=0; i<el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i=0; i<mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() < 2) return idx;

    // Pick one l+ and one l−
    const V4* lplus = nullptr;
    const V4* lminus = nullptr;
    for (const auto& lep : leptons){
        if (lep.charge > 0.f && !lplus)  lplus = &lep.p4;
        if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
    if (!lplus || !lminus) return idx;

    // Build jets
    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 constants
    constexpr float mu_lpb = 97.79f, mu_lmbb = 97.95f;
    constexpr float s2_lp = 29.95f*29.95f, s2_lm = 29.99f*29.99f;
    auto chi2 = [&](float a, float b){
        return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs
    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (*lplus + jets[i]).M(), (*lminus + jets[j]).M() );
            float chiB = chi2( (*lplus + jets[j]).M(), (*lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }
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
  ){
    constexpr float GeV = 1.f/1000.f;
    RVec<int> idx(2,-1);

    // Collect all leptons
    RVec<Lepton> leptons;
    for (size_t i=0; i<el_pt.size(); ++i)
        leptons.push_back(Lepton{V4(el_pt[i]*GeV, el_eta[i], el_phi[i], el_e[i]*GeV), el_charge[i]});
    for (size_t i=0; i<mu_pt.size(); ++i)
        leptons.push_back(Lepton{V4(mu_pt[i]*GeV, mu_eta[i], mu_phi[i], mu_e[i]*GeV), mu_charge[i]});

    if (leptons.size() < 2) return idx;

    // Pick one l+ and one l−
    const V4* lplus = nullptr;
    const V4* lminus = nullptr;
    for (const auto& lep : leptons){
        if (lep.charge > 0.f && !lplus)  lplus = &lep.p4;
        if (lep.charge < 0.f && !lminus) lminus = &lep.p4;
    }
    if (!lplus || !lminus) return idx;

    // Build jets
    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
        jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);
    if (jets.size() < 2) return idx;

    // chi2 constants
    constexpr float mu_lpb = 97.79f, mu_lmbb = 97.95f;
    constexpr float s2_lp = 29.95f*29.95f, s2_lm = 29.99f*29.99f;
    auto chi2 = [&](float a, float b){
        return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    // Loop over all unique jet pairs
    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (*lplus + jets[i]).M(), (*lminus + jets[j]).M() );
            float chiB = chi2( (*lplus + jets[j]).M(), (*lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }

  // Per branch
  // ============================================================
  // Per-branch enum: l+ ↔ b
  // returns 0=invalid, 1=wrong, 2=correct
  // ============================================================
  int chi2_vs_dR_enum_lpb(
    const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  ){
    if (truth_lp_lm.size() != 2 || chi2_lp_lm.size() != 2) return 0;
    if (truth_lp_lm[0] < 0 || chi2_lp_lm[0] < 0) return 0;

    return (chi2_lp_lm[0] == truth_lp_lm[0]) ? 2 : 1;
  }

  // ============================================================
  // Per-branch enum: l- ↔ bbar
  // returns 0=invalid, 1=wrong, 2=correct
  // ============================================================
  int chi2_vs_dR_enum_lmbb(
    const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  ){
    if (truth_lp_lm.size() != 2 || chi2_lp_lm.size() != 2) return 0;
    if (truth_lp_lm[1] < 0 || chi2_lp_lm[1] < 0) return 0;

    return (chi2_lp_lm[1] == truth_lp_lm[1]) ? 2 : 1;
  }


// ============================================================
// Minimum Invariant Squared Mass Sum (MISMS)
// pairing in the (l+, l-) basis
// Minimises:  M(l+ b)^2 + M(l- b)^2
// ============================================================

RVec<int> misms_pairing_min_mlb2_by_charge(
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
  RVec<int> idx(2,-1);

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
        (*lplus  + jets[i]).M2() +
        (*lminus + jets[j]).M2();

      // Assignment B: l+ → j, l− → i
      float sumB =
        (*lplus  + jets[j]).M2() +
        (*lminus + jets[i]).M2();

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
  if (best_i >= 0 && best_j >= 0){
    idx[0] = best_i;  // l+
    idx[1] = best_j;  // l−
  }
  return idx;
}
// Per branch
// ============================================================
// Per-branch enum: l+ ↔ b
// returns 0=invalid, 1=wrong, 2=correct
// ============================================================
int misms_pairing_min_mlb2_lpb(
  const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
  const RVec<int>& misms_lp_lm   // [jet_for_lplus, jet_for_lminus]
){
  if (truth_lp_lm.size() != 2 || misms_lp_lm.size() != 2) return 0;
  if (truth_lp_lm[0] < 0 || misms_lp_lm[0] < 0) return 0;

  return (misms_lp_lm[0] == truth_lp_lm[0]) ? 2 : 1;
}
// ============================================================
// Per-branch enum: l- ↔ bbar
// returns 0=invalid, 1=wrong, 2=correct
// ============================================================
int misms_pairing_min_mlb2_lmbb(
  const RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
  const RVec<int>& misms_lp_lm   // [jet_for_lplus, jet_for_lminus]
){
  if (truth_lp_lm.size() != 2 || misms_lp_lm.size() != 2) return 0;
  if (truth_lp_lm[1] < 0 || misms_lp_lm[1] < 0) return 0;

  return (misms_lp_lm[1] == truth_lp_lm[1]) ? 2 : 1;
}

} // namespace ttZ