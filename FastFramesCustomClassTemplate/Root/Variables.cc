#include "MyCustomFrame/Variables.h"
#include "FastFrames/DefineHelpers.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include "ROOT/RVec.hxx"
#include <Math/VectorUtil.h>

#include "Math/Vector4D.h"
#include "TLorentzVector.h"

#include "TH1F.h"
#include "TF1.h"
#include <tuple>

using ROOT::Math::PtEtaPhiMVector;
using ROOT::VecOps::RVec;

namespace ttZ{ //GPT aid

  // Usefull formulae for the analysis
  // Invariant mass calculator

  // Settings changer - for HyPER comparison ////////////////////////
  // Number of jets /////////////////////////////////////////////////
  // -----------------------------------------------------------------------------
  // C2.5 — Exactly two reconstructed jets -changed from 1 to 13
  // -----------------------------------------------------------------------------
  bool cutC25_exactly2jets(const RVec<float>& j_pt)
  {
      return (j_pt.size() <= 13); // CHANGED ACCORDING TO THE NEW REQUIREMENT
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.1 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Compute ΔR between two objects
  float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return ROOT::Math::VectorUtil::DeltaR(
        ROOT::Math::PtEtaPhiEVector(1, eta1, phi1, 1),
        ROOT::Math::PtEtaPhiEVector(1, eta2, phi2, 1)
    );
  }

  // Reject jets if ΔR(jet, electron) < 0.2
  ROOT::VecOps::RVec<int> jets_clean_from_e(
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi,
    const RVec<float>& el_eta,
    const RVec<float>& el_phi)
  {
    size_t nJets = jet_eta.size();
    ROOT::VecOps::RVec<int> keep(nJets, 1);

    if (el_eta.empty()) return keep;

    float e_eta = el_eta[0];
    float e_phi = el_phi[0];

    for (size_t j = 0; j < nJets; ++j) {
        float dR = deltaR(jet_eta[j], jet_phi[j], e_eta, e_phi);
        if (dR < 0.2) keep[j] = 0;  // reject jet
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
    if (el_eta.empty()) return 0;
    float e_eta = el_eta[0];
    float e_phi = el_phi[0];

    for (size_t j = 0; j < jet_eta.size(); ++j) {
        float dR = deltaR(e_eta, e_phi, jet_eta[j], jet_phi[j]);
        if (dR < 0.4) return 0;   // remove electron
    }
    return 1;  // electron survives cleaning
  }

  // Reject muons if ΔR(μ, jet) < 0.4
  int muon_clean_from_jets(
    const RVec<float>& mu_eta,
    const RVec<float>& mu_phi,
    const RVec<float>& jet_eta,
    const RVec<float>& jet_phi)
  {
    if (mu_eta.empty()) return 0;
    float m_eta = mu_eta[0];
    float m_phi = mu_phi[0];

    for (size_t j = 0; j < jet_eta.size(); ++j) {
        float dR = deltaR(m_eta, m_phi, jet_eta[j], jet_phi[j]);
        if (dR < 0.4) return 0;  // reject muon
    }
    return 1;  // muon survives
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
    return (!el_pt.empty() && el_pt[0] >= 28000);
  }
  // -----------------------------------------------------------------------------
  // A2 — Electron pseudorapidity acceptance
  // Requirement: |η(e)| < 2.47
  // -----------------------------------------------------------------------------
  bool cutA2_el_eta(const RVec<float>& el_eta) {
    return (!el_eta.empty() && std::abs(el_eta[0]) < 2.47);
  }
  // -----------------------------------------------------------------------------
  // A3 — Electron calorimeter crack veto
  // Requirement: Reject 1.37 ≤ |η| < 1.52
  // -----------------------------------------------------------------------------
  bool cutA3_el_crack(const RVec<float>& el_eta) {
    float eta = std::abs(el_eta[0]);
    return !(eta >= 1.37 && eta < 1.52);
  }
  // -----------------------------------------------------------------------------
  // A4 — Electron identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA4_el_tight(const RVec<char>& el_tight) {
    return (!el_tight.empty() && el_tight[0] == 1);
  }
  // -----------------------------------------------------------------------------
  // A5 — Muon transverse momentum requirement
  // Requirement: pT(μ) ≥ 28 GeV
  // -----------------------------------------------------------------------------
  bool cutA5_mu_pt(const RVec<float>& mu_pt) {
    return (!mu_pt.empty() && mu_pt[0] >= 28000);
  }
  // -----------------------------------------------------------------------------
  // A6 — Muon pseudorapidity acceptance
  // Requirement: |η(μ)| < 2.5
  // -----------------------------------------------------------------------------
  bool cutA6_mu_eta(const RVec<float>& mu_eta) {
    return (!mu_eta.empty() && std::abs(mu_eta[0]) < 2.5);
  }
  // -----------------------------------------------------------------------------
  // A7 — Muon identification quality
  // Requirement: tight ID flag == 1
  // -----------------------------------------------------------------------------
  bool cutA7_mu_tight(const RVec<char>& mu_tight) {
    return (!mu_tight.empty() && mu_tight[0] == 1);
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
    return (el_charge.size() == 1 &&
    mu_charge.size() == 1 &&
    el_charge[0] * mu_charge[0] == -1);
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

  // Compute m(e,μ) in GeV
  float inv_mass_elmu(const RVec<float>& el_pt, const RVec<float>& mu_pt,
    const RVec<float>& el_eta, const RVec<float>& mu_eta,
    const RVec<float>& el_phi, const RVec<float>& mu_phi)
  {
    float pt1 = el_pt[0] / 1000.f;
    float pt2 = mu_pt[0] / 1000.f;
    float dη = el_eta[0] - mu_eta[0];
    float dφ = el_phi[0] - mu_phi[0];

    return std::sqrt(2 * pt1 * pt2 * (std::cosh(dη) - std::cos(dφ)));
  }
  // -----------------------------------------------------------------------------
  // C2.4 — Invariant mass threshold m(eμ) ≥ 50 GeV
  // Requirement: m(eμ) ≥ 50 GeV
  // -----------------------------------------------------------------------------
  bool cutC24_meemu_gt50(const RVec<float>& el_pt,  const RVec<float>& mu_pt,
    const RVec<float>& el_eta, const RVec<float>& mu_eta,
    const RVec<float>& el_phi, const RVec<float>& mu_phi)
  {
    return (inv_mass_elmu(el_pt, mu_pt, el_eta, mu_eta, el_phi, mu_phi) >= 50.0f);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Section 6.3 /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    const RVec<float>& j_phi,  const RVec<float>& j_e)
  {
    PairingResult R;
    R.valid = false;

    // Need exactly 2 jets and 1 of each lepton j_pt.size() < 2 ||
    if (el_pt.empty() || mu_pt.empty())
        return R;

    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float g = 1.f / 1000.f; // MeV → GeV

    // Build 4-vectors
    V4 e (el_pt[0]*g, el_eta[0], el_phi[0], el_e[0]*g);
    V4 mu(mu_pt[0]*g, mu_eta[0], mu_phi[0], mu_e[0]*g);
    V4 j1(j_pt[0]*g, j_eta[0], j_phi[0], j_e[0]*g); // Leading jet
    V4 j2(j_pt[1]*g, j_eta[1], j_phi[1], j_e[1]*g); // Subleading jet

    // Masses for Pairing A
    float m_ej1 = (e + j1).M();
    float m_mj2 = (mu + j2).M();
    float sumA  = m_ej1*m_ej1 + m_mj2*m_mj2;

    // Masses for Pairing B
    float m_ej2 = (e + j2).M();
    float m_mj1 = (mu + j1).M();
    float sumB  = m_ej2*m_ej2 + m_mj1*m_mj1;

    // Choose minimal-sum pairing
    if (sumA < sumB) {
        R.mj1 = m_ej1;  // jet1–leptonA
        R.mj2 = m_mj2;  // jet2–leptonB
    }
    else {
        R.mj1 = m_mj1;  // jet1–leptonB
        R.mj2 = m_ej2;  // jet2–leptonA
    }

    R.valid = true;
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
    const RVec<char>&  mu_tight
  )
  {
      return (
        cutA1_el_pt(el_pt)     &&
        cutA2_el_eta(el_eta)   &&
        cutA3_el_crack(el_eta) &&
        cutA4_el_tight(el_tight) &&
        cutA5_mu_pt(mu_pt)     &&
        cutA6_mu_eta(mu_eta)   &&
        cutA7_mu_tight(mu_tight)
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
    const RVec<char>&  mu_tight
  )
  {
    return (
      section_6_1(el_pt, el_eta, el_tight, mu_pt, mu_eta, mu_tight) &&
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
  ROOT::VecOps::RVec<int> dR_matched(
    const ROOT::VecOps::RVec<float>& b_pt,
    const ROOT::VecOps::RVec<float>& b_eta,
    const ROOT::VecOps::RVec<float>& b_phi,
    const ROOT::VecOps::RVec<float>& b_e,
    const ROOT::VecOps::RVec<float>& bbar_pt,
    const ROOT::VecOps::RVec<float>& bbar_eta,
    const ROOT::VecOps::RVec<float>& bbar_phi,
    const ROOT::VecOps::RVec<float>& bbar_e,
    const ROOT::VecOps::RVec<float>& jet_pt,
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi,
    const ROOT::VecOps::RVec<float>& jet_e,
    const ROOT::VecOps::RVec<float>& el_pt,
    const ROOT::VecOps::RVec<float>& el_eta,
    const ROOT::VecOps::RVec<float>& el_phi,
    const ROOT::VecOps::RVec<float>& el_e,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_pt,
    const ROOT::VecOps::RVec<float>& mu_eta,
    const ROOT::VecOps::RVec<float>& mu_phi,
    const ROOT::VecOps::RVec<float>& mu_e,
    const ROOT::VecOps::RVec<float>& mu_charge
  ){
    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float GeV = 1.f / 1000.f;

    ROOT::VecOps::RVec<float> out;

    if (el_pt.empty() || mu_pt.empty() || el_charge.empty() || mu_charge.empty()) return out;
    if (b_pt.empty() || bbar_pt.empty() || jet_pt.size() < 2) return out;

    const V4 e  (el_pt[0]*GeV, el_eta[0], el_phi[0], el_e[0]*GeV);
    const V4 mu (mu_pt[0]*GeV, mu_eta[0], mu_phi[0], mu_e[0]*GeV);

    const V4 lplus  = (el_charge[0] > 0.f) ? e  : mu;
    const V4 lminus = (el_charge[0] > 0.f) ? mu : e;

    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);

    const V4 b    (b_pt[0]*GeV,    b_eta[0],    b_phi[0],    b_e[0]*GeV);
    const V4 bbar (bbar_pt[0]*GeV, bbar_eta[0], bbar_phi[0], bbar_e[0]*GeV);

    auto dR = [](const V4& a,const V4& b){ return ROOT::Math::VectorUtil::DeltaR(a,b); };

    int ib=-1, ibbar=-1;
    float minDRb=1e9f, minDRbb=1e9f;

    for (size_t i=0;i<jets.size();++i){
      float drb  = dR(jets[i], b);
      float drbb = dR(jets[i], bbar);

      if (drb < 0.4f && drb < minDRb)   { minDRb  = drb;  ib    = i; }
      if (drbb< 0.4f && drbb< minDRbb)  { minDRbb = drbb; ibbar = i; }
    }

    if (ib < 0 || ibbar < 0 || ib == ibbar) return out;

    out.push_back( (lplus  + jets[ib]).M() );
    out.push_back( (lminus + jets[ibbar]).M() );
    return out;
  }

  // Safe scalar extractors
  float truth_m_lpb (const ROOT::VecOps::RVec<int>& v){ return (v.size()>0 ? v[0] : -1.f); }
  float truth_m_lmbb(const ROOT::VecOps::RVec<int>& v){ return (v.size()>1 ? v[1] : -1.f); }

  // ============================================================
  // dR truth pairing indices in the (l+, l-) basis
  // returns [jet_for_lplus, jet_for_lminus] else [-1,-1]
  // ============================================================
  ROOT::VecOps::RVec<int> dR_truth_pairing_idx_lp_lm(
    const ROOT::VecOps::RVec<float>& b_pt,
    const ROOT::VecOps::RVec<float>& b_eta,
    const ROOT::VecOps::RVec<float>& b_phi,
    const ROOT::VecOps::RVec<float>& b_e,
    const ROOT::VecOps::RVec<float>& bbar_pt,
    const ROOT::VecOps::RVec<float>& bbar_eta,
    const ROOT::VecOps::RVec<float>& bbar_phi,
    const ROOT::VecOps::RVec<float>& bbar_e,
    const ROOT::VecOps::RVec<float>& jet_pt,
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi,
    const ROOT::VecOps::RVec<float>& jet_e
  ){
    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float GeV = 1.f/1000.f;
    ROOT::VecOps::RVec<int> idx(2,-1);

    if (b_pt.empty() || bbar_pt.empty() || jet_pt.size()<2) return idx;

    const V4 b    (b_pt[0]*GeV,    b_eta[0],    b_phi[0],    b_e[0]*GeV);
    const V4 bbar (bbar_pt[0]*GeV, bbar_eta[0], bbar_phi[0], bbar_e[0]*GeV);

    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);

    auto dR = [](const V4& a,const V4& b){ return ROOT::Math::VectorUtil::DeltaR(a,b); };

    float minDRb=1e9f, minDRbb=1e9f;
    int ib=-1, ibbar=-1;

    for (size_t i=0;i<jets.size();++i){
      float drb  = dR(jets[i], b);
      float drbb = dR(jets[i], bbar);

      if (drb < 0.4f && drb < minDRb)   { minDRb  = drb;  ib    = i; }
      if (drbb< 0.4f && drbb< minDRbb)  { minDRbb = drbb; ibbar = i; }
    }

    if (ib>=0 && ibbar>=0 && ib!=ibbar){
      idx[0]=ib;
      idx[1]=ibbar;
    }
    return idx;
  }

  // ============================================================
  // chi2 pairing in the (l+, l-) basis
  // ============================================================
  ROOT::VecOps::RVec<int> chi2_pairing_min_mlb_by_charge(
    const ROOT::VecOps::RVec<float>& jet_pt,
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi,
    const ROOT::VecOps::RVec<float>& jet_e,
    const ROOT::VecOps::RVec<float>& el_pt,
    const ROOT::VecOps::RVec<float>& el_eta,
    const ROOT::VecOps::RVec<float>& el_phi,
    const ROOT::VecOps::RVec<float>& el_e,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_pt,
    const ROOT::VecOps::RVec<float>& mu_eta,
    const ROOT::VecOps::RVec<float>& mu_phi,
    const ROOT::VecOps::RVec<float>& mu_e,
    const ROOT::VecOps::RVec<float>& mu_charge
  ){
    using V4 = ROOT::Math::PtEtaPhiEVector;
    constexpr float GeV = 1.f/1000.f;
    ROOT::VecOps::RVec<int> idx(2,-1);

    if (el_pt.empty() || mu_pt.empty() || el_charge.empty() || mu_charge.empty()) return idx;
    if (el_charge[0]*mu_charge[0]>=0.f || jet_pt.size()<2) return idx;

    const V4 lplus  = (el_charge[0]>0.f)
      ? V4(el_pt[0]*GeV, el_eta[0], el_phi[0], el_e[0]*GeV)
      : V4(mu_pt[0]*GeV, mu_eta[0], mu_phi[0], mu_e[0]*GeV);

    const V4 lminus = (el_charge[0]>0.f)
      ? V4(mu_pt[0]*GeV, mu_eta[0], mu_phi[0], mu_e[0]*GeV)
      : V4(el_pt[0]*GeV, el_eta[0], el_phi[0], el_e[0]*GeV);

    std::vector<V4> jets;
    for (size_t i=0;i<jet_pt.size();++i)
      jets.emplace_back(jet_pt[i]*GeV, jet_eta[i], jet_phi[i], jet_e[i]*GeV);

    constexpr float mu_lpb=97.82f, mu_lmbb=97.93f;
    constexpr float s2_lp=29.95f*29.95f, s2_lm=29.98f*29.98f;

    auto chi2=[&](float a,float b){
      return (a-mu_lpb)*(a-mu_lpb)/s2_lp + (b-mu_lmbb)*(b-mu_lmbb)/s2_lm;
    };

    float min_chi2 = 1e12f;
    int best_i=-1, best_j=-1;

    // Loop over all unique jet pairs
    for (size_t i=0; i<jets.size(); ++i){
        for (size_t j=i+1; j<jets.size(); ++j){
            float chiA = chi2( (lplus + jets[i]).M(), (lminus + jets[j]).M() );
            float chiB = chi2( (lplus + jets[j]).M(), (lminus + jets[i]).M() );

            if (chiA < min_chi2) { min_chi2 = chiA; best_i = i; best_j = j; }
            if (chiB < min_chi2) { min_chi2 = chiB; best_i = j; best_j = i; }
        }
    }

    if (best_i >= 0 && best_j >= 0) { idx[0] = best_i; idx[1] = best_j; }
    return idx;
  }

  // ============================================================
  // Enum comparison: truth (l+,l-) vs chi2 (l+,l-)
  // returns 0=invalid, 1=wrong, 2=correct
  // ============================================================
  int chi2_vs_dR_enum(
    const ROOT::VecOps::RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const ROOT::VecOps::RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  ){
    if (truth_lp_lm.size() != 2 || chi2_lp_lm.size() != 2) return 0;
    if (truth_lp_lm[0] < 0 || truth_lp_lm[1] < 0) return 0;
    if (chi2_lp_lm[0]  < 0 || chi2_lp_lm[1]  < 0) return 0;

    const bool correct =
      (chi2_lp_lm[0] == truth_lp_lm[0]) &&
      (chi2_lp_lm[1] == truth_lp_lm[1]);

    return correct ? 2 : 1;
  }

  // Per branch
  // ============================================================
  // Per-branch enum: l+ ↔ b
  // returns 0=invalid, 1=wrong, 2=correct
  // ============================================================
  int chi2_vs_dR_enum_lpb(
    const ROOT::VecOps::RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const ROOT::VecOps::RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
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
    const ROOT::VecOps::RVec<int>& truth_lp_lm, // [jet_for_lplus, jet_for_lminus]
    const ROOT::VecOps::RVec<int>& chi2_lp_lm   // [jet_for_lplus, jet_for_lminus]
  ){
    if (truth_lp_lm.size() != 2 || chi2_lp_lm.size() != 2) return 0;
    if (truth_lp_lm[1] < 0 || chi2_lp_lm[1] < 0) return 0;

    return (chi2_lp_lm[1] == truth_lp_lm[1]) ? 2 : 1;
  }
} // namespace ttZ