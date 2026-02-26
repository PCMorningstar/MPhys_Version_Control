// m_lb(+)
chi += std::pow(
    (vis_plus.M() - obs_map["mlb_plus"].mean)
    / obs_map["mlb_plus"].sigma,
    2.0
);

// m_lb(-)
chi += std::pow(
    (vis_minus.M() - obs_map["mlb_minus"].mean)
    / obs_map["mlb_minus"].sigma,
    2.0
);

// pTdiff  (pT(l+b) difference)
double ptDiff = vis_plus.Pt() - vis_minus.Pt();

chi += std::pow(
    (ptDiff - obs_map["pTdiff"].mean)
    / obs_map["pTdiff"].sigma,
    2.0
);

// Sum ΔR
chi += std::pow(
    (ROOT::Math::VectorUtil::DeltaR(lplus, jet1)
   + ROOT::Math::VectorUtil::DeltaR(lminus, jet2)
   - obs_map["sum_deltaR"].mean)
   / obs_map["sum_deltaR"].sigma,
    2.0
);

// mllbb
chi += std::pow(
  (vis_tot.M() - obs_map["mllbb"].mean) / obs_map["mllbb"].sigma,
  2.0
);

// mT_ttbar  (transverse mass of the (l+l-+b+bbar) system)
chi += std::pow(
  (vis_tot.Mt() - obs_map["mT_ttbar"].mean) / obs_map["mT_ttbar"].sigma,
  2.0
);