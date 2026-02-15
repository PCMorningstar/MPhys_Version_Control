#include "MyCustomFrame/MyCustomFrame.h"
#include "MyCustomFrame/Variables.h"
#include "FastFrames/DefineHelpers.h"
#include "FastFrames/UniqueSampleID.h"

ROOT::RDF::RNode MyCustomFrame::defineVariables(ROOT::RDF::RNode mainNode,
                                                const std::shared_ptr<Sample>& /*sample*/,
                                                const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom defione
  // based on the sample and the subsample
  //   sample->name(): is the name of the sample defined in the config
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  // Add a message to be sure we are running this custom class.
  //LOG(INFO) << "You are running custom class: MyCustomFrame"<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Sorting functions //////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Sorting pT columns
    LOG(INFO) << "Adding variable: el_pt_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_pt_new_NOSYS",
        ttZ::pt_order,
        {"el_pt_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_pt_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_pt_new_NOSYS",
        ttZ::pt_order,
        {"mu_pt_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_pt_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_pt_new_NOSYS",
        ttZ::pt_order,
        {"jet_pt_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting eta columns
    LOG(INFO) << "Adding variable: el_eta_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_eta_new_NOSYS",
        ttZ::pt_order,
        {"el_eta", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_eta_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_eta_new_NOSYS",
        ttZ::pt_order,
        {"mu_eta", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_eta_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_eta_new_NOSYS",
        ttZ::pt_order,
        {"jet_eta", "jet_pt_NOSYS"}
    );

    // Sorting phi columns
    LOG(INFO) << "Adding variable: el_phi_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_phi_new_NOSYS",
        ttZ::pt_order,
        {"el_phi", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_phi_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_phi_new_NOSYS",
        ttZ::pt_order,
        {"mu_phi", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_phi_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_phi_new_NOSYS",
        ttZ::pt_order,
        {"jet_phi", "jet_pt_NOSYS"}
    );

    // Sorting E columns
    LOG(INFO) << "Adding variable: el_e_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_e_new_NOSYS",
        ttZ::pt_order,
        {"el_e_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_e_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_e_new_NOSYS",
        ttZ::pt_order,
        {"mu_e_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_e_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_e_new_NOSYS",
        ttZ::pt_order,
        {"jet_e_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting "tight_ID" columns
    LOG(INFO) << "Adding variable: el_tight_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_tight_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"el_select_tight_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_tight_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_tight_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"mu_select_tight_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_jvt_new_NOSYS" << std::endl; // Jets
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_jvt_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"jet_select_baselineJvt_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting charge columns
    LOG(INFO) << "Adding variable: el_charge_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_charge_new_NOSYS",
        ttZ::pt_order,
        {"el_charge", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_charge_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_charge_new_NOSYS",
        ttZ::pt_order,
        {"mu_charge", "mu_pt_NOSYS"}
    );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.1 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: cutA1_el_pt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA1_el_pt_NOSYS",
        ttZ::cutA1_el_pt,
        {"el_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA2_el_eta_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA2_el_eta_NOSYS",
        ttZ::cutA2_el_eta,
        {"el_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA3_el_crack_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA3_el_crack_NOSYS",
        ttZ::cutA3_el_crack,
        {"el_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA4_el_tight_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA4_el_tight_NOSYS",
        ttZ::cutA4_el_tight,
        {"el_tight_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA5_mu_pt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA5_mu_pt_NOSYS",
        ttZ::cutA5_mu_pt,
        {"mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA6_mu_eta_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA6_mu_eta_NOSYS",
        ttZ::cutA6_mu_eta,
        {"mu_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA7_mu_tight_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "cutA7_mu_tight_NOSYS",
        ttZ::cutA7_mu_tight,
        {"mu_tight_new_NOSYS"}
    );

    LOG(INFO) << "Adding variable: cutA7_1_jet_jvt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "cutA7_1_jet_jvt_NOSYS",
        ttZ::cutA7_1_jet_jvt,
        {"jet_jvt_new_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A8 — Jet cleaning (jets removed if ΔR < 0.2 w.r.t electron)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: jets_clean_from_e_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "jets_clean_from_e_NOSYS",
        ttZ::jets_clean_from_e,
        {"jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "el_eta_new_NOSYS", "el_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA8_jet_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA8_jet_clean_NOSYS",
        ttZ::cutA8_jet_clean,
        {"jets_clean_from_e_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A9 — Electron cleaning (electron removed if ΔR < 0.4 to any jet)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: electron_clean_from_jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "electron_clean_from_jets_NOSYS",
        ttZ::electron_clean_from_jets,
        {"el_eta_new_NOSYS", "el_phi_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA9_el_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA9_el_clean_NOSYS",
        ttZ::cutA9_el_clean,
        {"electron_clean_from_jets_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A10 — Muon cleaning (muon removed if ΔR < 0.4 to any jet)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: muon_clean_from_jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "muon_clean_from_jets_NOSYS",
        ttZ::muon_clean_from_jets,
        {"mu_eta_new_NOSYS", "mu_phi_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA10_mu_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA10_mu_clean_NOSYS",
        ttZ::cutA10_mu_clean,
        {"muon_clean_from_jets_NOSYS"}
    );
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.2 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: cutC21_exactly2leptons_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC21_exactly2leptons_NOSYS",
        ttZ::cutC21_exactly2leptons,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC22_opposite_charge_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC22_opposite_charge_NOSYS",
        ttZ::cutC22_opposite_charge,
        {"el_charge_new_NOSYS", "mu_charge_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC23_one_el_one_mu_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC23_one_el_one_mu_NOSYS",
        ttZ::cutC23_one_el_one_mu,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC24_meemu_gt50_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC24_meemu_gt50_NOSYS",
        ttZ::cutC24_meemu_gt50,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS",
        "el_eta_new_NOSYS", "mu_eta_new_NOSYS", "el_phi_new_NOSYS", "mu_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC25_exactly2jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC25_exactly2jets_NOSYS",
        ttZ::cutC25_exactly2jets,
        {"jet_pt_new_NOSYS"}
    );
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.3 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: pairing_S6_3_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_S6_3_NOSYS",
        ttZ::compute_pairing_S6_3,
        {"el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: pairing_valid_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_valid_NOSYS",
        ttZ::cutD31_pairing_valid,
        {"pairing_S6_3_NOSYS"}   // see note below
    );
  
    LOG(INFO) << "Adding variable: pairing_mj1_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_mj1_NOSYS",
        ttZ::cutD32_mj1,
        {"pairing_S6_3_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: pairing_mj2_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_mj2_NOSYS",
        ttZ::cutD33_mj2,
        {"pairing_S6_3_NOSYS"}
    );
  
    // Altogether
    LOG(INFO) << "Adding variable: section_6_1_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_1_NOSYS",
        ttZ::section_6_1,
      {
          "el_pt_new_NOSYS",
          "el_eta_new_NOSYS",
          "el_tight_new_NOSYS",
          "mu_pt_new_NOSYS",
          "mu_eta_new_NOSYS",
          "mu_tight_new_NOSYS",
          "jet_jvt_new_NOSYS"
      }
    );
    LOG(INFO) << "Adding variable: section_6_2_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_2_NOSYS",
        ttZ::section_6_2,
        {
            "el_pt_new_NOSYS",
            "mu_pt_new_NOSYS",
            "el_eta_new_NOSYS",
            "mu_eta_new_NOSYS",
            "el_phi_new_NOSYS",
            "mu_phi_new_NOSYS",
            "el_charge_new_NOSYS",
            "mu_charge_new_NOSYS",
            "jet_pt_new_NOSYS"
        }
    );
    LOG(INFO) << "Adding variable: section_6_3_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_3_NOSYS",
        ttZ::section_6_3,
        {
          "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
          "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
          "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS"
        }
    );
  
  ////////////////////// Altogether - Selection Cuts ///////////////////////////////////////////////////
  
  LOG(INFO) << "Adding variable: selection_cuts_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "selection_cuts_NOSYS",
        ttZ::selection_cuts,
        {
          "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
          "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
          "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
          "el_charge_new_NOSYS", "mu_charge_new_NOSYS", 
          "jets_clean_from_e_NOSYS", "electron_clean_from_jets_NOSYS", "muon_clean_from_jets_NOSYS",
          "el_tight_new_NOSYS", "mu_tight_new_NOSYS", "jet_jvt_new_NOSYS"
        }
    );

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Matching Algorithms (Section – Truth & Reconstruction Pairing Diagnostic Tools)
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  // dR_truth_pairing_idx_lp_lm
  LOG(INFO) << "Adding variable: dr_truth_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dr_truth_NOSYS",
      ttZ::dr_truth,
      {}
  );

  // dR matching (calibration)
  LOG(INFO) << "Adding variable: dR_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dR_NOSYS",
      ttZ::dR_matched,
    { 
       // Truth b
       "Ttbar_History_MC_b_beforeFSR_from_t_pt",
       "Ttbar_History_MC_b_beforeFSR_from_t_eta",
       "Ttbar_History_MC_b_beforeFSR_from_t_phi",
       "Ttbar_History_MC_b_beforeFSR_from_t_m",

       // Truth bbar
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_pt",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_eta",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_phi",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_m",

       // Reco jets
       "jet_pt_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "jet_e_new_NOSYS",
        "el_pt_new_NOSYS", "el_eta_new_NOSYS", "el_phi_new_NOSYS", "el_e_new_NOSYS", "el_charge_new_NOSYS",
        "mu_pt_new_NOSYS", "mu_eta_new_NOSYS", "mu_phi_new_NOSYS", "mu_e_new_NOSYS", "mu_charge_new_NOSYS",
        "dr_truth_NOSYS"
    }
  );

  // Splitting dR outputs (plotting reasons)
  LOG(INFO) << "Adding variable: truth_m_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "truth_m_lpb_NOSYS",
      ttZ::truth_m_lpb,
    {
     "dR_NOSYS"
    }
  );
  LOG(INFO) << "Adding variable: truth_m_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "truth_m_lmbb_NOSYS",
      ttZ::truth_m_lmbb,
    {
    "dR_NOSYS"
    }
  );

  LOG(INFO) << "Adding variable: dR_truth_pairing_idx_lp_lm_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dR_truth_pairing_idx_lp_lm_NOSYS",
      ttZ::dR_truth_pairing_idx_lp_lm,
      {
        // Truth b
        "Ttbar_History_MC_b_beforeFSR_from_t_pt",
        "Ttbar_History_MC_b_beforeFSR_from_t_eta",
        "Ttbar_History_MC_b_beforeFSR_from_t_phi",
        "Ttbar_History_MC_b_beforeFSR_from_t_m",
      
        // Truth bbar
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_pt",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_eta",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_phi",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_m",
      
        // Reco jets
        "jet_pt_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "jet_e_new_NOSYS", "dr_truth_NOSYS"
      }
  );

    //chi2 Indexed
  LOG(INFO) << "Adding variable: chi_indexed_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_indexed_NOSYS",
      ttZ::chi2_pairing_min_mlb_by_charge,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
      }
  );

    // chi2 v dR enum per branch
  // l+ ↔ b
  LOG(INFO) << "Adding variable: chi_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_lpb_NOSYS",
      ttZ::chi2_vs_dR_enum_lpb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "chi_indexed_NOSYS"
      }
  );

  // l- ↔ bbar
  LOG(INFO) << "Adding variable: chi_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_lmbb_NOSYS",
      ttZ::chi2_vs_dR_enum_lmbb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "chi_indexed_NOSYS"
      }
  );

  // Jet number condition (efficiency v. jet multiplicity)
  LOG(INFO) << "Adding variable: jet_size_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_size_NOSYS",
      ttZ::jet_size,
      {
        "jet_pt_new_NOSYS"
      }
  );


  //MISMS Indexed
  LOG(INFO) << "Adding variable: misms_indexed_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_indexed_NOSYS",
      ttZ::misms_pairing_min_mlb2_by_charge,
    {
      "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
      "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
      "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
    }
  );

  // MISMS v dR enum per branch
  // l+ ↔ b
  LOG(INFO) << "Adding variable: misms_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_lpb_NOSYS",
      ttZ::misms_pairing_min_mlb2_lpb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "misms_indexed_NOSYS"
      }
  );

  // l- ↔ bbar
  LOG(INFO) << "Adding variable: misms_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_lmbb_NOSYS",
      ttZ::misms_pairing_min_mlb2_lmbb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "misms_indexed_NOSYS"
      }
  );

  return mainNode;
}

ROOT::RDF::RNode MyCustomFrame::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                      const std::shared_ptr<Sample>& /*sample*/,
                                                      const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom defione
  // based on the sample and the subsample
  //   sample->name(): is the name of the sample defined in the config
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  // Add a message to be sure we are running this custom class.
  //LOG(INFO) << "You are running custom class: MyCustomFrame"<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Sorting functions //////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Sorting pT columns
    LOG(INFO) << "Adding variable: el_pt_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_pt_new_NOSYS",
        ttZ::pt_order,
        {"el_pt_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_pt_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_pt_new_NOSYS",
        ttZ::pt_order,
        {"mu_pt_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_pt_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_pt_new_NOSYS",
        ttZ::pt_order,
        {"jet_pt_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting eta columns
    LOG(INFO) << "Adding variable: el_eta_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_eta_new_NOSYS",
        ttZ::pt_order,
        {"el_eta", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_eta_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_eta_new_NOSYS",
        ttZ::pt_order,
        {"mu_eta", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_eta_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_eta_new_NOSYS",
        ttZ::pt_order,
        {"jet_eta", "jet_pt_NOSYS"}
    );

    // Sorting phi columns
    LOG(INFO) << "Adding variable: el_phi_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_phi_new_NOSYS",
        ttZ::pt_order,
        {"el_phi", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_phi_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_phi_new_NOSYS",
        ttZ::pt_order,
        {"mu_phi", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_phi_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_phi_new_NOSYS",
        ttZ::pt_order,
        {"jet_phi", "jet_pt_NOSYS"}
    );

    // Sorting E columns
    LOG(INFO) << "Adding variable: el_e_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_e_new_NOSYS",
        ttZ::pt_order,
        {"el_e_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_e_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_e_new_NOSYS",
        ttZ::pt_order,
        {"mu_e_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_e_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_e_new_NOSYS",
        ttZ::pt_order,
        {"jet_e_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting "tight_ID" columns
    LOG(INFO) << "Adding variable: el_tight_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_tight_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"el_select_tight_NOSYS", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_tight_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_tight_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"mu_select_tight_NOSYS", "mu_pt_NOSYS"}
    );
    LOG(INFO) << "Adding variable: jet_jvt_new_NOSYS" << std::endl; // Jets
    mainNode = MainFrame::systematicDefine(mainNode,
        "jet_jvt_new_NOSYS",
        ttZ::pt_order_nonfloat,
        {"jet_select_baselineJvt_NOSYS", "jet_pt_NOSYS"}
    );

    // Sorting charge columns
    LOG(INFO) << "Adding variable: el_charge_new_NOSYS" << std::endl; // Electrons
    mainNode = MainFrame::systematicDefine(mainNode,
        "el_charge_new_NOSYS",
        ttZ::pt_order,
        {"el_charge", "el_pt_NOSYS"} // {to-be-sorted, sorted idx reference}
    );
    LOG(INFO) << "Adding variable: mu_charge_new_NOSYS" << std::endl; // Muons
    mainNode = MainFrame::systematicDefine(mainNode,
        "mu_charge_new_NOSYS",
        ttZ::pt_order,
        {"mu_charge", "mu_pt_NOSYS"}
    );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.1 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: cutA1_el_pt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA1_el_pt_NOSYS",
        ttZ::cutA1_el_pt,
        {"el_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA2_el_eta_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA2_el_eta_NOSYS",
        ttZ::cutA2_el_eta,
        {"el_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA3_el_crack_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA3_el_crack_NOSYS",
        ttZ::cutA3_el_crack,
        {"el_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA4_el_tight_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA4_el_tight_NOSYS",
        ttZ::cutA4_el_tight,
        {"el_tight_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA5_mu_pt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA5_mu_pt_NOSYS",
        ttZ::cutA5_mu_pt,
        {"mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA6_mu_eta_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA6_mu_eta_NOSYS",
        ttZ::cutA6_mu_eta,
        {"mu_eta_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA7_mu_tight_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "cutA7_mu_tight_NOSYS",
        ttZ::cutA7_mu_tight,
        {"mu_tight_new_NOSYS"}
    );

    LOG(INFO) << "Adding variable: cutA7_1_jet_jvt_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "cutA7_1_jet_jvt_NOSYS",
        ttZ::cutA7_1_jet_jvt,
        {"jet_jvt_new_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A8 — Jet cleaning (jets removed if ΔR < 0.2 w.r.t electron)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: jets_clean_from_e_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "jets_clean_from_e_NOSYS",
        ttZ::jets_clean_from_e,
        {"jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "el_eta_new_NOSYS", "el_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA8_jet_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA8_jet_clean_NOSYS",
        ttZ::cutA8_jet_clean,
        {"jets_clean_from_e_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A9 — Electron cleaning (electron removed if ΔR < 0.4 to any jet)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: electron_clean_from_jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "electron_clean_from_jets_NOSYS",
        ttZ::electron_clean_from_jets,
        {"el_eta_new_NOSYS", "el_phi_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA9_el_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA9_el_clean_NOSYS",
        ttZ::cutA9_el_clean,
        {"electron_clean_from_jets_NOSYS"}
    );
  
    // -----------------------------------------------------------------------------
    // A10 — Muon cleaning (muon removed if ΔR < 0.4 to any jet)
    // -----------------------------------------------------------------------------
  
    LOG(INFO) << "Adding variable: muon_clean_from_jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "muon_clean_from_jets_NOSYS",
        ttZ::muon_clean_from_jets,
        {"mu_eta_new_NOSYS", "mu_phi_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutA10_mu_clean_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutA10_mu_clean_NOSYS",
        ttZ::cutA10_mu_clean,
        {"muon_clean_from_jets_NOSYS"}
    );
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.2 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: cutC21_exactly2leptons_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC21_exactly2leptons_NOSYS",
        ttZ::cutC21_exactly2leptons,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC22_opposite_charge_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC22_opposite_charge_NOSYS",
        ttZ::cutC22_opposite_charge,
        {"el_charge_new_NOSYS", "mu_charge_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC23_one_el_one_mu_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC23_one_el_one_mu_NOSYS",
        ttZ::cutC23_one_el_one_mu,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC24_meemu_gt50_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC24_meemu_gt50_NOSYS",
        ttZ::cutC24_meemu_gt50,
        {"el_pt_new_NOSYS", "mu_pt_new_NOSYS",
        "el_eta_new_NOSYS", "mu_eta_new_NOSYS", "el_phi_new_NOSYS", "mu_phi_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: cutC25_exactly2jets_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(mainNode,
        "cutC25_exactly2jets_NOSYS",
        ttZ::cutC25_exactly2jets,
        {"jet_pt_new_NOSYS"}
    );
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Section 6.3 ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    LOG(INFO) << "Adding variable: pairing_S6_3_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_S6_3_NOSYS",
        ttZ::compute_pairing_S6_3,
        {"el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: pairing_valid_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_valid_NOSYS",
        ttZ::cutD31_pairing_valid,
        {"pairing_S6_3_NOSYS"}   // see note below
    );
  
    LOG(INFO) << "Adding variable: pairing_mj1_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_mj1_NOSYS",
        ttZ::cutD32_mj1,
        {"pairing_S6_3_NOSYS"}
    );
  
    LOG(INFO) << "Adding variable: pairing_mj2_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "pairing_mj2_NOSYS",
        ttZ::cutD33_mj2,
        {"pairing_S6_3_NOSYS"}
    );
  
    // Altogether
    LOG(INFO) << "Adding variable: section_6_1_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_1_NOSYS",
        ttZ::section_6_1,
      {
          "el_pt_new_NOSYS",
          "el_eta_new_NOSYS",
          "el_tight_new_NOSYS",
          "mu_pt_new_NOSYS",
          "mu_eta_new_NOSYS",
          "mu_tight_new_NOSYS",
          "jet_jvt_new_NOSYS"
      }
    );
    LOG(INFO) << "Adding variable: section_6_2_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_2_NOSYS",
        ttZ::section_6_2,
        {
            "el_pt_new_NOSYS",
            "mu_pt_new_NOSYS",
            "el_eta_new_NOSYS",
            "mu_eta_new_NOSYS",
            "el_phi_new_NOSYS",
            "mu_phi_new_NOSYS",
            "el_charge_new_NOSYS",
            "mu_charge_new_NOSYS",
            "jet_pt_new_NOSYS"
        }
    );
    LOG(INFO) << "Adding variable: section_6_3_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "section_6_3_NOSYS",
        ttZ::section_6_3,
        {
          "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
          "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
          "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS"
        }
    );
  
  ////////////////////// Altogether - Selection Cuts ///////////////////////////////////////////////////
  
  LOG(INFO) << "Adding variable: selection_cuts_NOSYS" << std::endl;
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "selection_cuts_NOSYS",
        ttZ::selection_cuts,
        {
          "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
          "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
          "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
          "el_charge_new_NOSYS", "mu_charge_new_NOSYS", 
          "jets_clean_from_e_NOSYS", "electron_clean_from_jets_NOSYS", "muon_clean_from_jets_NOSYS",
          "el_tight_new_NOSYS", "mu_tight_new_NOSYS", "jet_jvt_new_NOSYS"
        }
    );

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Matching Algorithms (Section – Truth & Reconstruction Pairing Diagnostic Tools)
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  // dR_truth_pairing_idx_lp_lm
  LOG(INFO) << "Adding variable: dr_truth_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dr_truth_NOSYS",
      ttZ::dr_truth,
      {}
  );

  // dR matching (calibration)
  LOG(INFO) << "Adding variable: dR_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dR_NOSYS",
      ttZ::dR_matched,
    { 
       // Truth b
       "Ttbar_History_MC_b_beforeFSR_from_t_pt",
       "Ttbar_History_MC_b_beforeFSR_from_t_eta",
       "Ttbar_History_MC_b_beforeFSR_from_t_phi",
       "Ttbar_History_MC_b_beforeFSR_from_t_m",

       // Truth bbar
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_pt",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_eta",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_phi",
       "Ttbar_History_MC_bbar_beforeFSR_from_tbar_m",

       // Reco jets
       "jet_pt_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "jet_e_new_NOSYS",
        "el_pt_new_NOSYS", "el_eta_new_NOSYS", "el_phi_new_NOSYS", "el_e_new_NOSYS", "el_charge_new_NOSYS",
        "mu_pt_new_NOSYS", "mu_eta_new_NOSYS", "mu_phi_new_NOSYS", "mu_e_new_NOSYS", "mu_charge_new_NOSYS",
        "dr_truth_NOSYS"
    }
  );

  // Splitting dR outputs (plotting reasons)
  LOG(INFO) << "Adding variable: truth_m_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "truth_m_lpb_NOSYS",
      ttZ::truth_m_lpb,
    {
     "dR_NOSYS"
    }
  );
  LOG(INFO) << "Adding variable: truth_m_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "truth_m_lmbb_NOSYS",
      ttZ::truth_m_lmbb,
    {
    "dR_NOSYS"
    }
  );

  LOG(INFO) << "Adding variable: dR_truth_pairing_idx_lp_lm_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "dR_truth_pairing_idx_lp_lm_NOSYS",
      ttZ::dR_truth_pairing_idx_lp_lm,
      {
        // Truth b
        "Ttbar_History_MC_b_beforeFSR_from_t_pt",
        "Ttbar_History_MC_b_beforeFSR_from_t_eta",
        "Ttbar_History_MC_b_beforeFSR_from_t_phi",
        "Ttbar_History_MC_b_beforeFSR_from_t_m",
      
        // Truth bbar
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_pt",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_eta",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_phi",
        "Ttbar_History_MC_bbar_beforeFSR_from_tbar_m",
      
        // Reco jets
        "jet_pt_new_NOSYS", "jet_eta_new_NOSYS", "jet_phi_new_NOSYS", "jet_e_new_NOSYS", "dr_truth_NOSYS"
      }
  );

    //chi2 Indexed
  LOG(INFO) << "Adding variable: chi_indexed_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_indexed_NOSYS",
      ttZ::chi2_pairing_min_mlb_by_charge,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
      }
  );

    // chi2 v dR enum per branch
  // l+ ↔ b
  LOG(INFO) << "Adding variable: chi_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_lpb_NOSYS",
      ttZ::chi2_vs_dR_enum_lpb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "chi_indexed_NOSYS"
      }
  );

  // l- ↔ bbar
  LOG(INFO) << "Adding variable: chi_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "chi_lmbb_NOSYS",
      ttZ::chi2_vs_dR_enum_lmbb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "chi_indexed_NOSYS"
      }
  );

  // Jet number condition (efficiency v. jet multiplicity)
  LOG(INFO) << "Adding variable: jet_size_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_size_NOSYS",
      ttZ::jet_size,
      {
        "jet_pt_new_NOSYS"
      }
  );

  //MISMS Indexed
  LOG(INFO) << "Adding variable: misms_indexed_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_indexed_NOSYS",
      ttZ::misms_pairing_min_mlb2_by_charge,
    {
      "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
      "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
      "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
    }
  );

  // MISMS v dR enum per branch
  // l+ ↔ b
  LOG(INFO) << "Adding variable: misms_lpb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_lpb_NOSYS",
      ttZ::misms_pairing_min_mlb2_lpb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "misms_indexed_NOSYS"
      }
  );

  // l- ↔ bbar
  LOG(INFO) << "Adding variable: misms_lmbb_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "misms_lmbb_NOSYS",
      ttZ::misms_pairing_min_mlb2_lmbb,
      {
        "dR_truth_pairing_idx_lp_lm_NOSYS",
        "misms_indexed_NOSYS"
      }
  );

  return mainNode;
}

ROOT::RDF::RNode MyCustomFrame::defineVariablesTruth(ROOT::RDF::RNode node,
                                                     const std::string& /*sample*/,
                                                     const std::shared_ptr<Sample>& /*sample*/,
                                                     const UniqueSampleID& /*sampleID*/) {

  return node;

}

ROOT::RDF::RNode MyCustomFrame::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}