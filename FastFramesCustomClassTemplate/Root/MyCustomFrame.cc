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
  

  // New attempt at Chi2 analysis (no external dR func idx needed BUT means and stand devs!)
  // Ordering truth jet idx & candidates
  // b-jet
  LOG(INFO) << "Adding variable: bjet_new_idx_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "bjet_new_idx_NOSYS",
      ttZ::b_selector,
      {"event_jet_truth_idx"} 
  );
  LOG(INFO) << "Adding variable: bjet_new_candicate_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "bjet_new_candicate_NOSYS",
      ttZ::b_selector,
      {"event_jet_truth_candidates"}
  );
  // bb-jet
  LOG(INFO) << "Adding variable: bbarjet_new_idx_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "bbarjet_new_idx_NOSYS",
      ttZ::bbar_selector,
      {"event_jet_truth_idx"} 
  );
  LOG(INFO) << "Adding variable: bbarjet_new_candicate_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "bbarjet_new_candicate_NOSYS",
      ttZ::bbar_selector,
      {"event_jet_truth_candidates"}
  );

  // New detailed truth distributions
  // 2 jets
  LOG(INFO) << "Adding variable: new_detailed_truth_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_detailed_truth_NOSYS",
      ttZ::new_detailed_truth,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "met_met_NOSYS", "met_phi_NOSYS", "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS"
      }
  );

  // truth distribution extraction
  LOG(INFO) << "Adding variable: new_truth_mlpb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mlpb_NOSYS",
      ttZ::new_truth_mlpb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mlmbb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mlmbb_NOSYS",
      ttZ::new_truth_mlmbb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_pTdiff_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_pTdiff_NOSYS",
      ttZ::new_truth_pTdiff,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_sum_deltaR_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_sum_deltaR_NOSYS",
      ttZ::new_truth_sum_deltaR,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mllbb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mllbb_NOSYS",
      ttZ::new_truth_mllbb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mT_ttbar_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mT_ttbar_NOSYS",
      ttZ::new_truth_mT_ttbar,
      {"new_detailed_truth_NOSYS"}
  );

  
  // Chi2 
  LOG(INFO) << "Adding variable: new_chi_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_chi_indexed_jets_NOSYS",
      ttZ::new_chi_indexed,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "met_met_NOSYS", "met_phi_NOSYS", "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
      }
  );
  // MISMS
  LOG(INFO) << "Adding variable: new_misms_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_misms_indexed_jets_NOSYS",
      ttZ::new_misms_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
      }
  );
  // MDRS
  LOG(INFO) << "Adding variable: new_mdrs_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_mdrs_indexed_jets_NOSYS",
      ttZ::new_MDRS_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
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

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Selections for cut-flow /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: electron_selections_paper_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "electron_selections_paper_NOSYS",
      ttZ::electron_selections_paper,
      {"el_pt_new_NOSYS","el_eta_new_NOSYS","el_tight_new_NOSYS", "electron_clean_from_jets_NOSYS"}
  );
  LOG(INFO) << "Adding variable: muon_selections_paper_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "muon_selections_paper_NOSYS",
      ttZ::muon_selections_paper,
      {"mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_tight_new_NOSYS", "muon_clean_from_jets_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_selections_paper_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "jet_selections_paper_NOSYS",
      ttZ::jet_selections_paper,
      {"jet_jvt_new_NOSYS", "jets_clean_from_e_NOSYS", "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS",
      "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS",
      "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS"}
  );
  LOG(INFO) << "Adding variable: dilepton_selections_paper_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "dilepton_selections_paper_NOSYS",
      ttZ::dilepton_selections_paper,
      {"el_pt_new_NOSYS","mu_pt_new_NOSYS","el_eta_new_NOSYS","mu_eta_new_NOSYS","el_phi_new_NOSYS",
        "mu_phi_new_NOSYS","el_charge_new_NOSYS","mu_charge_new_NOSYS"}
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Matching Algorithms (Section – Truth & Reconstruction Pairing Diagnostic Tools)
  //////////////////////////////////////////////////////////////////////////////////////////////////////


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
  

  // New attempt at Chi2 analysis (no external dR func idx needed BUT means and stand devs!)
  // Ordering truth jet idx & candidates
  // b-jet
  LOG(INFO) << "Adding variable: bjet_new_idx_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "bjet_new_idx_NOSYS",
      ttZ::b_selector,
      {"event_jet_truth_idx"} 
  );
  LOG(INFO) << "Adding variable: bjet_new_candicate_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "bjet_new_candicate_NOSYS",
      ttZ::b_selector,
      {"event_jet_truth_candidates"}
  );
  // bb-jet
  LOG(INFO) << "Adding variable: bbarjet_new_idx_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(mainNode,
      "bbarjet_new_idx_NOSYS",
      ttZ::bbar_selector,
      {"event_jet_truth_idx"} 
  );
  LOG(INFO) << "Adding variable: bbarjet_new_candicate_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "bbarjet_new_candicate_NOSYS",
      ttZ::bbar_selector,
      {"event_jet_truth_candidates"}
  );

  // New detailed truth distributions
  // 2 jets
  LOG(INFO) << "Adding variable: new_detailed_truth_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_detailed_truth_NOSYS",
      ttZ::new_detailed_truth,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "met_met_NOSYS", "met_phi_NOSYS", "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS"
      }
  );

  // truth distribution extraction
  LOG(INFO) << "Adding variable: new_truth_mlpb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mlpb_NOSYS",
      ttZ::new_truth_mlpb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mlmbb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mlmbb_NOSYS",
      ttZ::new_truth_mlmbb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_pTdiff_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_pTdiff_NOSYS",
      ttZ::new_truth_pTdiff,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_sum_deltaR_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_sum_deltaR_NOSYS",
      ttZ::new_truth_sum_deltaR,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mllbb_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mllbb_NOSYS",
      ttZ::new_truth_mllbb,
      {"new_detailed_truth_NOSYS"}
  );
  LOG(INFO) << "Adding variable: new_truth_mT_ttbar_NOSYS" << std::endl; 
  mainNode = MainFrame::systematicDefine(mainNode,
      "new_truth_mT_ttbar_NOSYS",
      ttZ::new_truth_mT_ttbar,
      {"new_detailed_truth_NOSYS"}
  );

  
  // Chi2 
  LOG(INFO) << "Adding variable: new_chi_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_chi_indexed_jets_NOSYS",
      ttZ::new_chi_indexed,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "met_met_NOSYS", "met_phi_NOSYS", "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
      }
  );
  // MISMS
  LOG(INFO) << "Adding variable: new_misms_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_misms_indexed_jets_NOSYS",
      ttZ::new_misms_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
      }
  );

  // MDRS
  LOG(INFO) << "Adding variable: new_mdrs_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "new_mdrs_indexed_jets_NOSYS",
      ttZ::new_MDRS_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "bjet_new_idx_NOSYS", "bbarjet_new_idx_NOSYS",
        "bjet_new_candicate_NOSYS", "bbarjet_new_candicate_NOSYS"
      }
  );

  // ============================================================
  // Raw pairing algorithms in the (l+, l-) basis
  // Returns:
  //  best_i  -> index of the best jet for l+ (for b)
  //  best_j  -> index of the best jet for l- (for bbar)
  // -1      -> no valid pairing
  // ============================================================
  LOG(INFO) << "Adding variable: raw_chi2_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_chi2_indexed_jets_NOSYS",
      ttZ::raw_chi2_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
      }
  );
  LOG(INFO) << "Adding variable: raw_mdrs_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_mdrs_indexed_jets_NOSYS",
      ttZ::raw_MDRS_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
      }
  );
  LOG(INFO) << "Adding variable: raw_misms_indexed_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_misms_indexed_jets_NOSYS",
      ttZ::raw_MISMS_pairing,
      {
        "jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"
      }
  );

  // Ordered jet truth flavour vector
  LOG(INFO) << "Adding variable: ordered_jet_truth_flavour_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ordered_jet_truth_flavour_NOSYS",
      ttZ::pt_order_int,
      {"jet_TruthFlavour", "jet_pt_NOSYS"}
  );

  // Ordered jet truth flavour vector - extended
  LOG(INFO) << "Adding variable: ordered_jet_truth_flavour_extended_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ordered_jet_truth_flavour_extended_NOSYS",
      ttZ::pt_order_int,
      {"jet_TruthFlavouExtended", "jet_pt_NOSYS"}
  );

  // Raw pairing indices
  LOG(INFO) << "Adding variable: raw_chi2_pairing_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_chi2_pairing_NOSYS",
      ttZ::raw_chi2_pairing,
      {"jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"}
  );
  LOG(INFO) << "Adding variable: raw_mdrs_pairing_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_mdrs_pairing_NOSYS",
      ttZ::raw_MDRS_pairing,
      {"jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"}
  );
  LOG(INFO) << "Adding variable: raw_misms_pairing_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_misms_pairing_NOSYS",
      ttZ::raw_MISMS_pairing,
      {"jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS"}
  );

  // Chi2 essential outputs - algorithm analysis
  LOG(INFO) << "Adding variable: raw_chi2_minval_truthall_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "raw_chi2_minval_truthall_NOSYS",
      ttZ::raw_chi2_minval_truthall,
      {"jet_pt_new_NOSYS","jet_eta_new_NOSYS","jet_phi_new_NOSYS","jet_e_new_NOSYS",
        "el_pt_new_NOSYS","el_eta_new_NOSYS","el_phi_new_NOSYS","el_e_new_NOSYS","el_charge_new_NOSYS",
        "mu_pt_new_NOSYS","mu_eta_new_NOSYS","mu_phi_new_NOSYS","mu_e_new_NOSYS","mu_charge_new_NOSYS",
        "event_jet_truth_idx", "event_jet_truth_candidates"}
    );

  // SV1 - ordering - pt_order_float
  LOG(INFO) << "Adding variable: sv1_ordered_jets_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv1_ordered_jets_NOSYS",
      ttZ::pt_order,
      {"jet_SV1_masssvx", "jet_pt_new_NOSYS"}
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Jet pT selection (changed in yaml) /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: jet_pt_region_0to30_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_0to30_GeV_NOSYS",
      ttZ::jet_pt_region_0to30_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_30to60_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_30to60_GeV_NOSYS",
      ttZ::jet_pt_region_30to60_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_60to90_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_60to90_GeV_NOSYS",
      ttZ::jet_pt_region_60to90_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_90to120_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_90to120_GeV_NOSYS",
      ttZ::jet_pt_region_90to120_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_120to150_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_120to150_GeV_NOSYS",
      ttZ::jet_pt_region_120to150_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_150to180_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_150to180_GeV_NOSYS",
      ttZ::jet_pt_region_150to180_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_180to210_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_180to210_GeV_NOSYS",
      ttZ::jet_pt_region_180to210_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_210to240_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_210to240_GeV_NOSYS",
      ttZ::jet_pt_region_210to240_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_240to270_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_240to270_GeV_NOSYS",
      ttZ::jet_pt_region_240to270_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_270to300_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_270to300_GeV_NOSYS",
      ttZ::jet_pt_region_270to300_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_300to360_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(   
      mainNode,
      "jet_pt_region_300to360_GeV_NOSYS",
      ttZ::jet_pt_region_300to360_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: jet_pt_region_360to900_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_region_360to900_GeV_NOSYS",
      ttZ::jet_pt_region_360to900_GeV,
      {"jet_pt_new_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// SV Invariant Mass selection (changed in yaml) /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: sv_invariant_mass_region_0to0point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_0to0point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_0to0point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_0point5to1_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_0point5to1_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_0point5to1_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_1to1point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_1to1point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_1to1point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_1point5to2_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_1point5to2_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_1point5to2_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_2to2point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_2to2point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_2to2point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_2point5to3_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_2point5to3_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_2point5to3_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_3to3point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_3to3point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_3to3point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_3point5to4_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_3point5to4_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_3point5to4_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_4to4point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_4to4point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_4to4point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_4point5to5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_4point5to5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_4point5to5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_5to5point5_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_5to5point5_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_5to5point5_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
  );
  LOG(INFO) << "Adding variable: sv_invariant_mass_region_5point5to6_GeV_NOSYS" << std::endl;
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sv_invariant_mass_region_5point5to6_GeV_NOSYS",
      ttZ::sv_invariant_mass_region_5point5to6_GeV,
      {"sv1_ordered_jets_NOSYS", "raw_chi2_minval_truthall_NOSYS"}
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