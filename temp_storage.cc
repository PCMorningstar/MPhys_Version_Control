// New Chi2
if (jet_size == 2) {
    obs_map["mlb_plus"]   = {98.26f, 30.42f};
    obs_map["mlb_minus"]  = {98.34f, 30.49f};
    obs_map["pTdiff"]     = {0.00f, 50.81f};
    obs_map["sum_deltaR"] = {3.60f, 1.43f};
    obs_map["mllbb"]      = {303.21f, 85.56f};
    obs_map["mT_ttbar"]   = {377.03f, 94.96f};}
  else if (jet_size == 3) {
    obs_map["mlb_plus"]   = {97.23f, 31.28f};
    obs_map["mlb_minus"]  = {97.32f, 31.34f};
    obs_map["pTdiff"]     = {0.00f, 59.38f};
    obs_map["sum_deltaR"] = {3.46f, 1.39f};
    obs_map["mllbb"]      = {306.21f, 89.54f};
    obs_map["mT_ttbar"]   = {381.30f, 99.61f};}
  else if (jet_size == 4) {
    obs_map["mlb_plus"]   = {96.46f, 32.05f};
    obs_map["mlb_minus"]  = {96.84f, 31.92f};
    obs_map["pTdiff"]     = {0.00f, 68.03f};
    obs_map["sum_deltaR"] = {3.35f, 1.35f};
    obs_map["mllbb"]      = {311.15f, 94.40f};
    obs_map["mT_ttbar"]   = {388.19f, 105.37f};}
  else if (jet_size == 5) {
    obs_map["mlb_plus"]   = {96.53f, 32.30f};
    obs_map["mlb_minus"]  = {96.69f, 32.26f};
    obs_map["pTdiff"]     = {0.00f, 77.32f};
    obs_map["sum_deltaR"] = {3.26f, 1.32f};
    obs_map["mllbb"]      = {317.27f, 99.47f};
    obs_map["mT_ttbar"]   = {395.96f, 110.84f};}
  else if (jet_size == 6) {
    obs_map["mlb_plus"]   = {96.40f, 32.60f};
    obs_map["mlb_minus"]  = {96.45f, 32.56f};
    obs_map["pTdiff"]     = {0.00f, 87.36f};
    obs_map["sum_deltaR"] = {3.17f, 1.30f};
    obs_map["mllbb"]      = {322.65f, 104.40f};
    obs_map["mT_ttbar"]   = {403.99f, 117.36f};}
  else if (jet_size == 7) {
    obs_map["mlb_plus"]   = {96.54f, 32.92f};
    obs_map["mlb_minus"]  = {96.44f, 33.14f};
    obs_map["pTdiff"]     = {0.00f, 94.02f};
    obs_map["sum_deltaR"] = {3.08f, 1.28f};
    obs_map["mllbb"]      = {329.41f, 108.16f};
    obs_map["mT_ttbar"]   = {413.13f, 121.96f};}
  else if (jet_size == 8) {
    obs_map["mlb_plus"]   = {97.62f, 33.56f};
    obs_map["mlb_minus"]  = {96.62f, 33.30f};
    obs_map["pTdiff"]     = {0.00f, 108.20f};
    obs_map["sum_deltaR"] = {3.04f, 1.27f};
    obs_map["mllbb"]      = {338.37f, 117.51f};
    obs_map["mT_ttbar"]   = {423.95f, 131.37f};}
  else if (jet_size == 9) {
    obs_map["mlb_plus"]   = {97.62f, 33.56f};
    obs_map["mlb_minus"]  = {96.62f, 33.30f};
    obs_map["pTdiff"]     = {0.00f, 108.20f};
    obs_map["sum_deltaR"] = {3.04f, 1.27f};
    obs_map["mllbb"]      = {338.37f, 117.51f};
    obs_map["mT_ttbar"]   = {423.95f, 131.37f};}
  else if (jet_size == 10) {
    obs_map["mlb_plus"]   = {101.64f, 29.87f};
    obs_map["mlb_minus"]  = {99.40f, 36.47f};
    obs_map["pTdiff"]     = {0.00f, 124.86f};
    obs_map["sum_deltaR"] = {2.94f, 1.26f};
    obs_map["mllbb"]      = {357.22f, 131.86f};
    obs_map["mT_ttbar"]   = {450.76f, 151.15f};}

// Old Chi2
if (jet_size == 2) {
    mean_m_lpb  = 98.10f;
    mean_m_lmbb = 98.27f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 3.58f;
    mean_mllbb = 302.99f;
    mean_mT_ttbar = 376.69f;

    sigma_m_lpb = 29.89f;
    sigma_m_lmbb =  29.92f;
    sigma_pTdiff =  50.74f;
    sigma_sum_deltaR = 1.42f;
    sigma_mllbb =  85.26f;
    sigma_mT_ttbar =  94.61f;
  }
  else if (jet_size == 3) {
    mean_m_lpb  = 97.47f;
    mean_m_lmbb = 97.54f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 3.42f;
    mean_mllbb = 306.08f;
    mean_mT_ttbar = 381.11f;

    sigma_m_lpb = 30.09f;
    sigma_m_lmbb =  30.17f;
    sigma_pTdiff =  59.51f;
    sigma_sum_deltaR = 1.36f;
    sigma_mllbb =  88.19f;
    sigma_mT_ttbar =  98.19f;
  }
  else if (jet_size == 4) {
    mean_m_lpb  = 96.85f;
    mean_m_lmbb = 97.18f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 3.28f;
    mean_mllbb = 310.44f;
    mean_mT_ttbar = 387.35f;

    sigma_m_lpb = 30.41f;
    sigma_m_lmbb =  30.35f;
    sigma_pTdiff =  68.50f;
    sigma_sum_deltaR = 1.30f;
    sigma_mllbb =  92.31f;
    sigma_mT_ttbar =  102.98f;
  }
  else if (jet_size == 5) {
    mean_m_lpb  = 96.73f;
    mean_m_lmbb = 96.97f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 3.16f;
    mean_mllbb = 315.27f;
    mean_mT_ttbar = 393.89f;

    sigma_m_lpb = 30.43f;
    sigma_m_lmbb =  30.39f;
    sigma_pTdiff =  78.57f;
    sigma_sum_deltaR = 1.25f;
    sigma_mllbb =  95.95f;
    sigma_mT_ttbar =  107.37f;
  }
  else if (jet_size == 6) {
    mean_m_lpb  = 96.63f;
    mean_m_lmbb = 96.64f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 3.04f;
    mean_mllbb = 319.54f;
    mean_mT_ttbar = 399.76f;

    sigma_m_lpb =  30.52f;
    sigma_m_lmbb =  30.59f;
    sigma_pTdiff =  89.33f;
    sigma_sum_deltaR = 1.21f;
    sigma_mllbb =  99.70f;
    sigma_mT_ttbar =  111.71f;
  }
  else if (jet_size == 7) {
    mean_m_lpb  = 96.47f;
    mean_m_lmbb = 96.51f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 2.92f;
    mean_mllbb = 323.03f;
    mean_mT_ttbar = 405.69f;

    sigma_m_lpb =  30.60f;
    sigma_m_lmbb =  30.81f;
    sigma_pTdiff =  97.78f;
    sigma_sum_deltaR = 1.17f;
    sigma_mllbb =  101.49f;
    sigma_mT_ttbar =  114.18f;
  }
  else if (jet_size == 8) {
    mean_m_lpb  = 96.75f;
    mean_m_lmbb = 96.34f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 2.83f;
    mean_mllbb = 328.79f;
    mean_mT_ttbar = 413.22f;

    sigma_m_lpb =  30.58f;
    sigma_m_lmbb =  31.34f;
    sigma_pTdiff =  109.03f;
    sigma_sum_deltaR = 1.12f;
    sigma_mllbb =  109.02f;
    sigma_mT_ttbar =  120.78f;
  }
  else if (jet_size == 9) {
    mean_m_lpb  = 98.13f;
    mean_m_lmbb = 98.47f;
    mean_pTdiff =  0.00f;
    mean_sum_deltaR = 2.75f;
    mean_mllbb = 326.38f;
    mean_mT_ttbar = 410.57f;

    sigma_m_lpb =  30.78f;
    sigma_m_lmbb =  31.79f;
    sigma_pTdiff =  110.81f;
    sigma_sum_deltaR = 1.08f;
    sigma_mllbb =  100.17f;
    sigma_mT_ttbar =  116.48f;
  }
  else if (jet_size == 10) {
    mean_m_lpb  = 100.64f;
    mean_m_lmbb = 96.27f;
    mean_pTdiff = 0.00f;
    mean_sum_deltaR = 2.55f;
    mean_mllbb = 343.23f;
    mean_mT_ttbar = 436.73f;

    sigma_m_lpb =  25.88f;
    sigma_m_lmbb =  32.44f;
    sigma_pTdiff =  130.77f;
    sigma_sum_deltaR = 1.02f;
    sigma_mllbb =  125.58f;
    sigma_mT_ttbar =  147.79f;
  }