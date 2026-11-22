  for (m in 1:N_market) {
    int start = market_obs_start[m];
    int L = market_obs_n[m];
    array[L] int rows = market_obs_order[start:(start + L - 1)];
    int G = quality_scope == 2 ? quality_block_n[m] : N_quality;
    market_panel[m] = panel_design(quality_local[rows], G, relative_cell_local[rows],
      n_mky_per_market[m], margin_obs[rows], rho_active);
  }
  // A shared quality effect spans markets. Check the entire incidence matrix
  // and evaluate one joint density in the original global coordinates.
  if (quality_scope == 1 && N_market > 1 && quality_parameterization == 2)
    global_panel = panel_design(obs_to_quality, N_quality, market_year_idx,
      N_market_year, margin_obs, rho_active);
