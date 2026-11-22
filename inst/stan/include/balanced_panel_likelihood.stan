  // One product per ranked firm and period. The R dispatcher verifies the
  // rectangular quality incidence. Indices are built once in transformed data.
  for (m in 1:N_market) {
    int T = market_mky_n[m];
    int J = firm_market_n[m];
    int first_t = market_mky_start[m];
    int u = mky_to_partition_unit[market_mky_order[first_t]];
    int nk = partition_unit_count[u];
    vector[nk] score;
    real v = square(sigma_share_market[m]) * (1 - square(rho_market[m]));
    real mean_var = v + T * square(sigma_b_strat_economic);
    real slope = rho_market[m] * sigma_share_market[m] / sigma_margin_market[m];
    // Each of the J-1 orthonormal product contrasts has T-1 time contrasts
    // and one normalized time average.
    real log_norm = (J - 1) * (T * log(2 * pi()) + (T - 1) * log(v) + log(mean_var));
    matrix[J, T] share_resid_panel;
    matrix[J, T] follower_margin_resid;
    matrix[J, T] leader_margin_resid;
    matrix[J, T] shares;
    matrix[J + 1, T] follower_sum = rep_matrix(0.0, J + 1, T);
    for (t in 1:T) {
      int mky = market_mky_order[first_t + t - 1];
      for (back in 1:J) {
        int j = J + 1 - back;
        int fj = firm_mky_start[mky] + j - 1;
        int i = panel_firm_obs[fj];
        real S = firm_mky_shareIn[fj] * (1 - s0_prob[mky]);
        real ap = alpha_structural_market[m] * price[i];
        real strategic_h = supply_model == 1 ? 1 / (1 - S) : 1 + S / s0_prob[mky];
        shares[j, t] = S;
        share_resid_panel[j, t] = log_shareIn[i] - relative_share_predictor[i];
        leader_margin_resid[j, t] = log_marginNormInv[i] - log(ap / strategic_h + 1e-9);
        if (structure == 1) {
          follower_margin_resid[j, t] = log_marginNormInv[i] - log(ap + 1e-9);
        } else {
          follower_margin_resid[j, t] = leader_margin_resid[j, t];
          follower_sum[j, t] = follower_sum[j + 1, t]
            + (supply_model == 1 ? square(S) / (1 - S + square(S)) : S);
        }
      }
    }
    for (kk in 1:nk) {
      int cc = partition_unit_start[u] + kk - 1;
      int k = partition_component_k[cc];
      matrix[J, T] residual = follower_margin_resid;
      matrix[J, T] z;
      vector[J] time_mean;
      for (t in 1:T) {
        int mky = market_mky_order[first_t + t - 1];
        if (k > 0) {
          if (structure == 1) residual[1:k, t] = leader_margin_resid[1:k, t];
          else for (j in 1:k) {
            int i = panel_firm_obs[firm_mky_start[mky] + j - 1];
            real S = shares[j, t];
            real f = follower_sum[k + 1, t];
            real h = supply_model == 1 ? 1 / (1 - S / (1 - f))
              : 1 + S / (s0_prob[mky] + f);
            residual[j, t] = log_marginNormInv[i]
              - log(price[i] * alpha_structural_market[m] / h + 1e-9);
          }
        }
        z[, t] = share_resid_panel[, t] - slope * residual[, t];
        z[, t] -= mean(z[, t]);
      }
      time_mean = z * rep_vector(1.0 / T, T);
      // The unconditional share term cancels exactly in the joint mixture.
      score[kk] = partition_prior_probability[cc] > 0
        ? partition_prior_logprob[cc] + normal_lpdf(to_vector(residual) | 0, sigma_margin_market[m])
          - 0.5 * (log_norm + dot_self(to_vector(z - rep_matrix(time_mean, T))) / v
                   + T * dot_self(time_mean) / mean_var)
        : negative_infinity();
    }
    target += log_sum_exp(score);
  }
