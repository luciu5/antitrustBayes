// Bayesian CES Demand + Supply Model for Merger Simulation
// Hierarchical model with market-level random effects
// Supports 4 supply models: Bertrand, Monopolistic Competition, Cournot, 2nd Score Auction
// CES: revenue shares, log-prices, gamma = elasticity of substitution

functions {
  real bounded_inv_logit(real eta) {
    real eps = 1e-6;
    return eps + (1.0 - 2.0 * eps) * inv_logit(eta);
  }

  real soft_hmt_penalty(real agg_elast, real max_elast, real scale) {
    if (agg_elast > max_elast) {
      return normal_lpdf(agg_elast | max_elast, scale) - normal_lpdf(max_elast | max_elast, scale);
    }
    return 0.0;
  }

  real partial_sum_fn(array[] int slice_idx, int start, int end,
                      vector lernerInv, vector revShare, vector firm_revShare, vector log_revShare,
                      vector log_lernerInv, vector log_price_std,
                      array[] int market, array[] int group, array[] int year,
                      array[] int market_year_idx, array[] int margin_obs,
                      vector a_market, vector gamma_structural_market,
                      vector b_market, vector b_quality_all,
                      real sigma_b_strat, real mu_b_fringe,
                      vector s0_prob,
                      vector cutoff_share,
                      vector sigma_share_market, vector sigma_margin_market,
                      real rho,
                      int supply_model, int market_sign, int use_cutoff,
                      int is_single_market) {

    real lp = 0;
    int N_slice = end - start + 1;
    real k = 0.2;

    for (n in 1:N_slice) {
      int i = slice_idx[n];
      int m = market[i];
      int mky = market_year_idx[i];
      real sigma_share = sigma_share_market[m];
      real sigma_margin = sigma_margin_market[m];
      real sigma_cond = sigma_margin * sqrt(1.0 - square(rho));
      real slope_cond = rho * (sigma_margin / sigma_share);

      real R_0 = s0_prob[mky];

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(revShare[i] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[i]] + (w * mu_b_fringe);
      // CES: negative sign on price (gamma > 1 requires negative coeff on log price)
      real delta = b_market[m] + b_quality +
        (market_sign * a_market[m] * log_price_std[i]);

      // revShare is normalized over inside goods. For structural outside
      // revenue share R_0, log(r_j | inside) = logit(R_0) + delta_j.
      real pred_log_r = logit(R_0) + delta;

      real gamma_structural = gamma_structural_market[m];

      // Unconditional revenue share: outside option disciplines the markup
      // (aggregative game). revShare is the inside/conditional share.
      real pred_lernerInv;
      if (supply_model == 1) {
        // Bertrand: lernerInv = gamma - (gamma-1)*r_unc
        pred_lernerInv = w * gamma_structural + (1-w) *
          (gamma_structural - (gamma_structural - 1.0) *
           (firm_revShare[i] * (1.0 - R_0)));
      } else if (supply_model == 2) {
        // Monopolistic competition: lernerInv = gamma
        pred_lernerInv = gamma_structural;
      } else if (supply_model == 3) {
        // Cournot: L = (1 + (gamma-1)*r_unc) / gamma => lernerInv = gamma / (1 + (gamma-1)*r_unc)
        real cournot_lernerInv = gamma_structural /
          (1.0 + (gamma_structural - 1.0) *
           (firm_revShare[i] * (1.0 - R_0)));
        pred_lernerInv = w * gamma_structural + (1-w) * cournot_lernerInv;
      } else {
        // 2nd Score Auction: p/mc = (1 - R_F)^(-1/(gamma-1))
        // lernerInv = 1 / (1 - (1 - r_unc)^(1/(gamma-1)))
        real r_unc_auc = fmin(0.99, firm_revShare[i] * (1.0 - R_0));
        real gamma_minus_one = fmax(gamma_structural - 1.0, 1e-6);
        real base_lernerInv = 1.0 / (1.0 - pow(1.0 - r_unc_auc, 1.0 / gamma_minus_one));
        pred_lernerInv = w * gamma_structural + (1-w) * base_lernerInv;
      }

      real mu2 = pred_lernerInv;
      real pred_log_lernerInv = log(mu2 + 1e-9);
      real err_share = log_revShare[i] - pred_log_r;
      real mu2_cond = pred_log_lernerInv + slope_cond * err_share;

      lp += normal_lpdf(log_revShare[i] | pred_log_r, sigma_share);
      if (margin_obs[i] == 1) {
        lp += normal_lpdf(log_lernerInv[i] | mu2_cond, sigma_cond);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  vector[N] revShare;
  vector[N] firm_revShare;
  vector[N] lernerInv;
  vector[N] log_revShare;
  vector[N] log_lernerInv;
  vector[N] log_price_std;
  real<lower=0> log_price_sd;
  int<lower=-1, upper=1> market_sign;
  int<lower=1> N_market;
  array[N] int<lower=1> market;
  int<lower=1> N_group;
  array[N] int<lower=1> group;
  int<lower=1> N_year;
  array[N] int<lower=1> year;
  int<lower=1, upper=4> supply_model;
  int<lower=0, upper=1> use_cutoff;
  int<lower=1> grainsize;
  int<lower=1> N_market_year;
  array[N] int<lower=1> market_year_idx;
  array[N] int<lower=0, upper=1> margin_obs;
  int<lower=0> K_s0;
  matrix[N_market_year, K_s0] X_s0;
  array[N_market_year] int<lower=1> mky_to_market;
  array[N_market_year] int<lower=1> mky_to_year;
  int<lower=0, upper=1> is_single_market;
  int<lower=0, upper=1> use_rho;
  int<lower=0, upper=1> use_market_resid_scale;
  int<lower=0, upper=1> use_alpha_market_pool;
  int<lower=0, upper=1> use_s0_market_year_pool;
  int<lower=0, upper=1> use_centered_alpha;
  int<lower=0, upper=1> use_centered_s0_pool;
  int<lower=0, upper=1> use_hmt;
  real avg_price_hmt;
  real avg_margin_hmt;
  real ssnip_hmt;
  real<lower=0> prior_sigma_share;
  real<lower=0> prior_sigma_margin;
  real<lower=0> prior_sigma_meanval_strat;
  real<lower=0> prior_sigma_meanval_fringe;
  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  real s0_prior_mean;
  real<lower=0> s0_prior_sd;
  real<lower=0> prior_sigma_alpha;
  real<lower=0> prior_sigma_s0_pool;
  real<lower=0> prior_sigma_beta_s0;
  real<lower=0> prior_sigma_resid_market;
  real<lower=0> prior_lkj;
  real<lower=0> prior_cutoff_alpha;
  real<lower=0> prior_cutoff_beta;
}

transformed data {
  array[N] int seq;
  int pool_alpha_by_market = (is_single_market == 0 && use_alpha_market_pool == 1) ? 1 : 0;
  array[N_market] int n_mky_per_market = rep_array(0, N_market);
  for (i in 1:N) seq[i] = i;
  for (mky in 1:N_market_year) {
    n_mky_per_market[mky_to_market[mky]] += 1;
  }
}

parameters {
  vector[is_single_market == 0 ? 1 : 0] mu_b_market_vec;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;
  vector[is_single_market == 0 ? N_market : 0] b_market_rot;
  vector[N_group] b_group_raw;
  real<lower=0> sigma_b_strat;
  vector[use_cutoff == 1 ? 1 : 0] mu_b_fringe_vec;
  real mu_log_a;
  vector<lower=0>[pool_alpha_by_market == 1 ? 1 : 0] sigma_log_a_vec;
  vector[(pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
  vector[(pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
  vector[use_s0_market_year_pool == 1 ? 0 : N_market_year] s0_logit;
  vector[use_s0_market_year_pool == 1 ? 1 : 0] mu_s0_pool_vec;
  vector<lower=0>[use_s0_market_year_pool == 1 ? 1 : 0] sigma_s0_pool_vec;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 0) ? N_market_year : 0] s0_pool_raw;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 1) ? N_market_year : 0] eta_s0_pool;
  vector<lower=0, upper=1>[(use_cutoff == 1) ? N_market : 0] cutoff_share;
  vector[use_rho == 1 ? 1 : 0] rho_raw;
  real log_sigma_share;
  real log_sigma_margin;
  vector[use_market_resid_scale == 1 ? N_market : 0] log_sigma_share_market_raw;
  vector[use_market_resid_scale == 1 ? N_market : 0] log_sigma_margin_market_raw;
  vector<lower=0>[use_market_resid_scale == 1 ? 1 : 0] sigma_log_sigma_share_market;
  vector<lower=0>[use_market_resid_scale == 1 ? 1 : 0] sigma_log_sigma_margin_market;
  vector[K_s0] beta_s0;
}

transformed parameters {
  real sigma_share = exp(log_sigma_share);
  real sigma_margin = exp(log_sigma_margin);
  vector[N_market] sigma_share_market;
  vector[N_market] sigma_margin_market;
  real rho = 0.0;
  real sigma_log_a = 0.0;
  vector[N_market] a_market;
  vector[N_market] b_market;
  vector[N_group] b_quality_all = sigma_b_strat * (b_group_raw - mean(b_group_raw));
  real mu_b_fringe = 0.0;
  vector[N_market_year] s0;
  vector[N_market_year] s0_prob;
  vector[N_market] mean_logit_s0_market = rep_vector(0.0, N_market);
  vector[N_market] gamma_structural_market;

  if (use_market_resid_scale == 1) {
    for (m in 1:N_market) {
      sigma_share_market[m] = exp(
        log_sigma_share +
        sigma_log_sigma_share_market[1] * log_sigma_share_market_raw[m]
      );
      sigma_margin_market[m] = exp(
        log_sigma_margin +
        sigma_log_sigma_margin_market[1] * log_sigma_margin_market_raw[m]
      );
    }
  } else {
    sigma_share_market = rep_vector(sigma_share, N_market);
    sigma_margin_market = rep_vector(sigma_margin, N_market);
  }

  if (use_rho == 1) rho = 0.95 * tanh(rho_raw[1]);

  if (pool_alpha_by_market == 1) {
    sigma_log_a = sigma_log_a_vec[1];
  }

  if (use_cutoff == 1) {
    mu_b_fringe = mu_b_fringe_vec[1];
  }

  for (m in 1:N_market_year) {
     real base_logit;
     real x_beta = (K_s0 > 0) ? (X_s0[m] * beta_s0) : 0.0;
     if (use_s0_market_year_pool == 1) {
       base_logit = (use_centered_s0_pool == 1)
         ? eta_s0_pool[m]
         : mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * s0_pool_raw[m];
     } else {
       base_logit = s0_logit[m];
     }
     s0[m] = base_logit + x_beta;
     s0_prob[m] = bounded_inv_logit(s0[m]);
  }

  // Rotate the market intercept by mean logit(s0). This removes arbitrary
  // outside-share level information from conditional inside revenue shares.
  for (mky in 1:N_market_year) {
    mean_logit_s0_market[mky_to_market[mky]] += logit(s0_prob[mky]);
  }
  for (mm in 1:N_market) {
    mean_logit_s0_market[mm] /= n_mky_per_market[mm];
  }

  if (is_single_market == 1) {
    // A single market has no cross-market intercept to speak of; normalize
    // it to zero and let s0's level be identified directly off the share
    // equation, rather than sampling a free intercept that would compete
    // with s0 for the same variation.
    b_market = rep_vector(0.0, N_market);
  } else {
    for (mm in 1:N_market) {
      if (n_mky_per_market[mm] > 1) {
        b_market[mm] = b_market_rot[mm] - mean_logit_s0_market[mm];
      } else {
        // One market-year: nothing to rotate against, so the intercept is
        // just b_market_rot directly and logit(R_0) survives uncancelled.
        b_market[mm] = b_market_rot[mm];
      }
    }
  }

  if (is_single_market == 1 || use_alpha_market_pool == 0) {
    a_market = rep_vector(exp(mu_log_a), N_market);
  } else if (use_centered_alpha == 1) {
    a_market = exp(log_a_market_pool);
  } else {
    a_market = exp(mu_log_a + sigma_log_a * r_market_a_raw);
  }

  // gamma = 1 + a_market / log_price_sd = sigma (elasticity of substitution)
  gamma_structural_market = 1.0 + a_market / log_price_sd;
}

model {
  mu_log_a ~ normal(prior_alpha_mean, prior_alpha_sd);
  if (pool_alpha_by_market == 1) {
    sigma_log_a_vec ~ normal(0, prior_sigma_alpha);
    if (use_centered_alpha == 1) {
      log_a_market_pool ~ normal(mu_log_a, sigma_log_a);
    } else {
      r_market_a_raw ~ std_normal();
    }
  }
  if (is_single_market == 0) {
    // No market intercept is sampled when N_market == 1 (b_market is fixed
    // at 0 in that case; see transformed parameters).
    mu_b_market_vec ~ normal(0, 1.0);
    sigma_b_market_vec ~ normal(0, 0.5);
    // Prior is on the derived (possibly rotated) b_market, not the raw
    // b_market_rot, so it targets the intended economic quantity even when
    // b_market = b_market_rot - mean_logit_s0_market(s0).
    b_market ~ normal(mu_b_market_vec[1], sigma_b_market_vec[1]);
  }
  b_group_raw ~ std_normal();

  // Informative prior pushes away from 0/infinity for better HMC geometry
  sigma_b_strat ~ normal(1.0, prior_sigma_meanval_strat);
  if (use_cutoff == 1) {
    cutoff_share ~ beta(prior_cutoff_alpha, prior_cutoff_beta);
    mu_b_fringe_vec ~ normal(0, prior_sigma_meanval_fringe);
  }

  if (use_s0_market_year_pool == 1) {
    mu_s0_pool_vec ~ normal(s0_prior_mean, s0_prior_sd);
    sigma_s0_pool_vec ~ normal(0, prior_sigma_s0_pool);
    if (use_centered_s0_pool == 1) {
      eta_s0_pool ~ normal(mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
    } else {
      s0_pool_raw ~ std_normal();
    }
  } else {
    s0_logit ~ normal(s0_prior_mean, s0_prior_sd);
  }
  beta_s0 ~ normal(0, prior_sigma_beta_s0);

  if (use_rho == 1) {
    // Two-dimensional LKJ prior, with the Jacobian for rho = 0.95*tanh(raw).
    target += (prior_lkj - 1.0) * log1m(square(rho)) +
      log1m(square(tanh(rho_raw[1])));
  }
  // Log-space priors prevent sigma collapse toward 0,
  // which causes extreme curvature and divergences via exp() share link
  log_sigma_share ~ normal(log(prior_sigma_share), 0.5);
  log_sigma_margin ~ normal(log(prior_sigma_margin), 1.0);
  if (use_market_resid_scale == 1) {
    log_sigma_share_market_raw ~ std_normal();
    log_sigma_margin_market_raw ~ std_normal();
    sigma_log_sigma_share_market ~ normal(0, prior_sigma_resid_market);
    sigma_log_sigma_margin_market ~ normal(0, prior_sigma_resid_market);
  }

  if (use_hmt == 1) {
     // For CES, aggregate elasticity uses gamma and revenue-weighted outside share
     real gamma_check = 1.0 + exp(mu_log_a) / log_price_sd;
     real R_0_sum = 0.0;
     real R_0_check;
     for (m in 1:N_market_year) {
       R_0_sum += bounded_inv_logit(s0[m]);
     }
     R_0_check = R_0_sum / N_market_year;
     real agg_elasticity = abs((gamma_check - 1.0) * R_0_check);
     real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
     target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  target += reduce_sum(partial_sum_fn, seq, grainsize,
                       lernerInv, revShare, firm_revShare, log_revShare, log_lernerInv,
                       log_price_std, market, group, year,
                       market_year_idx, margin_obs,
                       a_market, gamma_structural_market, b_market, b_quality_all, sigma_b_strat,
                       mu_b_fringe,
                       s0_prob, cutoff_share,
                       sigma_share_market, sigma_margin_market,
                       rho, supply_model, market_sign,
                       use_cutoff, is_single_market);
}

generated quantities {
  vector[N] pred_revShare;
  vector[N] pred_lernerInv;
  vector[N] log_lik;

  {
    real k = 0.2;

    for (n in 1:N) {
      int m = market[n];
      int mky = market_year_idx[n];
      real sigma_share_n = sigma_share_market[m];
      real sigma_margin_n = sigma_margin_market[m];
      real sigma_cond = sigma_margin_n * sqrt(1.0 - square(rho));
      real slope_cond = rho * (sigma_margin_n / sigma_share_n);

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(revShare[n] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[n]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality +
        (market_sign * a_market[m] * log_price_std[n]);

      real R_0 = s0_prob[mky];
      real pred_log_r = logit(R_0) + delta;
      pred_revShare[n] = exp(pred_log_r);

      real gamma_structural = gamma_structural_market[m];

      real r_unc = firm_revShare[n] * (1.0 - R_0);

      real pred_lInv;
      if (supply_model == 1) {
        // Bertrand: lernerInv = gamma - (gamma-1)*r_unc
        pred_lInv = w * gamma_structural + (1-w) * (gamma_structural - (gamma_structural - 1.0) * r_unc);
      } else if (supply_model == 2) {
        // Monopolistic competition: lernerInv = gamma
        pred_lInv = gamma_structural;
      } else if (supply_model == 3) {
        // Cournot: L = (1 + (gamma-1)*r_unc) / gamma => lernerInv = gamma / (1 + (gamma-1)*r_unc)
        real cournot_lInv = gamma_structural / (1.0 + (gamma_structural - 1.0) * r_unc);
        pred_lInv = w * gamma_structural + (1-w) * cournot_lInv;
      } else {
        // 2nd Score Auction: lernerInv = 1 / (1 - (1-r_unc)^(1/(gamma-1)))
        real r_unc_auc = fmin(0.99, r_unc);
        real gamma_minus_one = fmax(gamma_structural - 1.0, 1e-6);
        real base_lInv = 1.0 / (1.0 - pow(1.0 - r_unc_auc, 1.0 / gamma_minus_one));
        pred_lInv = w * gamma_structural + (1-w) * base_lInv;
      }

      pred_lernerInv[n] = pred_lInv;

      // Pointwise log-likelihood for LOO/WAIC. The margin term is gated by
      // margin_obs so masked (fringe/NA-margin) rows contribute only their
      // share likelihood, matching partial_sum_fn's training objective.
      real err_share = log_revShare[n] - pred_log_r;
      real mu2_cond = log(pred_lInv + 1e-9) + slope_cond * err_share;
      log_lik[n] = normal_lpdf(log_revShare[n] | pred_log_r, sigma_share_n);
      if (margin_obs[n] == 1) {
        log_lik[n] += normal_lpdf(log_lernerInv[n] | mu2_cond, sigma_cond);
      }
    }
  }
}
