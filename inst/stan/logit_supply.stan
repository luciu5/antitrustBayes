// Bayesian Logit Demand + Supply Model for Merger Simulation
// Hierarchical model with market-level random effects
// Supports 4 supply models: Bertrand, Monopolistic Competition, Cournot, 2nd Score Auction

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
                      vector marginInv, vector shareIn, vector log_shareIn,
                      vector log_marginInv, vector price_std,
                      array[] int market, array[] int group, array[] int year,
                      array[] int market_year_idx,
                      vector a_market, vector b_market, vector b_quality_all,
                      real sigma_b_strat, real mu_b_fringe,
                      vector s0,
                      vector cutoff_share,
                      vector sigma_share_market, vector sigma_margin_market,
                      real rho,
                      real price_sd, int supply_model, int use_cutoff,
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

      real s_out = bounded_inv_logit(s0[mky]);

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(shareIn[i] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[i]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality + (a_market[m] * price_std[i]);

      real pred_log_s = log(s_out) + delta;

      real alpha_structural = a_market[m] / price_sd;

      real inv_markup;
      if (supply_model == 1) {
        // Bertrand: margin = 1 / (alpha * (1 - s_uncond))
        inv_markup = alpha_structural * (w + (1-w) * (1.0 - shareIn[i] * (1.0 - s_out)));
      } else if (supply_model == 2) {
        // Monopolistic competition: alpha = 1/(p-c)
        inv_markup = alpha_structural;
      } else if (supply_model == 3) {
        // Cournot: margin = (s0 + s) / (alpha * s0)
        inv_markup = alpha_structural * (w + (1-w) * (s_out / (s_out + shareIn[i])));
      } else {
        // 2nd Score Auction: margin = -log(1 - s) / (alpha * s)
        real denom = fmin(0.99, shareIn[i] * (1.0 - s_out));
        inv_markup = alpha_structural * (w + (1-w) * (denom / -log(1.0 - denom)));
      }

      real mu2 = inv_markup;
      real pred_log_marginInv = log(mu2 + 1e-9);
      real err_share = log_shareIn[i] - pred_log_s;
      real mu2_cond = pred_log_marginInv + slope_cond * err_share;

      lp += normal_lpdf(log_shareIn[i] | pred_log_s, sigma_share);
      lp += normal_lpdf(log_marginInv[i] | mu2_cond, sigma_cond);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  vector[N] shareIn;
  vector[N] marginInv;
  vector[N] log_shareIn;
  vector[N] log_marginInv;
  vector[N] price_std;
  real<lower=0> price_sd;
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
  int<lower=0> K_s0;
  matrix[N_market_year, K_s0] X_s0;
  array[N_market_year] int<lower=1> mky_to_market;
  array[N_market_year] int<lower=1> mky_to_year;
  int<lower=0, upper=1> is_single_market;
  int<lower=0, upper=1> use_rho;
  int<lower=0, upper=1> use_market_resid_scale;
  int<lower=0, upper=1> use_alpha_market_pool;
  int<lower=0, upper=1> use_s0_market_year_pool;
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
  vector<lower=0, upper=1>[N_market_year] implied_s0;
}

transformed data {
  array[N] int seq;
  vector[N_market_year] implied_s0_logit;
  real implied_s0_logit_mean = 0.0;
  for (i in 1:N) seq[i] = i;
  for (m in 1:N_market_year) {
    implied_s0_logit[m] = logit(implied_s0[m]);
    implied_s0_logit_mean += implied_s0_logit[m] / N_market_year;
  }
}

parameters {
  vector[is_single_market == 0 ? 1 : 0] mu_b_market_vec;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;
  vector[is_single_market == 0 ? N_market : 0] b_market_raw;
  vector[N_group] b_group_raw;
  real<lower=0> sigma_b_strat;
  vector[use_cutoff == 1 ? 1 : 0] mu_b_fringe_vec;
  real mu_log_a;
  vector<lower=0>[(is_single_market == 0 && use_alpha_market_pool == 1) ? 1 : 0] sigma_log_a_vec;
  vector[(is_single_market == 0 && use_alpha_market_pool == 1) ? N_market : 0] r_market_a_raw;
  vector<lower=logit(0.01), upper=logit(0.99)>[use_s0_market_year_pool == 1 ? 0 : N_market_year] s0_logit;
  vector[use_s0_market_year_pool == 1 ? 1 : 0] mu_s0_pool_vec;
  vector<lower=0>[use_s0_market_year_pool == 1 ? 1 : 0] sigma_s0_pool_vec;
  vector[use_s0_market_year_pool == 1 ? N_market_year : 0] s0_pool_raw;
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
  vector[N_group] b_quality_all = sigma_b_strat * b_group_raw;
  real mu_b_fringe = 0.0;
  vector[N_market_year] s0;

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

  if (is_single_market == 0 && use_alpha_market_pool == 1) {
    sigma_log_a = sigma_log_a_vec[1];
  }

  if (is_single_market == 1 || use_alpha_market_pool == 0) {
    a_market = rep_vector(exp(mu_log_a), N_market);
  } else {
    a_market = exp(mu_log_a + sigma_log_a * r_market_a_raw);
  }

  if (is_single_market == 1) {
    b_market = rep_vector(0.0, N_market);
  } else {
    b_market = mu_b_market_vec[1] + sigma_b_market_vec[1] * b_market_raw;
  }

  if (use_cutoff == 1) {
    mu_b_fringe = mu_b_fringe_vec[1];
  }

  for (m in 1:N_market_year) {
     real base_logit;
     real x_beta = (K_s0 > 0) ? (X_s0[m] * beta_s0) : 0.0;
     if (use_s0_market_year_pool == 1) {
       base_logit = mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * s0_pool_raw[m];
     } else {
       base_logit = s0_logit[m];
     }
     s0[m] = base_logit + x_beta;
  }
}

model {
  mu_log_a ~ normal(prior_alpha_mean, prior_alpha_sd);
  if (is_single_market == 0 && use_alpha_market_pool == 1) {
    sigma_log_a_vec ~ normal(0, prior_sigma_alpha);
    r_market_a_raw ~ std_normal();
  }
  if (is_single_market == 0) {
    mu_b_market_vec ~ normal(0, 1.0);
    b_market_raw ~ std_normal();
    sigma_b_market_vec ~ normal(0, 0.5);
  }
  b_group_raw ~ std_normal();

  // Informative prior pushes away from 0/infinity for better HMC geometry
  sigma_b_strat ~ normal(1.0, 0.5);
  if (use_cutoff == 1) {
    cutoff_share ~ beta(prior_cutoff_alpha, prior_cutoff_beta);
    mu_b_fringe_vec ~ normal(0, 2);
  }

  if (use_s0_market_year_pool == 1) {
    mu_s0_pool_vec ~ normal(implied_s0_logit_mean, s0_prior_sd);
    sigma_s0_pool_vec ~ normal(0, prior_sigma_s0_pool);
    s0_pool_raw ~ std_normal();
  } else {
    if (N_market_year > 1) {
      s0_logit ~ normal(implied_s0_logit, s0_prior_sd);
    } else {
      // If substantial outside share (> 2%), center prior on data.
      // Otherwise fall back to user-specified informative prior.
      if (implied_s0[1] > 0.02) {
        s0_logit ~ normal(implied_s0_logit, 1.0);
      } else {
        s0_logit ~ normal(s0_prior_mean, s0_prior_sd);
      }
    }
  }
  beta_s0 ~ normal(0, prior_sigma_beta_s0);

  if (use_rho == 1) {
    rho_raw ~ normal(0, 0.75);
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
     real alpha_check = exp(mu_log_a) / price_sd;
     real s0_sum = 0.0;
     real s0_check;
     for (m in 1:N_market_year) {
       s0_sum += bounded_inv_logit(s0[m]);
     }
     s0_check = s0_sum / N_market_year;
     real agg_elasticity = abs(alpha_check * avg_price_hmt * s0_check);
     real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
     target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  target += reduce_sum(partial_sum_fn, seq, grainsize,
                       marginInv, shareIn, log_shareIn, log_marginInv,
                       price_std, market, group, year,
                       market_year_idx,
                       a_market, b_market, b_quality_all, sigma_b_strat,
                       mu_b_fringe,
                       s0, cutoff_share,
                       sigma_share_market, sigma_margin_market,
                       rho, price_sd, supply_model,
                       use_cutoff, is_single_market);
}

generated quantities {
  vector[N] pred_shareIn;
  vector[N] pred_marginInv;
  vector[N] log_lik;
  vector[N_market_year] s0_prob;

  {
    real k = 0.2;
    for (m in 1:N_market_year) {
      s0_prob[m] = bounded_inv_logit(s0[m]);
    }

    for (n in 1:N) {
      int m = market[n];
      int mky = market_year_idx[n];
      real sigma_share_n = sigma_share_market[m];
      real sigma_margin_n = sigma_margin_market[m];
      real sigma_cond = sigma_margin_n * sqrt(1.0 - square(rho));
      real slope_cond = rho * (sigma_margin_n / sigma_share_n);

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(shareIn[n] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[n]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality + (a_market[m] * price_std[n]);

      real s_out = s0_prob[mky];
      real pred_log_s = log(s_out) + delta;
      pred_shareIn[n] = exp(pred_log_s);

      real alpha_structural = a_market[m] / price_sd;

      real inv_markup;
      if (supply_model == 1) {
        // Bertrand
        inv_markup = alpha_structural * (w + (1-w) * (1.0 - shareIn[n] * (1.0 - s_out)));
      } else if (supply_model == 2) {
        // Monopolistic competition: alpha = 1/(p-c)
        inv_markup = alpha_structural;
      } else if (supply_model == 3) {
        // Cournot
        inv_markup = alpha_structural * (w + (1-w) * (s_out / (s_out + shareIn[n])));
      } else {
        // 2nd Score Auction
        real denom = fmin(0.99, shareIn[n] * (1.0 - s_out));
        inv_markup = alpha_structural * (w + (1-w) * (denom / -log(1.0 - denom)));
      }

      pred_marginInv[n] = inv_markup;

      // Pointwise log-likelihood for LOO/WAIC
      real err_share = log_shareIn[n] - pred_log_s;
      real mu2_cond = log(inv_markup + 1e-9) + slope_cond * err_share;
      log_lik[n] = normal_lpdf(log_shareIn[n] | pred_log_s, sigma_share_n)
                 + normal_lpdf(log_marginInv[n] | mu2_cond, sigma_cond);
    }
  }
}
