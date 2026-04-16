// Bayesian Logit Demand + Supply Model for Merger Simulation
// Hierarchical model with market-level random effects
// Supports 4 supply models: Bertrand, Monopolistic Competition, Cournot, 2nd Score Auction

functions {
  real soft_hmt_penalty(real agg_elast, real max_elast, real scale) {
    if (agg_elast > max_elast) {
      return normal_lpdf(agg_elast | max_elast, scale) - normal_lpdf(max_elast | max_elast, scale);
    }
    return 0.0;
  }

  real partial_sum_fn(array[] int slice_idx, int start, int end,
                      vector marginInv, vector shareIn, vector price_std,
                      array[] int market, array[] int group, array[] int year,
                      array[] int market_year_idx,
                      vector a_market, vector b_market, vector b_quality_all,
                      real sigma_b_strat, real mu_b_fringe, real sigma_b_fringe,
                      vector s0,
                      vector cutoff_share, real sigma_share, real sigma_margin,
                      real rho,
                      real price_sd, int supply_model, int use_cutoff,
                      int is_single_market) {

    real lp = 0;
    int N_slice = end - start + 1;
    real sigma_cond = sigma_margin * sqrt(1.0 - square(rho));
    real slope_cond = rho * (sigma_margin / sigma_share);
    real k = 0.2;

    for (n in 1:N_slice) {
      int i = slice_idx[n];
      int m = market[i];
      int mky = market_year_idx[i];

      real s_out = inv_logit(s0[mky]);

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(shareIn[i] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[i]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality + (a_market[m] * price_std[i]);

      real mu1 = s0[mky] + delta;
      real pred_s = exp(mu1);

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
      real err_share = shareIn[i] - pred_s;
      real mu2_cond = mu2 + slope_cond * err_share;

      lp += normal_lpdf(shareIn[i] | pred_s, sigma_share);
      lp += normal_lpdf(marginInv[i] | mu2_cond, sigma_cond);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  vector[N] shareIn;
  vector[N] marginInv;
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
  real<lower=0> prior_sigma_beta_s0;
  real<lower=0> prior_lkj;
  real<lower=0> prior_cutoff_alpha;
  real<lower=0> prior_cutoff_beta;
  vector<lower=0, upper=1>[N_market_year] implied_s0;
}

transformed data {
  array[N] int seq;
  for (i in 1:N) seq[i] = i;
}

parameters {
  vector[is_single_market == 0 ? 1 : 0] mu_b_market_vec;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;
  vector[is_single_market == 0 ? N_market : 0] b_market_raw;
  vector[N_group] b_group_raw;
  real<lower=0> sigma_b_strat;
  vector[use_cutoff == 1 ? 1 : 0] mu_b_fringe_vec;
  vector<lower=0>[use_cutoff == 1 ? 1 : 0] sigma_b_fringe_vec;
  real mu_log_a;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_log_a_vec;
  vector[is_single_market == 0 ? N_market : 0] r_market_a_raw;
  vector<lower=logit(0.01), upper=logit(0.99)>[N_market_year] s0_logit;
  vector<lower=0, upper=1>[(use_cutoff == 1) ? N_market : 0] cutoff_share;
  vector<lower=-1, upper=1>[use_rho == 1 ? 1 : 0] rho_param;
  real log_sigma_share;
  real log_sigma_margin;
  vector<lower=0>[N_market_year > 1 ? 1 : 0] sigma_s0_vec;

  // Hierarchical s0
  vector[is_single_market == 0 ? N_market : 0] mu_s0_market_raw;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_s0_market;
  vector[N_year > 1 ? N_year : 0] year_effect_s0_raw;
  vector<lower=0>[N_year > 1 ? 1 : 0] sigma_year_s0;

  vector[K_s0] beta_s0;
}

transformed parameters {
  real sigma_share = exp(log_sigma_share);
  real sigma_margin = exp(log_sigma_margin);

  real rho = 0.0;
  if (use_rho == 1) rho = rho_param[1];

  real sigma_log_a = 0.0;
  if (is_single_market == 0) sigma_log_a = sigma_log_a_vec[1];

  vector[N_market] a_market;
  if (is_single_market == 1) {
    a_market = rep_vector(exp(mu_log_a), N_market);
  } else {
    a_market = exp(mu_log_a + sigma_log_a * r_market_a_raw);
  }

  vector[N_market] b_market;
  if (is_single_market == 1) {
    b_market = rep_vector(0.0, N_market);
  } else {
    b_market = mu_b_market_vec[1] + sigma_b_market_vec[1] * b_market_raw;
  }
  vector[N_group] b_quality_all = sigma_b_strat * b_group_raw;

  real mu_b_fringe = 0.0;
  real sigma_b_fringe = 0.0;
  if (use_cutoff == 1) {
    mu_b_fringe = mu_b_fringe_vec[1];
    sigma_b_fringe = sigma_b_fringe_vec[1];
  }

  vector[is_single_market == 0 ? N_market : 0] mu_s0_market;
  if (is_single_market == 1) {
    mu_s0_market = rep_vector(0.0, 0);
  } else {
    mu_s0_market = sigma_s0_market[1] * mu_s0_market_raw;
  }

  vector[N_year > 1 ? N_year : 0] year_effect_s0;
  if (N_year > 1) {
    year_effect_s0 = sigma_year_s0[1] * year_effect_s0_raw;
  } else {
    year_effect_s0 = rep_vector(0.0, 0);
  }

  vector[N_market_year] s0;
  for (m in 1:N_market_year) {
     int e = mky_to_market[m];
     int t = mky_to_year[m];
     real base_logit = s0_logit[m];
     real x_beta = (K_s0 > 0) ? (X_s0[m] * beta_s0) : 0.0;
     real mkt_effect = 0.0;
     if (is_single_market == 0) mkt_effect = mu_s0_market[e];

     real yr_effect = 0.0;
     if (N_year > 1) yr_effect = year_effect_s0[t];

     s0[m] = base_logit + x_beta + mkt_effect + yr_effect;
  }
}

model {
  mu_log_a ~ normal(log(prior_alpha_mean), prior_alpha_sd);
  if (is_single_market == 1) {
    // No prior needed for sigma_log_a
  } else {
    sigma_log_a_vec ~ normal(0, prior_sigma_alpha);
    r_market_a_raw ~ std_normal();
    mu_b_market_vec ~ normal(0, 1.0);
    b_market_raw ~ std_normal();
    sigma_b_market_vec ~ normal(0, 0.5);
  }
  b_group_raw ~ std_normal();

  // Informative prior pushes away from 0/infinity for better HMC geometry
  sigma_b_strat ~ normal(1.0, 0.5);
  if (use_cutoff == 1) {
    cutoff_share ~ beta(prior_cutoff_alpha, prior_cutoff_beta);
    sigma_b_fringe_vec ~ normal(0, prior_sigma_meanval_fringe);
    mu_b_fringe_vec ~ normal(0, 2);
  }

  if (N_market_year > 1) {
    s0_logit ~ normal(logit(implied_s0), sigma_s0_vec[1]);
    sigma_s0_vec ~ normal(0, 1.0);
  } else {
    // If substantial outside share (> 2%), center prior on data.
    // Otherwise fall back to user-specified informative prior.
    if (implied_s0[1] > 0.02) {
      s0_logit ~ normal(logit(implied_s0), 1.0);
    } else {
      s0_logit ~ normal(s0_prior_mean, s0_prior_sd);
    }
  }
  beta_s0 ~ normal(0, prior_sigma_beta_s0);

  if (is_single_market == 0) {
    mu_s0_market_raw ~ std_normal();
    sigma_s0_market ~ normal(0, 1.0);
  }
  if (N_year > 1) {
    year_effect_s0_raw ~ std_normal();
    sigma_year_s0 ~ normal(0, 1.0);
  }
  if (use_rho == 1) {
    target += (prior_lkj - 1) * log1m(square(rho_param[1]));
  }
  // Log-space priors prevent sigma collapse toward 0,
  // which causes extreme curvature and divergences via exp() share link
  log_sigma_share ~ normal(log(prior_sigma_share), 0.5);
  log_sigma_margin ~ normal(log(prior_sigma_margin), 1.0);

  if (use_hmt == 1) {
     real alpha_check = exp(mu_log_a) / price_sd;
     real s0_check = mean(inv_logit(s0));
     real agg_elasticity = abs(alpha_check * avg_price_hmt * s0_check);
     real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
     target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  target += reduce_sum(partial_sum_fn, seq, grainsize,
                       marginInv, shareIn, price_std, market, group, year,
                       market_year_idx,
                       a_market, b_market, b_quality_all, sigma_b_strat,
                       mu_b_fringe, sigma_b_fringe,
                       s0, cutoff_share,
                       sigma_share, sigma_margin, rho, price_sd, supply_model,
                       use_cutoff, is_single_market);
}

generated quantities {
  vector[N] pred_shareIn;
  vector[N] pred_marginInv;
  vector[N] log_lik;

  {
    real k = 0.2;
    real sigma_cond = sigma_margin * sqrt(1.0 - square(rho));
    real slope_cond = rho * (sigma_margin / sigma_share);

    for (n in 1:N) {
      int m = market[n];
      int mky = market_year_idx[n];

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(shareIn[n] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[n]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality + (a_market[m] * price_std[n]);

      real mu1 = s0[mky] + delta;
      pred_shareIn[n] = exp(mu1);

      real s_out = inv_logit(s0[mky]);
      if (s_out < 0.01) s_out = 0.01;
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
      real err_share = shareIn[n] - pred_shareIn[n];
      real mu2_cond = inv_markup + slope_cond * err_share;
      log_lik[n] = normal_lpdf(shareIn[n] | pred_shareIn[n], sigma_share)
                 + normal_lpdf(marginInv[n] | mu2_cond, sigma_cond);
    }
  }
}
