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
                      vector lernerInv, vector revShare, vector log_revShare,
                      vector log_lernerInv, vector log_price_std,
                      array[] int market, array[] int group, array[] int year,
                      array[] int market_year_idx,
                      vector a_market, vector b_market, vector b_quality_all,
                      real sigma_b_strat, real mu_b_fringe,
                      vector s0,
                      vector cutoff_share, real sigma_share, real sigma_margin,
                      real rho,
                      real log_price_sd, int supply_model, int use_cutoff,
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

      real R_0 = bounded_inv_logit(s0[mky]);

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(revShare[i] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[i]] + (w * mu_b_fringe);
      // CES: negative sign on price (gamma > 1 requires negative coeff on log price)
      real delta = b_market[m] + b_quality - (a_market[m] * log_price_std[i]);

      real pred_log_r = log(R_0) + delta;

      // gamma = 1 + a_market / log_price_sd = sigma (elasticity of substitution)
      real gamma_structural = 1.0 + a_market[m] / log_price_sd;

      real pred_lernerInv;
      if (supply_model == 1) {
        // Bertrand: lernerInv = gamma - (gamma-1)*r_i
        pred_lernerInv = w * gamma_structural + (1-w) * (gamma_structural - (gamma_structural - 1.0) * revShare[i]);
      } else if (supply_model == 2) {
        // Monopolistic competition: lernerInv = gamma
        pred_lernerInv = gamma_structural;
      } else if (supply_model == 3) {
        // Cournot: L = (1 + (gamma-1)*r) / gamma => lernerInv = gamma / (1 + (gamma-1)*r)
        real cournot_lernerInv = gamma_structural / (1.0 + (gamma_structural - 1.0) * revShare[i]);
        pred_lernerInv = w * gamma_structural + (1-w) * cournot_lernerInv;
      } else {
        // 2nd Score Auction: p/mc = (1 - R_F)^(-1/(gamma-1))
        // lernerInv = 1 / (1 - (1 - r_unc)^(1/(gamma-1)))
        real r_unc = fmin(0.99, revShare[i] * (1.0 - R_0));
        real gamma_minus_one = fmax(gamma_structural - 1.0, 1e-6);
        real base_lernerInv = 1.0 / (1.0 - pow(1.0 - r_unc, 1.0 / gamma_minus_one));
        pred_lernerInv = w * gamma_structural + (1-w) * base_lernerInv;
      }

      real mu2 = pred_lernerInv;
      real pred_log_lernerInv = log(mu2 + 1e-9);
      real err_share = log_revShare[i] - pred_log_r;
      real mu2_cond = pred_log_lernerInv + slope_cond * err_share;

      lp += normal_lpdf(log_revShare[i] | pred_log_r, sigma_share);
      lp += normal_lpdf(log_lernerInv[i] | mu2_cond, sigma_cond);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  vector[N] revShare;
  vector[N] lernerInv;
  vector[N] log_revShare;
  vector[N] log_lernerInv;
  vector[N] log_price_std;
  real<lower=0> log_price_sd;
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
  real mu_log_a;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_log_a_vec;
  vector[is_single_market == 0 ? N_market : 0] r_market_a_raw;
  vector<lower=logit(0.01), upper=logit(0.99)>[N_market_year] s0_logit;
  vector<lower=0, upper=1>[(use_cutoff == 1) ? N_market : 0] cutoff_share;
  vector[use_rho == 1 ? 1 : 0] rho_raw;
  real log_sigma_share;
  real log_sigma_margin;
  vector<lower=0>[N_market_year > 1 ? 1 : 0] sigma_s0_vec;

  // Hierarchical s0 (outside revenue share for CES)
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
  if (use_rho == 1) rho = 0.95 * tanh(rho_raw[1]);

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
  if (use_cutoff == 1) {
    mu_b_fringe = mu_b_fringe_vec[1];
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
  mu_log_a ~ normal(prior_alpha_mean, prior_alpha_sd);
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
    rho_raw ~ normal(0, 0.75);
  }
  // Log-space priors prevent sigma collapse toward 0,
  // which causes extreme curvature and divergences via exp() share link
  log_sigma_share ~ normal(log(prior_sigma_share), 0.5);
  log_sigma_margin ~ normal(log(prior_sigma_margin), 1.0);

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
                       lernerInv, revShare, log_revShare, log_lernerInv,
                       log_price_std, market, group, year,
                       market_year_idx,
                       a_market, b_market, b_quality_all, sigma_b_strat,
                       mu_b_fringe,
                       s0, cutoff_share,
                       sigma_share, sigma_margin, rho, log_price_sd, supply_model,
                       use_cutoff, is_single_market);
}

generated quantities {
  vector[N] pred_revShare;
  vector[N] pred_lernerInv;
  vector[N] log_lik;
  vector[N_market_year] s0_prob;

  {
    real k = 0.2;
    real sigma_cond = sigma_margin * sqrt(1.0 - square(rho));
    real slope_cond = rho * (sigma_margin / sigma_share);

    for (m in 1:N_market_year) {
      s0_prob[m] = bounded_inv_logit(s0[m]);
    }

    for (n in 1:N) {
      int m = market[n];
      int mky = market_year_idx[n];

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(revShare[n] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[n]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality - (a_market[m] * log_price_std[n]);

      real R_0 = s0_prob[mky];
      real pred_log_r = log(R_0) + delta;
      pred_revShare[n] = exp(pred_log_r);

      real gamma_structural = 1.0 + a_market[m] / log_price_sd;

      real pred_lInv;
      if (supply_model == 1) {
        // Bertrand: lernerInv = gamma - (gamma-1)*r_i
        pred_lInv = w * gamma_structural + (1-w) * (gamma_structural - (gamma_structural - 1.0) * revShare[n]);
      } else if (supply_model == 2) {
        // Monopolistic competition: lernerInv = gamma
        pred_lInv = gamma_structural;
      } else if (supply_model == 3) {
        // Cournot: L = (1 + (gamma-1)*r) / gamma => lernerInv = gamma / (1 + (gamma-1)*r)
        real cournot_lInv = gamma_structural / (1.0 + (gamma_structural - 1.0) * revShare[n]);
        pred_lInv = w * gamma_structural + (1-w) * cournot_lInv;
      } else {
        // 2nd Score Auction: lernerInv = 1 / (1 - (1-r_unc)^(1/(gamma-1)))
        real r_unc = fmin(0.99, revShare[n] * (1.0 - R_0));
        real gamma_minus_one = fmax(gamma_structural - 1.0, 1e-6);
        real base_lInv = 1.0 / (1.0 - pow(1.0 - r_unc, 1.0 / gamma_minus_one));
        pred_lInv = w * gamma_structural + (1-w) * base_lInv;
      }

      pred_lernerInv[n] = pred_lInv;

      // Pointwise log-likelihood for LOO/WAIC
      real err_share = log_revShare[n] - pred_log_r;
      real mu2_cond = log(pred_lInv + 1e-9) + slope_cond * err_share;
      log_lik[n] = normal_lpdf(log_revShare[n] | pred_log_r, sigma_share)
                 + normal_lpdf(log_lernerInv[n] | mu2_cond, sigma_cond);
    }
  }
}
