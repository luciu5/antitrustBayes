// Bayesian CES demand and supply model with the same observed-outcome contract
// as logit_supply.stan. Conditional quantity shares and inverse normalized
// margins are coupled per observation through a free cross-equation
// correlation rho, using a bivariate Student-t factored as
// p(share) * p(margin | share). The latent outside share R0/s0 is identified
// jointly through the share equation (log conditional share ~
// logit(s0) + delta) and the equilibrium margin condition.
// likelihood_mode == 2 (share-only) omits the margin equation entirely, so
// it lacks the equilibrium-margin information that principally identifies
// s0's level; s0 may be weakly identified or prior-dominated in that mode,
// depending on market/time/intercept structure. This is an expected
// property of the model and is not compensated for with extra pooling.

functions {
  real bounded_inv_logit(real eta) {
    real eps = 1e-6;
    return eps + (1.0 - 2.0 * eps) * inv_logit(eta);
  }

  real soft_hmt_penalty(real agg_elast, real max_elast, real scale) {
    if (agg_elast > max_elast) {
      return normal_lpdf(agg_elast | max_elast, scale)
        - normal_lpdf(max_elast | max_elast, scale);
    }
    return 0.0;
  }

  real conduct_weight_ces(real firm_share_inside, real cutoff, int use_cutoff) {
    real x;
    if (use_cutoff == 0) return 0.0;
    x = -(firm_share_inside - cutoff) / 0.2;
    if (x > 10) return 1.0;
    if (x < -10) return 0.0;
    return inv_logit(x);
  }

  real ces_inv_margin(real firm_rev_share_inside, real R0, real gamma,
                      real fringe_weight, int supply_model,
                      real bargaining_power) {
    real r = fmin(0.999999, firm_rev_share_inside * (1.0 - R0));
    real conduct;
    if (supply_model == 1) {
      conduct = gamma - (gamma - 1.0) * r;
    } else if (supply_model == 2) {
      conduct = gamma;
    } else if (supply_model == 3) {
      real lerner = 1.0 / gamma
        + (gamma - 1.0) * r
          / (gamma * (1.0 + (gamma - 1.0) * R0));
      conduct = 1.0 / lerner;
    } else if (supply_model == 4) {
      conduct = 1.0 / (1.0 - pow(1.0 - r, 1.0 / (gamma - 1.0)));
    } else {
      conduct = (gamma - (gamma - 1.0) * r) / (1.0 - bargaining_power);
    }
    return fringe_weight * gamma + (1.0 - fringe_weight) * conduct;
  }
}

data {
  int<lower=1> N;
  vector[N] shareIn;
  vector[N] revShare;
  vector[N] firm_revShare;
  vector[N] clr_shareIn;
  vector[N] log_shareIn;
  vector[N] log_marginNormInv;
  vector[N] price;
  vector[N] log_price_std;
  real<lower=0> log_price_sd;
  int<lower=-1, upper=1> market_sign;
  int<lower=1> N_market;
  array[N] int<lower=1, upper=N_market> market;
  int<lower=1> N_group;
  array[N] int<lower=1, upper=N_group> group;
  int<lower=1> N_year;
  array[N] int<lower=1, upper=N_year> year;
  int<lower=1, upper=5> supply_model;
  int<lower=0, upper=1> use_cutoff;
  int<lower=1> N_market_year;
  vector<lower=0>[N_market_year] outside_price;
  vector<lower=0>[N_market_year] inside_price_index;
  array[N] int<lower=1, upper=N_market_year> market_year_idx;
  array[N_market_year] int<lower=1> mky_n;
  array[N_market_year] int<lower=1> mky_start;
  array[N] int<lower=1, upper=N> obs_order;
  int<lower=1> N_cluster;
  array[N] int<lower=1, upper=N_cluster> cluster;
  int<lower=1, upper=3> likelihood_mode;
  array[N] int<lower=0, upper=1> margin_obs;
  int<lower=0> K_s0;
  matrix[N_market_year, K_s0] X_s0;
  array[N_market_year] int<lower=1, upper=N_market> mky_to_market;
  array[N_market_year] int<lower=1, upper=N_year> mky_to_year;
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
  real<lower=0> prior_log_sigma_share_sd;
  real<lower=0> prior_sigma_margin;
  real<lower=0> prior_log_sigma_margin_sd;
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
  real<lower=0> prior_cutoff_alpha;
  real<lower=0> prior_cutoff_beta;
  int<lower=0, upper=1> estimate_df;
  real<lower=2> resid_df_fixed;
  real<lower=0> prior_nu_shape;
  real<lower=0> prior_nu_rate;
  real<lower=0> prior_lkj;
  real<lower=0, upper=1> bargaining_prior_lower;
  real<lower=0, upper=1> bargaining_prior_upper;
  real<lower=0> bargaining_prior_alpha;
  real<lower=0> bargaining_prior_beta;
}

transformed data {
  int is_bargaining = supply_model == 5;
  int pool_alpha_by_market = is_single_market == 0 && use_alpha_market_pool == 1;
  array[N_market] int n_mky_per_market = rep_array(0, N_market);
  int any_paired = 0;
  int rho_active;
  for (mky in 1:N_market_year) {
    n_mky_per_market[mky_to_market[mky]] += 1;
  }
  for (i in 1:N) {
    if (margin_obs[i] == 1) any_paired = 1;
  }
  rho_active = (use_rho == 1 && likelihood_mode == 1 && any_paired == 1) ? 1 : 0;
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
  vector<lower=0, upper=1>[use_cutoff == 1 ? N_market : 0] cutoff_share;
  vector[rho_active == 1 ? 1 : 0] rho_raw;
  real log_sigma_share;
  real log_sigma_margin;
  vector[use_market_resid_scale == 1 ? N_market : 0] log_sigma_share_market_raw;
  vector[use_market_resid_scale == 1 ? N_market : 0] log_sigma_margin_market_raw;
  vector<lower=0>[use_market_resid_scale == 1 ? 1 : 0] sigma_log_sigma_share_market;
  vector<lower=0>[use_market_resid_scale == 1 ? 1 : 0] sigma_log_sigma_margin_market;
  vector[K_s0] beta_s0;
  vector<lower=0, upper=1>[is_bargaining] bargaining_unit;
  vector<lower=0>[estimate_df == 1 ? 1 : 0] nu_minus_2;
}

transformed parameters {
  real sigma_share = exp(log_sigma_share);
  real sigma_margin = exp(log_sigma_margin);
  vector[N_market] sigma_share_market;
  vector[N_market] sigma_margin_market;
  real sigma_log_a = 0.0;
  vector[N_market] a_market;
  vector[N_market] b_market;
  vector[N_group] b_quality_all = sigma_b_strat * (b_group_raw - mean(b_group_raw));
  real mu_b_fringe = 0.0;
  vector[N_market_year] s0;
  vector[N_market_year] s0_prob;
  vector[N_market_year] R0_prob;
  vector[N_market] mean_logit_s0_market = rep_vector(0.0, N_market);
  vector[N_market] gamma_structural_market;
  real bargaining_power = 0.0;
  real rho = 0.0;
  real nu = estimate_df == 1 ? (2.0 + nu_minus_2[1]) : resid_df_fixed;

  if (use_market_resid_scale == 1) {
    for (m in 1:N_market) {
      sigma_share_market[m] = exp(log_sigma_share
        + sigma_log_sigma_share_market[1] * log_sigma_share_market_raw[m]);
      sigma_margin_market[m] = exp(log_sigma_margin
        + sigma_log_sigma_margin_market[1] * log_sigma_margin_market_raw[m]);
    }
  } else {
    sigma_share_market = rep_vector(sigma_share, N_market);
    sigma_margin_market = rep_vector(sigma_margin, N_market);
  }
  if (pool_alpha_by_market == 1) sigma_log_a = sigma_log_a_vec[1];
  if (use_cutoff == 1) mu_b_fringe = mu_b_fringe_vec[1];
  if (is_bargaining == 1) {
    bargaining_power = bargaining_prior_lower
      + (bargaining_prior_upper - bargaining_prior_lower) * bargaining_unit[1];
  }
  if (rho_active == 1) rho = 0.95 * tanh(rho_raw[1]);

  for (mky in 1:N_market_year) {
    real base_logit;
    real x_beta = K_s0 > 0 ? X_s0[mky] * beta_s0 : 0.0;
    if (use_s0_market_year_pool == 1) {
      base_logit = use_centered_s0_pool == 1
        ? eta_s0_pool[mky]
        : mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * s0_pool_raw[mky];
    } else base_logit = s0_logit[mky];
    s0[mky] = base_logit + x_beta;
    s0_prob[mky] = bounded_inv_logit(s0[mky]);
    R0_prob[mky] = bounded_inv_logit(
      logit(s0_prob[mky]) + log(outside_price[mky])
      - log(inside_price_index[mky])
    );
    mean_logit_s0_market[mky_to_market[mky]] += logit(s0_prob[mky]);
  }
  for (m in 1:N_market) {
    mean_logit_s0_market[m] /= n_mky_per_market[m];
  }

  if (is_single_market == 1) {
    b_market = rep_vector(0.0, N_market);
  } else {
    for (m in 1:N_market) {
      if (n_mky_per_market[m] > 1) {
        b_market[m] = b_market_rot[m] - mean_logit_s0_market[m];
      } else {
        b_market[m] = b_market_rot[m];
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
  gamma_structural_market = 1.0 + a_market / log_price_sd;
}

model {
  mu_log_a ~ normal(prior_alpha_mean, prior_alpha_sd);
  if (pool_alpha_by_market == 1) {
    sigma_log_a_vec ~ normal(0, prior_sigma_alpha);
    if (use_centered_alpha == 1) log_a_market_pool ~ normal(mu_log_a, sigma_log_a);
    else r_market_a_raw ~ std_normal();
  }
  if (is_single_market == 0) {
    mu_b_market_vec ~ normal(0, 1);
    sigma_b_market_vec ~ normal(0, 0.5);
    b_market ~ normal(mu_b_market_vec[1], sigma_b_market_vec[1]);
  }
  b_group_raw ~ std_normal();
  sigma_b_strat ~ normal(1, prior_sigma_meanval_strat);
  if (use_cutoff == 1) {
    cutoff_share ~ beta(prior_cutoff_alpha, prior_cutoff_beta);
    mu_b_fringe_vec ~ normal(0, prior_sigma_meanval_fringe);
  }
  if (use_s0_market_year_pool == 1) {
    mu_s0_pool_vec ~ normal(s0_prior_mean, s0_prior_sd);
    sigma_s0_pool_vec ~ normal(0, prior_sigma_s0_pool);
    if (use_centered_s0_pool == 1) {
      eta_s0_pool ~ normal(mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
    } else s0_pool_raw ~ std_normal();
  } else s0_logit ~ normal(s0_prior_mean, s0_prior_sd);
  beta_s0 ~ normal(0, prior_sigma_beta_s0);
  if (is_bargaining == 1) {
    bargaining_unit ~ beta(bargaining_prior_alpha, bargaining_prior_beta);
  }
  log_sigma_share ~ normal(log(prior_sigma_share), prior_log_sigma_share_sd);
  log_sigma_margin ~ normal(log(prior_sigma_margin), prior_log_sigma_margin_sd);
  if (use_market_resid_scale == 1) {
    log_sigma_share_market_raw ~ std_normal();
    log_sigma_margin_market_raw ~ std_normal();
    sigma_log_sigma_share_market ~ normal(0, prior_sigma_resid_market);
    sigma_log_sigma_margin_market ~ normal(0, prior_sigma_resid_market);
  }
  if (estimate_df == 1) nu_minus_2 ~ gamma(prior_nu_shape, prior_nu_rate);
  if (rho_active == 1) {
    target += (prior_lkj - 1.0) * log1m(square(rho))
      + log1m(square(tanh(rho_raw[1])));
  }

  if (use_hmt == 1 && likelihood_mode != 2) {
    real gamma_check = 1.0 + exp(mu_log_a) / log_price_sd;
    real R0_check = mean(R0_prob);
    real agg_elasticity = 1.0 + (gamma_check - 1.0) * R0_check;
    real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
    target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  for (i in 1:N) {
    int m = market[i];
    int mky = market_year_idx[i];
    real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
    real w = conduct_weight_ces(firm_revShare[i], cutoff, use_cutoff);
    real b_quality = b_quality_all[group[i]] + w * mu_b_fringe;
    real delta = b_market[m] + b_quality
      + market_sign * (a_market[m] + log_price_sd) * log_price_std[i]
      + gamma_structural_market[m] * log(outside_price[mky]);
    real pred_log_q = logit(s0_prob[mky]) + delta;
    real sigma_share_eff = sigma_share_market[m];
    real sigma_margin_eff = sigma_margin_market[m];

    if (likelihood_mode != 3) {
      real scale_share_marginal = sigma_share_eff * sqrt((nu - 2) / nu);
      target += student_t_lpdf(log_shareIn[i] | nu, pred_log_q, scale_share_marginal);
    }

    if (likelihood_mode != 2 && margin_obs[i] == 1) {
      real pred_inv_margin = ces_inv_margin(
        firm_revShare[i], R0_prob[mky], gamma_structural_market[m], w,
        supply_model, bargaining_power
      );
      real pred_log_margin = log(pred_inv_margin + 1e-9);
      if (rho_active == 1) {
        real err_share = log_shareIn[i] - pred_log_q;
        real z_s = err_share / sigma_share_eff;
        real slope_cond = rho * (sigma_margin_eff / sigma_share_eff);
        real mu2_cond = pred_log_margin + slope_cond * err_share;
        real scale2_cond = sigma_margin_eff * sqrt(1 - square(rho))
          * sqrt((nu - 2 + square(z_s)) / (nu + 1));
        target += student_t_lpdf(log_marginNormInv[i] | nu + 1, mu2_cond, scale2_cond);
      } else {
        real scale_margin_marginal = sigma_margin_eff * sqrt((nu - 2) / nu);
        target += student_t_lpdf(log_marginNormInv[i] | nu, pred_log_margin, scale_margin_marginal);
      }
    }
  }
}

generated quantities {
  vector[N] pred_shareIn = rep_vector(0.0, N);
  vector[N] pred_shareUncond = rep_vector(0.0, N);
  vector[N] pred_marginNormInv = rep_vector(0.0, N);
  vector[N] log_lik_share_marginal = rep_vector(0.0, N);
  vector[N] log_lik_margin = rep_vector(0.0, N);
  vector[N] log_lik = rep_vector(0.0, N);
  vector[N_market_year] log_lik_period = rep_vector(0.0, N_market_year);
  vector[N_market] log_lik_market = rep_vector(0.0, N_market);
  vector[N_cluster] log_lik_cluster = rep_vector(0.0, N_cluster);

  for (i in 1:N) {
    int m = market[i];
    int mky = market_year_idx[i];
    real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
    real w = conduct_weight_ces(firm_revShare[i], cutoff, use_cutoff);
    real b_quality = b_quality_all[group[i]] + w * mu_b_fringe;
    real delta = b_market[m] + b_quality
      + market_sign * (a_market[m] + log_price_sd) * log_price_std[i]
      + gamma_structural_market[m] * log(outside_price[mky]);
    real pred_log_q = logit(s0_prob[mky]) + delta;
    real sigma_share_eff = sigma_share_market[m];
    real sigma_margin_eff = sigma_margin_market[m];
    real pred_inv_margin = ces_inv_margin(
      firm_revShare[i], R0_prob[mky], gamma_structural_market[m], w,
      supply_model, bargaining_power
    );
    real pred_log_margin = log(pred_inv_margin + 1e-9);

    // pred_shareIn is the predicted *conditional inside* share (matches
    // log_shareIn); s0_prob[mky]*pred_shareIn is the predicted
    // *unconditional* share, not pred_shareIn itself.
    pred_shareIn[i] = exp(pred_log_q);
    pred_shareUncond[i] = s0_prob[mky] * pred_shareIn[i];
    pred_marginNormInv[i] = pred_inv_margin;

    if (likelihood_mode != 3) {
      real scale_share_marginal = sigma_share_eff * sqrt((nu - 2) / nu);
      log_lik_share_marginal[i] = student_t_lpdf(
        log_shareIn[i] | nu, pred_log_q, scale_share_marginal
      );
      log_lik[i] += log_lik_share_marginal[i];
    }

    if (likelihood_mode != 2 && margin_obs[i] == 1) {
      if (rho_active == 1) {
        real err_share = log_shareIn[i] - pred_log_q;
        real z_s = err_share / sigma_share_eff;
        real slope_cond = rho * (sigma_margin_eff / sigma_share_eff);
        real mu2_cond = pred_log_margin + slope_cond * err_share;
        real scale2_cond = sigma_margin_eff * sqrt(1 - square(rho))
          * sqrt((nu - 2 + square(z_s)) / (nu + 1));
        log_lik_margin[i] = student_t_lpdf(
          log_marginNormInv[i] | nu + 1, mu2_cond, scale2_cond
        );
      } else {
        real scale_margin_marginal = sigma_margin_eff * sqrt((nu - 2) / nu);
        log_lik_margin[i] = student_t_lpdf(
          log_marginNormInv[i] | nu, pred_log_margin, scale_margin_marginal
        );
      }
      log_lik[i] += log_lik_margin[i];
    }

    log_lik_period[mky] += log_lik[i];
    log_lik_market[m] += log_lik[i];
    log_lik_cluster[cluster[i]] += log_lik[i];
  }
}
