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
                      vector marginInv, vector shareIn, vector firm_shareIn, vector log_shareIn,
                      vector log_marginInv, vector price_std,
                      array[] int market, array[] int group, array[] int year,
                      array[] int market_year_idx,
                      vector a_market, vector alpha_structural_market, vector b_market, vector b_quality_all,
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

      real s_out = s0_prob[mky];

      real w = 0;
      if (use_cutoff == 1) {
        real x = -(shareIn[i] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[i]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality +
        (market_sign * a_market[m] * price_std[i]);

      // shareIn is normalized over inside goods. For structural outside share
      // s_out, log(s_j | inside) = logit(s_out) + delta_j.
      real pred_log_s = logit(s_out) + delta;

      real alpha_structural = alpha_structural_market[m];

      real inv_markup;
      if (supply_model == 1) {
        // Bertrand: margin = 1 / (alpha * (1 - s_uncond))
        inv_markup = alpha_structural * (w + (1-w) *
          (1.0 - firm_shareIn[i] * (1.0 - s_out)));
      } else if (supply_model == 2) {
        // Monopolistic competition: alpha = 1/(p-c)
        inv_markup = alpha_structural;
      } else if (supply_model == 3) {
        // Cournot: margin = (s0 + s) / (alpha * s0)
        inv_markup = alpha_structural * (w + (1-w) *
          (s_out / (s_out + firm_shareIn[i] * (1.0 - s_out))));
      } else {
        // 2nd Score Auction: margin = -log(1 - s) / (alpha * s)
        real denom = fmin(0.99, firm_shareIn[i] * (1.0 - s_out));
        inv_markup = alpha_structural * (w + (1-w) * (denom / -log1m(denom)));
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
  vector[N] firm_shareIn;
  vector[N] marginInv;
  vector[N] log_shareIn;
  vector[N] log_marginInv;
  vector[N] price_std;
  real<lower=0> price_sd;
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
  int is_cournot = (supply_model == 3) ? 1 : 0;
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
  vector[N_market] gamma_market;
  vector[N_group] b_group_raw;
  real<lower=0> sigma_b_strat;
  vector[use_cutoff == 1 ? 1 : 0] mu_b_fringe_vec;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 0 : 1] mu_log_a_free;
  vector<lower=0>[pool_alpha_by_market == 1 ? 1 : 0] sigma_log_a_vec;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 1 : 0] log_kappa_vec;
  vector[(is_cournot == 1 && pool_alpha_by_market == 1) ? N_market : 0] log_kappa_market;
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
  real mu_log_a;
  vector[N_market] a_market;
  vector[N_market] log_a_market_derived = rep_vector(0.0, N_market);
  vector[N_market] b_market;
  vector[N_group] b_quality_all = sigma_b_strat * (b_group_raw - mean(b_group_raw));
  real mu_b_fringe = 0.0;
  vector[N_market_year] s0;
  vector[N_market_year] s0_prob;
  vector[N_market] mean_log_s0_market = rep_vector(0.0, N_market);
  vector[N_market] mean_logit_s0_market = rep_vector(0.0, N_market);
  vector[N_market] alpha_structural_market;

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

  // Rotate the market intercept by mean logit(s0). This makes the level of
  // conditional inside shares neutral to the outside-share level while
  // retaining relative market-year movements in panel applications.
  for (mky in 1:N_market_year) {
    int mm = mky_to_market[mky];
    mean_log_s0_market[mm] += log(s0_prob[mky]);
    mean_logit_s0_market[mm] += logit(s0_prob[mky]);
  }
  for (mm in 1:N_market) {
    mean_log_s0_market[mm] /= n_mky_per_market[mm];
    mean_logit_s0_market[mm] /= n_mky_per_market[mm];
  }
  b_market = gamma_market - mean_logit_s0_market;

  if (is_cournot == 0) {
    mu_log_a = mu_log_a_free[1];
    if (is_single_market == 1 || use_alpha_market_pool == 0) {
      a_market = rep_vector(exp(mu_log_a), N_market);
    } else if (use_centered_alpha == 1) {
      a_market = exp(log_a_market_pool);
    } else {
      a_market = exp(mu_log_a + sigma_log_a * r_market_a_raw);
    }
  } else if (pool_alpha_by_market == 0) {
    real log_s0_ref_global = mean(log(s0_prob));
    mu_log_a = log_kappa_vec[1] - log_s0_ref_global;
    a_market = rep_vector(exp(mu_log_a), N_market);
  } else {
    log_a_market_derived = log_kappa_market - mean_log_s0_market;
    a_market = exp(log_a_market_derived);
    mu_log_a = mu_log_a_free[1];
  }

  alpha_structural_market = a_market / price_sd;
}

model {
  if (is_cournot == 1 && pool_alpha_by_market == 0) {
    mu_log_a ~ normal(prior_alpha_mean, prior_alpha_sd);
  } else {
    mu_log_a_free[1] ~ normal(prior_alpha_mean, prior_alpha_sd);
  }
  if (pool_alpha_by_market == 1) {
    sigma_log_a_vec ~ normal(0, prior_sigma_alpha);
  }
  if (is_cournot == 0 && pool_alpha_by_market == 1) {
    if (use_centered_alpha == 1) {
      log_a_market_pool ~ normal(mu_log_a, sigma_log_a);
    } else {
      r_market_a_raw ~ std_normal();
    }
  }
  if (is_cournot == 1 && pool_alpha_by_market == 1) {
    log_a_market_derived ~ normal(mu_log_a, sigma_log_a);
  }
  if (is_single_market == 1) {
    gamma_market ~ normal(0, 1.0);
  } else {
    mu_b_market_vec ~ normal(0, 1.0);
    sigma_b_market_vec ~ normal(0, 0.5);
    gamma_market ~ normal(mu_b_market_vec[1], sigma_b_market_vec[1]);
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
    // A two-dimensional LKJ prior has density proportional to
    // (1-rho^2)^(eta-1). Include the transform Jacobian for the bounded
    // rho = 0.95*tanh(rho_raw) parameterization.
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
                       marginInv, shareIn, firm_shareIn, log_shareIn, log_marginInv,
                       price_std, market, group, year,
                       market_year_idx,
                       a_market, alpha_structural_market, b_market, b_quality_all, sigma_b_strat,
                       mu_b_fringe,
                       s0_prob, cutoff_share,
                       sigma_share_market, sigma_margin_market,
                       rho, supply_model, market_sign,
                       use_cutoff, is_single_market);
}

generated quantities {
  vector[N] pred_shareIn;
  vector[N] pred_marginInv;
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
        real x = -(shareIn[n] - cutoff_share[m]) / k;
        if (x > 10) w = 1.0; else if (x < -10) w = 0.0; else w = inv_logit(x);
      }

      real b_quality = b_quality_all[group[n]] + (w * mu_b_fringe);
      real delta = b_market[m] + b_quality +
        (market_sign * a_market[m] * price_std[n]);

      real s_out = s0_prob[mky];
      real pred_log_s = logit(s_out) + delta;
      pred_shareIn[n] = exp(pred_log_s);

      real alpha_structural = alpha_structural_market[m];

      real inv_markup;
      if (supply_model == 1) {
        // Bertrand
        inv_markup = alpha_structural * (w + (1-w) *
          (1.0 - firm_shareIn[n] * (1.0 - s_out)));
      } else if (supply_model == 2) {
        // Monopolistic competition: alpha = 1/(p-c)
        inv_markup = alpha_structural;
      } else if (supply_model == 3) {
        // Cournot
        inv_markup = alpha_structural * (w + (1-w) *
          (s_out / (s_out + firm_shareIn[n] * (1.0 - s_out))));
      } else {
        // 2nd Score Auction
        real denom = fmin(0.99, firm_shareIn[n] * (1.0 - s_out));
        inv_markup = alpha_structural * (w + (1-w) * (denom / -log1m(denom)));
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
