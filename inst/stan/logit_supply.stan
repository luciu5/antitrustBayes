// Specification 9: centered orthonormal firm-quality coordinates.
// Bayesian Logit demand and supply model.
// Conditional inside shares and inverse normalized margins are coupled
// per observation through a free cross-equation correlation rho, using a
// bivariate Student-t factored as p(share) * p(margin | share). The
// latent outside share s0 is identified jointly through the share equation
// (log_shareIn ~ logit(s0) + delta) and the equilibrium margin condition.
// likelihood_mode == 2 (share-only) omits the margin equation entirely, so
// it lacks the equilibrium-margin information that principally identifies
// s0's level; s0 may be weakly identified or prior-dominated in that mode,
// depending on market/time/intercept structure. This is an expected
// property of the model and is not compensated for with extra pooling.
// Supports Bertrand, monopolistic competition, Cournot, second-score auction,
// and Nash bargaining conduct.

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

  real conduct_weight_logit(real firm_share_inside, real cutoff, int use_cutoff, real cutoff_scale) {
    real x;
    if (use_cutoff == 0) return 0.0;
    x = -(firm_share_inside - cutoff) / cutoff_scale;
    return inv_logit(x);
  }

  real logit_inv_markup(real firm_share_inside, real s0, real alpha,
                         real fringe_weight, int supply_model,
                         real bargaining_power) {
    real firm_share = firm_share_inside * (1.0 - s0);
    real conduct;
    if (supply_model == 1) {
      conduct = 1.0 - firm_share;
    } else if (supply_model == 2) {
      conduct = 1.0;
    } else if (supply_model == 3) {
      conduct = s0 / (s0 + firm_share);
    } else if (supply_model == 4) {
      if (firm_share < 1e-10) conduct = 1.0;
      else conduct = firm_share / -log1m(firm_share);
    } else {
      if (firm_share < 1e-10) {
        conduct = 1.0;
      } else {
        real kappa = bargaining_power / (1.0 - bargaining_power);
        conduct = (kappa * firm_share
          - (1.0 - firm_share) * log1m(firm_share)) / -log1m(firm_share);
      }
    }
    return alpha * (fringe_weight + (1.0 - fringe_weight) * conduct);
  }

  // Stable weighted within-group regression statistics avoid cancellation
  // when share noise is very small.  Every quadratic contribution is positive.
  real centered_gaussian_from_stats(real within, vector a, vector c,
                                    real logdet_v_diag, int L, int G,
                                    real tau) {
    vector[G] w;
    real logdet_v = logdet_v_diag;
    real q0 = within;
    real b = 0.0;
    real den = 0.0;
    real q;
    if (G == 1) {
      if (a[1] > 0) q0 += square(c[1] / a[1]) * a[1];
      return -0.5 * (L * log(2.0 * pi()) + logdet_v_diag + q0);
    }
    for (j in 1:G) {
      w[j] = 1.0 / (1.0 + square(tau) * a[j]);
      logdet_v += log1p(square(tau) * a[j]);
      den += w[j];
      b += c[j] * w[j];
      if (a[j] > 0) q0 += square(c[j] / a[j]) * a[j] * w[j];
    }
    q = q0 + square(tau) * square(b) / den;
    return -0.5 * (L * log(2.0 * pi()) + logdet_v + log(den / G) + q);
  }

  // Normalized O(L+G) density; one group loading per observation.
  real centered_gaussian_density(vector z, vector x, vector v,
                              array[] int g, int G, real tau) {
    int L = rows(z);
    vector[G] a = rep_vector(0.0, G);
    vector[G] c = rep_vector(0.0, G);
    real within = 0.0;
    real logdet_v = 0.0;
    for (l in 1:L) {
      real weight = square(x[l]) / v[l];
      if (a[g[l]] > 0) {
        real residual = z[l] - x[l] * c[g[l]] / a[g[l]];
        within += (a[g[l]] / (a[g[l]] + weight)) * square(residual) / v[l];
      } else if (x[l] == 0) within += square(z[l]) / v[l];
      a[g[l]] += weight;
      c[g[l]] += x[l] * z[l] / v[l];
      logdet_v += log(v[l]);
    }
    return centered_gaussian_from_stats(within, a, c, logdet_v, L, G, tau);
  }

  real collapsed_share_gaussian_density(vector y, vector mu,
                                     array[] int g, int G, real sigma,
                                     real tau) {
    int L = rows(y);
    vector[L] z;
    vector[L] x = rep_vector(1.0, L);
    vector[L] v = rep_vector(square(sigma), L);
    for (l in 1:L) z[l] = y[l] - mu[l];
    return centered_gaussian_density(z, x, v, g, G, tau);
  }

  real collapsed_joint_gaussian_density(vector share_resid,
                                     vector margin_resid,
                                     array[] int margin_idx,
                                     array[] int g, int G,
                                     real sigma_share, real sigma_margin,
                                     real rho, real tau) {
    int N0 = rows(share_resid);
    int Q = rows(margin_resid);
    int L = N0 + Q;
    vector[L] z;
    vector[L] x;
    vector[L] v;
    real slope = rho * sigma_margin / sigma_share;
    for (i in 1:N0) {
      z[i] = share_resid[i];
      x[i] = 1.0;
      v[i] = square(sigma_share);
    }
    for (j in 1:Q) {
      int i = margin_idx[j];
      z[N0 + j] = margin_resid[j] - slope * share_resid[i];
      x[N0 + j] = -slope;
      v[N0 + j] = square(sigma_margin) * (1.0 - square(rho));
    }
    {
      array[L] int g_aug;
      for (i in 1:N0) g_aug[i] = g[i];
      for (j in 1:Q) g_aug[N0 + j] = g[margin_idx[j]];
      return centered_gaussian_density(z, x, v, g_aug, G, tau);
    }
  }
}

data {
  int<lower=1> N;
  vector[N] shareIn;
  vector[N] firm_shareIn;
  vector[N] clr_shareIn;
  vector[N] log_shareIn;
  vector[N] log_marginNormInv;
  vector[N] price;
  vector[N] price_std;
  real<lower=0> price_sd;
  int<lower=-1, upper=1> market_sign;
  int<lower=1> N_market;
  array[N] int<lower=1, upper=N_market> market;
  int<lower=1> N_group;
  array[N] int<lower=1, upper=N_group> group;
  int<lower=1> N_year;
  array[N] int<lower=1, upper=N_year> year;
  int<lower=1, upper=5> supply_model;
  int<lower=0, upper=1> use_cutoff;
  real<lower=0> cutoff_scale;
  int<lower=1> N_market_year;
  array[N] int<lower=1, upper=N_market_year> market_year_idx;
  array[N_market_year] int<lower=1> mky_n;
  array[N_market_year] int<lower=1> mky_start;
  array[N] int<lower=1, upper=N> obs_order;
  int<lower=1> N_cluster;
  array[N] int<lower=1, upper=N_cluster> cluster;
  int<lower=1, upper=3> likelihood_mode;
  array[N] int<lower=0, upper=1> margin_obs;
  int<lower=0> N_margin;
  array[N_margin] int<lower=1, upper=N> margin_idx;
  int<lower=1, upper=2> share_likelihood;
  int<lower=1, upper=2> margin_likelihood;
  int<lower=1, upper=2> quality_parameterization;
  int<lower=0> K_s0;
  matrix[N_market_year, K_s0] X_s0;
  array[N_market_year] int<lower=1, upper=N_market> mky_to_market;
  array[N_market_year] int<lower=1, upper=N_year> mky_to_year;
  int<lower=0, upper=1> is_single_market;
  int<lower=0, upper=1> use_rho;
  int<lower=0, upper=1> use_market_resid_scale;
  int<lower=0, upper=1> use_alpha_market_pool;
  int<lower=0, upper=1> use_s0_market_year_pool;
  int<lower=0, upper=1> use_s0_market_pool;
  int<lower=1, upper=2> market_s0_parameterization;
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
  int is_cournot = supply_model == 3;
  int is_bargaining = supply_model == 5;
  int pool_alpha_by_market = is_single_market == 0 && use_alpha_market_pool == 1;
  array[N_market] int n_mky_per_market = rep_array(0, N_market);
  int any_paired = 0;
  int rho_active;
  int market_s0_rotation_active =
    (market_s0_parameterization == 2 && is_single_market == 0) ? 1 : 0;
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
  vector[is_single_market == 0 ? N_market : 0] b_market_raw;
  vector[quality_parameterization == 1 ? N_group - 1 : 0] b_group_centered;
  real<lower=0> sigma_b_strat;
  vector[use_cutoff == 1 ? 1 : 0] mu_b_fringe_vec;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 0 : 1] mu_log_a_free;
  vector<lower=0>[pool_alpha_by_market == 1 ? 1 : 0] sigma_log_a_vec;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 1 : 0] log_kappa_vec;
  vector[(is_cournot == 1 && pool_alpha_by_market == 1) ? N_market : 0] log_kappa_market;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1 ||
          market_s0_rotation_active == 1) ? 0 : N_market_year] s0_logit;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1) ? 1 : 0] mu_s0_pool_vec;
  vector<lower=0>[(use_s0_market_year_pool == 1 ||
                  (use_s0_market_pool == 1 && is_single_market == 0)) ? 1 : 0] sigma_s0_pool_vec;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 0 &&
          market_s0_rotation_active == 0) ? N_market_year : 0] s0_pool_raw;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 1 &&
          market_s0_rotation_active == 0) ? N_market_year : 0] eta_s0_pool;
  // Market-level s0 pooling: one outside share per market, shared across the
  // market's years (structurally parallel to a_market). Fed by N_market group
  // effects each identified by all of the market's observations, so its hyper-SD
  // is far better identified than the market-year pool's N_market_year effects.
  vector[(use_s0_market_pool == 1 && is_single_market == 0 &&
          use_centered_s0_pool == 0 && market_s0_rotation_active == 0) ? N_market : 0] s0_market_raw;
  vector[(use_s0_market_pool == 1 && is_single_market == 0 &&
          use_centered_s0_pool == 1 && market_s0_rotation_active == 0) ? N_market : 0] eta_s0_market;
  vector[(market_s0_rotation_active == 1 &&
          (use_centered_s0_pool == 1 ||
           (use_s0_market_year_pool == 0 && use_s0_market_pool == 0)))
         ? (use_s0_market_pool == 1 ? N_market : N_market_year) : 0] c_s0;
  vector[(market_s0_rotation_active == 1 && use_centered_s0_pool == 0 &&
          (use_s0_market_year_pool == 1 || use_s0_market_pool == 1))
         ? (use_s0_market_pool == 1 ? N_market : N_market_year) : 0] c_s0_raw;
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
  real mu_log_a;
  vector[N_market] a_market;
  vector[N_market] log_a_market_derived = rep_vector(0.0, N_market);
  vector[N_market] b_market;
  vector[N_market] gamma_market;
  vector[N_market_year] s0_base_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] s0_raw_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_unbounded_mky = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_effective_mky = rep_vector(0.0, N_market_year);
  vector[quality_parameterization == 1 ? N_group : 0] b_quality_all =
    rep_vector(0.0, quality_parameterization == 1 ? N_group : 0);
  vector[N] share_mean = rep_vector(0.0, N);
  real mu_b_fringe = 0.0;
  vector[N_market_year] s0;
  vector[N_market_year] s0_prob;
  vector[N_market] mean_log_s0_market = rep_vector(0.0, N_market);
  vector[N_market] mean_logit_s0_market = rep_vector(0.0, N_market);
  vector[N_market] alpha_structural_market;
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

  // Reconstruct the unchanged market effect before applying the outside-share
  // rotation. Single-market fits retain the historical zero normalization.
  if (is_single_market == 1) {
    b_market = rep_vector(0.0, N_market);
    gamma_market = rep_vector(0.0, N_market);
  } else {
    b_market = mu_b_market_vec[1]
      + sigma_b_market_vec[1] * b_market_raw;
  }

  for (mky in 1:N_market_year) {
    real base_logit;
    real x_beta = K_s0 > 0 ? X_s0[mky] * beta_s0 : 0.0;
    int m = mky_to_market[mky];
    if (use_s0_market_pool == 1) {
      // One s0 per market, applied to each of the market's market-years.
      if (is_single_market == 1) base_logit = mu_s0_pool_vec[1];
      else if (market_s0_rotation_active == 1) {
        if (use_centered_s0_pool == 1)
          base_logit = c_s0[m] - b_market[m];
        else {
          real w = c_s0_raw[m];
          base_logit = mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * w
            - b_market[m];
          s0_raw_reconstructed[mky] = w - b_market[m] / sigma_s0_pool_vec[1];
        }
      } else base_logit = use_centered_s0_pool == 1
        ? eta_s0_market[m]
        : mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * s0_market_raw[m];
    } else if (use_s0_market_year_pool == 1) {
      if (market_s0_rotation_active == 1) {
        if (use_centered_s0_pool == 1)
          base_logit = c_s0[mky] - b_market[m];
        else {
          real w = c_s0_raw[mky];
          base_logit = mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * w
            - b_market[m];
          s0_raw_reconstructed[mky] = w - b_market[m] / sigma_s0_pool_vec[1];
        }
      } else base_logit = use_centered_s0_pool == 1
        ? eta_s0_pool[mky]
        : mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * s0_pool_raw[mky];
    } else {
      base_logit = market_s0_rotation_active == 1
        ? c_s0[mky] - b_market[m] : s0_logit[mky];
    }
    if (market_s0_rotation_active == 0) {
      if (is_single_market == 0 && use_s0_market_pool == 1 &&
          use_centered_s0_pool == 0)
        s0_raw_reconstructed[mky] = s0_market_raw[m];
      else if (use_s0_market_year_pool == 1 && use_centered_s0_pool == 0)
        s0_raw_reconstructed[mky] = s0_pool_raw[mky];
    }
    s0_base_reconstructed[mky] = base_logit;
    s0[mky] = base_logit + x_beta;
    c_unbounded_mky[mky] = s0[mky] + b_market[m];
    s0_prob[mky] = bounded_inv_logit(s0[mky]);
    c_effective_mky[mky] = logit(s0_prob[mky]) + b_market[m];
    mean_log_s0_market[mky_to_market[mky]] += log(s0_prob[mky]);
    mean_logit_s0_market[mky_to_market[mky]] += logit(s0_prob[mky]);
  }
  for (m in 1:N_market) {
    mean_log_s0_market[m] /= n_mky_per_market[m];
    mean_logit_s0_market[m] /= n_mky_per_market[m];
    if (is_single_market == 0) {
      if (n_mky_per_market[m] > 1)
        gamma_market[m] = b_market[m] + mean_logit_s0_market[m];
      else gamma_market[m] = b_market[m];
    }
  }

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
  if (quality_parameterization == 1 && N_group > 1) {
    real suffix = 0.0;
    for (reverse_j in 1:(N_group - 1)) {
      int j = N_group + 1 - reverse_j;
      int k = j - 1;
      real v_k = b_group_centered[k] / sqrt(1.0 * k * (k + 1.0));
      b_quality_all[j] = suffix - k * v_k;
      suffix += v_k;
    }
    b_quality_all[1] = suffix;
  }
  for (i in 1:N) {
    int m = market[i];
    int mky = market_year_idx[i];
    real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
    real w = conduct_weight_logit(firm_shareIn[i], cutoff, use_cutoff, cutoff_scale);
    real b_quality = (quality_parameterization == 1 ? b_quality_all[group[i]] : 0.0)
      + w * mu_b_fringe;
    share_mean[i] = logit(s0_prob[mky]) + b_market[m] + b_quality
      + market_sign * a_market[m] * price_std[i];
  }
}

model {
  real collapsed_share_ll = 0.0;
  if (is_cournot == 1 && pool_alpha_by_market == 0) {
    target += normal_lpdf(mu_log_a | prior_alpha_mean, prior_alpha_sd);
  } else {
    target += normal_lpdf(mu_log_a_free[1] | prior_alpha_mean, prior_alpha_sd);
  }
  if (pool_alpha_by_market == 1) target += normal_lpdf(sigma_log_a_vec | 0, prior_sigma_alpha) - num_elements(sigma_log_a_vec) * normal_lccdf(0 | 0, prior_sigma_alpha);
  if (is_cournot == 0 && pool_alpha_by_market == 1) {
    if (use_centered_alpha == 1) target += normal_lpdf(log_a_market_pool | mu_log_a, sigma_log_a);
    else target += std_normal_lpdf(r_market_a_raw);
  }
  if (is_cournot == 1 && pool_alpha_by_market == 1) {
    target += normal_lpdf(log_a_market_derived | mu_log_a, sigma_log_a);
  }
  if (is_single_market == 0) {
    target += normal_lpdf(mu_b_market_vec | 0, 1);
    target += normal_lpdf(sigma_b_market_vec | 0, 0.5) - num_elements(sigma_b_market_vec) * normal_lccdf(0 | 0, 0.5);
    target += std_normal_lpdf(b_market_raw);
  }
  target += normal_lpdf(b_group_centered | 0, sigma_b_strat);
  target += normal_lpdf(sigma_b_strat | 1, prior_sigma_meanval_strat) - 1 * normal_lccdf(0 | 1, prior_sigma_meanval_strat);
  if (use_cutoff == 1) {
    target += beta_lpdf(cutoff_share | prior_cutoff_alpha, prior_cutoff_beta);
    target += normal_lpdf(mu_b_fringe_vec | 0, prior_sigma_meanval_fringe);
  }
  if (use_s0_market_pool == 1) {
    target += normal_lpdf(mu_s0_pool_vec | s0_prior_mean, s0_prior_sd);
    if (is_single_market == 0) {
      target += normal_lpdf(sigma_s0_pool_vec | 0, prior_sigma_s0_pool) - num_elements(sigma_s0_pool_vec) * normal_lccdf(0 | 0, prior_sigma_s0_pool);
      if (market_s0_rotation_active == 1) {
        if (use_centered_s0_pool == 1)
          target += normal_lpdf(c_s0 - b_market | mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
        else {
          vector[N_market] z_s0 = c_s0_raw - b_market / sigma_s0_pool_vec[1];
          target += std_normal_lpdf(z_s0);
        }
      } else if (use_centered_s0_pool == 1) {
        target += normal_lpdf(eta_s0_market | mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
      } else target += std_normal_lpdf(s0_market_raw);
    }
  } else if (use_s0_market_year_pool == 1) {
    target += normal_lpdf(mu_s0_pool_vec | s0_prior_mean, s0_prior_sd);
    target += normal_lpdf(sigma_s0_pool_vec | 0, prior_sigma_s0_pool) - num_elements(sigma_s0_pool_vec) * normal_lccdf(0 | 0, prior_sigma_s0_pool);
    if (market_s0_rotation_active == 1) {
      if (use_centered_s0_pool == 1) {
        vector[N_market_year] v_s0;
        for (mky in 1:N_market_year)
          v_s0[mky] = c_s0[mky] - b_market[mky_to_market[mky]];
        target += normal_lpdf(v_s0 | mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
      } else {
        vector[N_market_year] z_s0;
        for (mky in 1:N_market_year)
          z_s0[mky] = c_s0_raw[mky]
            - b_market[mky_to_market[mky]] / sigma_s0_pool_vec[1];
        target += std_normal_lpdf(z_s0);
      }
    } else if (use_centered_s0_pool == 1) {
      target += normal_lpdf(eta_s0_pool | mu_s0_pool_vec[1], sigma_s0_pool_vec[1]);
    } else target += std_normal_lpdf(s0_pool_raw);
  } else if (market_s0_rotation_active == 1) {
    vector[N_market_year] v_s0;
    for (mky in 1:N_market_year)
      v_s0[mky] = c_s0[mky] - b_market[mky_to_market[mky]];
    target += normal_lpdf(v_s0 | s0_prior_mean, s0_prior_sd);
  } else target += normal_lpdf(s0_logit | s0_prior_mean, s0_prior_sd);
  target += normal_lpdf(beta_s0 | 0, prior_sigma_beta_s0);
  if (is_bargaining == 1) {
    target += beta_lpdf(bargaining_unit | bargaining_prior_alpha, bargaining_prior_beta);
  }
  target += normal_lpdf(log_sigma_share | log(prior_sigma_share), prior_log_sigma_share_sd);
  target += normal_lpdf(log_sigma_margin | log(prior_sigma_margin), prior_log_sigma_margin_sd);
  if (use_market_resid_scale == 1) {
    target += std_normal_lpdf(log_sigma_share_market_raw);
    target += std_normal_lpdf(log_sigma_margin_market_raw);
    target += normal_lpdf(sigma_log_sigma_share_market | 0, prior_sigma_resid_market) - num_elements(sigma_log_sigma_share_market) * normal_lccdf(0 | 0, prior_sigma_resid_market);
    target += normal_lpdf(sigma_log_sigma_margin_market | 0, prior_sigma_resid_market) - num_elements(sigma_log_sigma_margin_market) * normal_lccdf(0 | 0, prior_sigma_resid_market);
  }
  if (estimate_df == 1) target += gamma_lpdf(nu_minus_2 | prior_nu_shape, prior_nu_rate);
  if (rho_active == 1) {
    // Two-dimensional LKJ-style prior, with the Jacobian for the bounded
    // rho = 0.95*tanh(rho_raw) parameterization.
    target += (prior_lkj - 1.0) * log1m(square(rho))
      + log1m(square(tanh(rho_raw[1])))
      - (lbeta(0.5, prior_lkj) + beta_lcdf(square(0.95) | 0.5, prior_lkj) - log(0.95));
  }

  if (use_hmt == 1 && likelihood_mode != 2) {
    real alpha_check = exp(mu_log_a) / price_sd;
    real s0_check = mean(s0_prob);
    real agg_elasticity = abs(alpha_check * avg_price_hmt * s0_check);
    real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
    target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  if (likelihood_mode != 3) {
    if (quality_parameterization == 2 && share_likelihood == 2) {
      collapsed_share_ll = collapsed_share_gaussian_density(
        log_shareIn, share_mean, group, N_group, sigma_share, sigma_b_strat);
      target += collapsed_share_ll;
    } else for (i in 1:N) {
      int m = market[i];
      if (share_likelihood == 2)
        target += normal_lpdf(log_shareIn[i] | share_mean[i], sigma_share_market[m]);
      else
        target += student_t_lpdf(log_shareIn[i] | nu, share_mean[i],
          sigma_share_market[m] * sqrt((nu - 2) / nu));
    }
  }

  if (likelihood_mode != 2) {
    if (quality_parameterization == 2 && likelihood_mode != 3 && share_likelihood == 2 &&
        margin_likelihood == 2 && rho_active == 1) {
      vector[N_margin] margin_resid;
      for (q in 1:N_margin) {
        int i = margin_idx[q];
        int m = market[i];
        int mky = market_year_idx[i];
        real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
        real w = conduct_weight_logit(firm_shareIn[i], cutoff, use_cutoff, cutoff_scale);
        real inv_markup = logit_inv_markup(firm_shareIn[i], s0_prob[mky],
          alpha_structural_market[m], w, supply_model, bargaining_power);
        margin_resid[q] = log_marginNormInv[i] - log(price[i] * inv_markup + 1e-9);
      }
      target += collapsed_joint_gaussian_density(
        log_shareIn - share_mean, margin_resid, margin_idx, group, N_group,
        sigma_share, sigma_margin, rho, sigma_b_strat) - collapsed_share_ll;
    } else for (i in 1:N) if (margin_obs[i] == 1) {
      int m = market[i];
      int mky = market_year_idx[i];
      real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
      real w = conduct_weight_logit(firm_shareIn[i], cutoff, use_cutoff, cutoff_scale);
      real inv_markup = logit_inv_markup(
        firm_shareIn[i], s0_prob[mky], alpha_structural_market[m], w,
        supply_model, bargaining_power
      );
      real pred_log_margin = log(price[i] * inv_markup + 1e-9);
      if (margin_likelihood == 2) {
        if (rho_active == 1) {
          real err_share = log_shareIn[i] - share_mean[i];
          target += normal_lpdf(log_marginNormInv[i] |
            pred_log_margin + rho * (sigma_margin_market[m] / sigma_share_market[m]) * err_share,
            sigma_margin_market[m] * sqrt(1 - square(rho)));
        } else target += normal_lpdf(log_marginNormInv[i] | pred_log_margin,
          sigma_margin_market[m]);
      } else if (rho_active == 1) {
        real err_share = log_shareIn[i] - share_mean[i];
        real z_s = err_share / sigma_share_market[m];
        real scale2_cond = sigma_margin_market[m] * sqrt(1 - square(rho))
          * sqrt((nu - 2 + square(z_s)) / (nu + 1));
        target += student_t_lpdf(log_marginNormInv[i] | nu + 1,
          pred_log_margin + rho * (sigma_margin_market[m] / sigma_share_market[m]) * err_share,
          scale2_cond);
      } else target += student_t_lpdf(log_marginNormInv[i] | nu, pred_log_margin,
        sigma_margin_market[m] * sqrt((nu - 2) / nu));
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
  real log_lik_share_collapsed = 0.0;
  real log_lik_joint_collapsed = 0.0;

  if (quality_parameterization == 2) {
    vector[N] margin_mean = rep_vector(0.0, N);
    vector[N_margin] margin_resid;
    vector[N_margin] margin_values;
    real share_ll = 0.0;
    real joint_ll = 0.0;
    for (i in 1:N) {
      int m = market[i];
      int mky = market_year_idx[i];
      real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
      real w = conduct_weight_logit(firm_shareIn[i], cutoff, use_cutoff, cutoff_scale);
      real inv_markup = logit_inv_markup(firm_shareIn[i], s0_prob[mky],
        alpha_structural_market[m], w, supply_model, bargaining_power);
      pred_shareIn[i] = exp(share_mean[i]);
      pred_shareUncond[i] = (1.0 - s0_prob[mky]) * pred_shareIn[i];
      pred_marginNormInv[i] = price[i] * inv_markup;
      margin_mean[i] = log(pred_marginNormInv[i] + 1e-9);
    }
    if (likelihood_mode != 3)
      share_ll = collapsed_share_gaussian_density(
        log_shareIn, share_mean, group, N_group, sigma_share, sigma_b_strat);
    log_lik_share_collapsed = share_ll;
    if (likelihood_mode != 2) {
      if (margin_likelihood == 2) {
        for (q in 1:N_margin) {
          int i = margin_idx[q];
          margin_values[q] = log_marginNormInv[i];
          margin_resid[q] = margin_values[q] - margin_mean[i];
        }
        if (likelihood_mode == 3) {
          joint_ll = 0.0;
          for (q in 1:N_margin)
            joint_ll += normal_lpdf(margin_values[q] | margin_mean[margin_idx[q]], sigma_margin);
        } else {
          joint_ll = collapsed_joint_gaussian_density(
            log_shareIn - share_mean, margin_resid, margin_idx, group, N_group,
            sigma_share, sigma_margin, rho, sigma_b_strat);
        }
      } else {
        joint_ll = share_ll;
        for (q in 1:N_margin) {
          int i = margin_idx[q];
          joint_ll += student_t_lpdf(log_marginNormInv[i] | nu, margin_mean[i],
            sigma_margin * sqrt((nu - 2) / nu));
        }
      }
    } else joint_ll = share_ll;
    log_lik_joint_collapsed = joint_ll;
  } else {
    for (i in 1:N) {
      int m = market[i];
      int mky = market_year_idx[i];
      real cutoff = use_cutoff == 1 ? cutoff_share[m] : 0.0;
      real w = conduct_weight_logit(firm_shareIn[i], cutoff, use_cutoff, cutoff_scale);
      real inv_markup = logit_inv_markup(
        firm_shareIn[i], s0_prob[mky], alpha_structural_market[m], w,
        supply_model, bargaining_power
      );
      pred_shareIn[i] = exp(share_mean[i]);
      pred_shareUncond[i] = (1.0 - s0_prob[mky]) * pred_shareIn[i];
      pred_marginNormInv[i] = price[i] * inv_markup;
      if (likelihood_mode != 3) {
        if (share_likelihood == 2)
          log_lik_share_marginal[i] = normal_lpdf(log_shareIn[i] |
            share_mean[i], sigma_share_market[m]);
        else log_lik_share_marginal[i] = student_t_lpdf(log_shareIn[i] | nu,
          share_mean[i], sigma_share_market[m] * sqrt((nu - 2) / nu));
        log_lik[i] += log_lik_share_marginal[i];
      }
      if (likelihood_mode != 2 && margin_obs[i] == 1) {
        real pred_log_margin = log(pred_marginNormInv[i] + 1e-9);
        if (margin_likelihood == 2) {
          if (rho_active == 1)
            log_lik_margin[i] = normal_lpdf(log_marginNormInv[i] |
              pred_log_margin + rho * (sigma_margin_market[m] / sigma_share_market[m])
                * (log_shareIn[i] - share_mean[i]),
              sigma_margin_market[m] * sqrt(1 - square(rho)));
          else log_lik_margin[i] = normal_lpdf(log_marginNormInv[i] |
            pred_log_margin, sigma_margin_market[m]);
        } else if (rho_active == 1) {
          real err = log_shareIn[i] - share_mean[i];
          real z = err / sigma_share_market[m];
          log_lik_margin[i] = student_t_lpdf(log_marginNormInv[i] | nu + 1,
            pred_log_margin + rho * (sigma_margin_market[m] / sigma_share_market[m]) * err,
            sigma_margin_market[m] * sqrt(1 - square(rho))
              * sqrt((nu - 2 + square(z)) / (nu + 1)));
        } else log_lik_margin[i] = student_t_lpdf(log_marginNormInv[i] | nu,
          pred_log_margin, sigma_margin_market[m] * sqrt((nu - 2) / nu));
        log_lik[i] += log_lik_margin[i];
      }
      log_lik_period[mky] += log_lik[i];
      log_lik_market[m] += log_lik[i];
      log_lik_cluster[cluster[i]] += log_lik[i];
    }
  }
}
