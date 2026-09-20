// Specification 10: scoped centered orthonormal firm-quality coordinates.
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

  // The exact programs use one supplied finite-state prior per market.  The
  // response statistic is phi for Stackelberg Bertrand and S_F for
  // Stackelberg Cournot; it is ignored for exact fringe states.
  real stack_logit_inv_margin(int i, int mky, int k, real s0, real alpha,
                              int conduct, int structure, real response,
                              array[] int firm_mky_start,
                              array[] int firm_mky_n, vector firm_mky_share,
                              array[] int firm_mky_rank,
                              array[] int obs_to_firm_mky, vector price) {
    int me = obs_to_firm_mky[i];
    real S = firm_mky_share[me] * (1 - s0);
    real h;
    if (structure == 1) {
      // BF/CF: k=0 is exact product-level MonCom (h=1).
      if (firm_mky_rank[me] > k) h = 1.0;
      else if (conduct == 1) h = 1 / (1 - S);
      else h = 1 + S / s0;
    } else if (conduct == 1) {
      if (firm_mky_rank[me] > k) h = 1 / (1 - S);
      else h = 1 / (1 - response * S);
    } else {
      if (firm_mky_rank[me] > k) h = 1 + S / s0;
      else h = 1 + S / (s0 + response);
    }
    // Return inverse Lerner in price units: alpha * p / h.
    return price[i] * alpha / h;
  }
  real stack_margin_density(real y, real pred_log, real log_share,
                            real pred_log_share, real sigma_share,
                            real sigma_margin, real rho, real nu,
                            int rho_active) {
    if (rho_active == 1) {
      real err = log_share - pred_log_share;
      real z = err / sigma_share;
      real mu = pred_log + rho * (sigma_margin / sigma_share) * err;
      real sc = sigma_margin * sqrt(1 - square(rho))
        * sqrt((nu - 2 + square(z)) / (nu + 1));
      return student_t_lpdf(y | nu + 1, mu, sc);
    }
    return student_t_lpdf(y | nu, pred_log, sigma_margin * sqrt((nu - 2) / nu));
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


  real collapsed_joint_from_share_stats(real share_within, vector share_a,
                                     vector share_c, real share_log_v,
                                     vector share_resid, vector margin_resid,
                                     array[] int margin_idx, array[] int g,
                                     int G, real sigma_share, real sigma_margin,
                                     real rho, real tau) {
    int N0 = rows(share_resid);
    int Q = rows(margin_resid);
    vector[N0] z = share_resid;
    vector[N0] x = rep_vector(1.0, N0);
    vector[N0] v = rep_vector(square(sigma_share), N0);
    real slope = rho * sigma_share / sigma_margin;
    real paired_variance = square(sigma_share) * (1.0 - square(rho));
    real margin_ll = normal_lpdf(margin_resid | 0, sigma_margin);
    for (j in 1:Q) {
      int i = margin_idx[j];
      z[i] -= slope * margin_resid[j];
      v[i] = paired_variance;
    }
    return margin_ll + centered_gaussian_density(z, x, v, g, G, tau);
  }

  // Market-local version.  margin_order contains positions in margin_idx;
  // share_n is explicit because share_resid is the full data vector.

  // Market-local margin-first factorization.  All share observations are
  // rebuilt once in market order; missing margins leave their rows at the
  // unconditional share variance.
  real collapsed_joint_from_local_share_stats(
      real share_within, vector share_a, vector share_c,
      real share_log_v, vector share_resid, vector margin_resid,
      array[] int margin_order, array[] int margin_idx, int margin_start,
      array[] int market_obs_order, array[] int market_obs_position,
      int obs_start, array[] int quality_local,
      int share_n, int margin_n, int G,
      real sigma_share, real sigma_margin, real rho, real tau) {
    vector[share_n] z;
    vector[share_n] x = rep_vector(1.0, share_n);
    vector[share_n] v = rep_vector(square(sigma_share), share_n);
    array[share_n] int g;
    real slope = rho * sigma_share / sigma_margin;
    real paired_variance = square(sigma_share) * (1.0 - square(rho));
    real margin_ll = normal_lpdf(margin_resid | 0, sigma_margin);
    for (pos in 1:share_n) {
      int i = market_obs_order[obs_start + pos - 1];
      z[pos] = share_resid[i];
      g[pos] = quality_local[i];
    }
    for (q in 1:margin_n) {
      int margin_pos = margin_order[margin_start + q - 1];
      int i = margin_idx[margin_pos];
      int pos = market_obs_position[i];
      z[pos] -= slope * margin_resid[q];
      v[pos] = paired_variance;
    }
    return margin_ll + centered_gaussian_density(z, x, v, g, G, tau);
  }


  // Exact margin-first factorization of the collapsed joint Gaussian.  It
  // avoids the equivalent margin-on-share coefficient rho*sigma_m/sigma_s,
  // which is unbounded as the share-noise scale approaches zero.
  real collapsed_joint_gaussian_density(vector share_resid,
                                     vector margin_resid,
                                     array[] int margin_idx,
                                     array[] int g, int G,
                                     real sigma_share, real sigma_margin,
                                     real rho, real tau) {
    int N0 = rows(share_resid);
    int Q = rows(margin_resid);
    vector[N0] z = share_resid;
    vector[N0] x = rep_vector(1.0, N0);
    vector[N0] v = rep_vector(square(sigma_share), N0);
    real slope = rho * sigma_share / sigma_margin;
    real paired_variance = square(sigma_share) * (1.0 - square(rho));
    real margin_ll = normal_lpdf(margin_resid | 0, sigma_margin);
    for (j in 1:Q) {
      int i = margin_idx[j];
      z[i] -= slope * margin_resid[j];
      v[i] = paired_variance;
    }
    return margin_ll + centered_gaussian_density(z, x, v, g, G, tau);
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
  int<lower=1, upper=2> quality_scale_parameterization;
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
  int<lower=1, upper=3> market_s0_parameterization;
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
  int<lower=0, upper=1> estimate_df;
  real<lower=2> resid_df_fixed;
  real<lower=0> prior_nu_shape;
  real<lower=0> prior_nu_rate;
  real<lower=0> prior_lkj;
  real<lower=0, upper=1> bargaining_prior_lower;
  real<lower=0, upper=1> bargaining_prior_upper;
  real<lower=0> bargaining_prior_alpha;
  real<lower=0> bargaining_prior_beta;
  int<lower=1, upper=2> structure;
  int<lower=0, upper=1> partition_prior_code;
  int<lower=0, upper=1> k_scope_code;
  int<lower=1> N_firm_mky;
  vector<lower=0>[N_firm_mky] firm_mky_shareIn;
  vector<lower=0>[N_firm_mky] firm_mky_revShare;
  array[N_firm_mky] int<lower=1> firm_mky_market;
  array[N_firm_mky] int<lower=1, upper=N_market_year> firm_mky_market_year;
  array[N_firm_mky] int<lower=1, upper=N_group> firm_mky_group;
  array[N_firm_mky] int<lower=1> firm_mky_rank;
  array[N_market_year] int<lower=1> firm_mky_start;
  array[N_market_year] int<lower=1> firm_mky_n;
  array[N] int<lower=1, upper=N_firm_mky> obs_to_firm_mky;
  int<lower=1> N_partition_unit;
  array[N_partition_unit] int<lower=1, upper=N_market> partition_unit_market;
  array[N_partition_unit] int<lower=1, upper=N_market_year> partition_unit_market_year;
  array[N_partition_unit] int<lower=1> partition_unit_start;
  array[N_partition_unit] int<lower=1> partition_unit_count;
  array[N_partition_unit] int<lower=0> partition_unit_k_min;
  array[N_partition_unit] int<lower=0> partition_unit_k_max;
  array[N_market_year] int<lower=1, upper=N_partition_unit> mky_to_partition_unit;
  int<lower=1> N_partition_component;
  array[N_partition_component] int<lower=1, upper=N_partition_unit> partition_component_unit;
  array[N_partition_component] int<lower=0> partition_component_k;
  vector[N_partition_component] partition_prior_logprob;
  vector<lower=0>[N_partition_component] partition_prior_probability;
  int<lower=1> N_firm_market;
  vector<lower=0>[N_firm_market] firm_market_shareIn;
  vector<lower=0>[N_firm_market] firm_market_revShare;
  array[N_firm_market] int<lower=1, upper=N_market> firm_market_market;
  array[N_firm_market] int<lower=1, upper=N_group> firm_market_group;
  array[N_firm_market] int<lower=1> firm_market_rank;
  array[N_market] int<lower=1> firm_market_start;
  array[N_market] int<lower=1> firm_market_n;
  array[N] int<lower=1, upper=N_firm_market> obs_to_firm_market;
  array[N] int<lower=1, upper=N> market_obs_order;
  array[N_market] int<lower=1> market_obs_start;
  array[N_market] int<lower=1> market_obs_n;
  int<lower=1, upper=2> quality_scope;
  int<lower=1> N_quality;
  array[N] int<lower=1, upper=N_quality> obs_to_quality;
  int<lower=1> N_quality_block;
  array[N_quality_block] int<lower=1, upper=N_quality> quality_block_start;
  array[N_quality_block] int<lower=1> quality_block_n;
  int<lower=0> N_quality_coord;
  array[N_quality_block] int<lower=1> quality_coord_start;
  array[N_margin] int<lower=1, upper=N_margin> market_margin_order;
  array[N_market] int<lower=1> market_margin_start;
  array[N_market] int<lower=0> market_margin_n;
  array[N_market_year] int<lower=1, upper=N_market_year> market_mky_order;
  array[N_market] int<lower=1> market_mky_start;
  array[N_market] int<lower=1> market_mky_n;
}

transformed data {
  int is_cournot = supply_model == 3;
  int is_bargaining = supply_model == 5;
  int pool_alpha_by_market = is_single_market == 0 && use_alpha_market_pool == 1;
  array[N_market] int n_mky_per_market = rep_array(0, N_market);
  array[N_market_year] int mky_pos_in_market = rep_array(0, N_market_year);
  int any_paired = 0;
  int rho_active;
  array[N] int quality_local;
  array[N] int market_obs_position;
  // Code 2 preserves the historical shifted raw rotation. Code 3 is
  // the joint orthogonal rotation, restricted to the noncentered market pool.
  int market_s0_rotation_active =
    (market_s0_parameterization == 2 && is_single_market == 0) ? 1 : 0;
  int market_s0_joint_active =
    (market_s0_parameterization == 3 && is_single_market == 0 &&
     use_s0_market_pool == 1 && use_centered_s0_pool == 0) ? 1 : 0;
  int market_s0_coordinate_active =
    (market_s0_rotation_active == 1 || market_s0_joint_active == 1) ? 1 : 0;
  for (mky in 1:N_market_year) {
    n_mky_per_market[mky_to_market[mky]] += 1;
  }
  for (m in 1:N_market)
    for (pos in 1:market_mky_n[m])
      mky_pos_in_market[market_mky_order[market_mky_start[m] + pos - 1]] = pos;
  if (structure < 1 || structure > 2)
    reject("exact ordered-k models require structure = 1 (fringe) or 2 (Stack)");
  if (supply_model != 1 && supply_model != 3)
    reject("exact ordered-k Logit model supports only Bertrand or Cournot");
  if (partition_prior_code != 1 || k_scope_code != 1)
    reject("exact ordered-k models require supplied market priors and market scope");
  if (structure == 2) for (cc in 1:N_partition_component)
    if (partition_component_k[cc] < 1)
      reject("Stackelberg ordered-k support must be strictly positive");
  for (i in 1:N) {
    if (margin_obs[i] == 1) any_paired = 1;
    quality_local[i] = quality_scope == 2
      ? obs_to_quality[i] - quality_block_start[market[i]] + 1
      : obs_to_quality[i];
  }
  for (m in 1:N_market) for (pos in 1:market_obs_n[m]) {
    int i = market_obs_order[market_obs_start[m] + pos - 1];
    market_obs_position[i] = pos;
  }
  if (market_s0_parameterization == 3 && market_s0_joint_active == 0)
    reject("joint market/s0 coordinates require multiple markets and a noncentered market-level s0 pool");
  if (quality_scale_parameterization == 2 &&
      (quality_parameterization != 2 || share_likelihood != 2 ||
       use_market_resid_scale == 1))
    reject("total/fraction quality-scale coordinates require collapsed Gaussian shares and a common residual scale");
  if (quality_scope == 1 && N_quality_block != 1)
    reject("global quality scope requires one quality block");
  if (quality_scope == 2 && N_quality_block != N_market)
    reject("market quality scope requires one quality block per market");
  if (quality_parameterization == 2 && share_likelihood == 2 &&
      margin_likelihood == 2 && use_rho == 1 && likelihood_mode == 1 &&
      any_paired == 1 && quality_scope == 1 && N_market > 1)
    reject("collapsed Gaussian rho exact-k likelihood currently requires one market");
  rho_active = (use_rho == 1 && likelihood_mode == 1 && any_paired == 1) ? 1 : 0;
}

parameters {
  vector[is_single_market == 0 ? 1 : 0] mu_b_market_vec;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;
  vector[is_single_market == 0 ? N_market : 0] b_market_raw;
  vector[quality_parameterization == 1 ? N_quality_coord : 0] b_group_centered;
  real<lower=0> sigma_b_strat;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 0 : 1] mu_log_a_free;
  vector<lower=0>[pool_alpha_by_market == 1 ? 1 : 0] sigma_log_a_vec;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
  vector[(is_cournot == 0 && pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
  vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 1 : 0] log_kappa_vec;
  vector[(is_cournot == 1 && pool_alpha_by_market == 1) ? N_market : 0] log_kappa_market;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1 ||
          market_s0_coordinate_active == 1) ? 0 : N_market_year] s0_logit;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1) ? 1 : 0] mu_s0_pool_vec;
  vector<lower=0>[(use_s0_market_year_pool == 1 ||
                  (use_s0_market_pool == 1 && is_single_market == 0)) ? 1 : 0] sigma_s0_pool_vec;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 0 &&
          market_s0_coordinate_active == 0) ? N_market_year : 0] s0_pool_raw;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 1 &&
          market_s0_coordinate_active == 0) ? N_market_year : 0] eta_s0_pool;
  // Market-level s0 pooling: one outside share per market, shared across the
  // market's years (structurally parallel to a_market). Fed by N_market group
  // effects each identified by all of the market's observations, so its hyper-SD
  // is far better identified than the market-year pool's N_market_year effects.
  vector[(use_s0_market_pool == 1 && is_single_market == 0 &&
          use_centered_s0_pool == 0 && market_s0_coordinate_active == 0) ? N_market : 0] s0_market_raw;
  vector[(use_s0_market_pool == 1 && is_single_market == 0 &&
          use_centered_s0_pool == 1 && market_s0_coordinate_active == 0) ? N_market : 0] eta_s0_market;
  vector[(market_s0_rotation_active == 1 &&
          (use_centered_s0_pool == 1 ||
           (use_s0_market_year_pool == 0 && use_s0_market_pool == 0)))
         ? (use_s0_market_pool == 1 ? N_market : N_market_year) : 0] c_s0;
  vector[(market_s0_coordinate_active == 1 && use_centered_s0_pool == 0 &&
          (use_s0_market_year_pool == 1 || use_s0_market_pool == 1))
         ? (use_s0_market_pool == 1 ? N_market : N_market_year) : 0] c_s0_raw;
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
  real share_variance_fraction = quality_scale_parameterization == 1
    ? inv_logit(2.0 * (log_sigma_share - log(sigma_b_strat)))
    : inv_logit(log_sigma_share);
  real log_sigma_share_economic = quality_scale_parameterization == 1
    ? log_sigma_share
    : log(sigma_b_strat) + 0.5 * log_inv_logit(log_sigma_share);
  real sigma_b_strat_economic = quality_scale_parameterization == 1
    ? sigma_b_strat
    : exp(log(sigma_b_strat) + 0.5 * log_inv_logit(-log_sigma_share));
  real sigma_share = exp(log_sigma_share_economic);
  real sigma_margin = exp(log_sigma_margin);
  vector[N_market] sigma_share_market;
  vector[N_market] sigma_margin_market;
  real sigma_log_a = 0.0;
  real mu_log_a;
  vector[N_market] a_market;
  vector[N_market] log_a_market_derived = rep_vector(0.0, N_market);
  vector[N_market] b_market;
  real market_s0_joint_scale = 0.0;
  vector[N_market] gamma_market;
  vector[N_market_year] s0_base_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] s0_raw_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_unbounded_mky = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_effective_mky = rep_vector(0.0, N_market_year);
  vector[quality_parameterization == 1 ? N_quality : 0] b_quality_all =
    rep_vector(0.0, quality_parameterization == 1 ? N_quality : 0);
  vector[N] share_mean = rep_vector(0.0, N);
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
  } else if (market_s0_joint_active == 1) {
    vector[N_market] u_market;
    market_s0_joint_scale = hypot(sigma_b_market_vec[1],
                                  sigma_s0_pool_vec[1]);
    u_market = (sigma_b_market_vec[1] * b_market_raw
      - sigma_s0_pool_vec[1] * c_s0_raw) / market_s0_joint_scale;
    b_market = mu_b_market_vec[1]
      + sigma_b_market_vec[1] * u_market;
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
      else if (market_s0_joint_active == 1) {
        real z = (sigma_s0_pool_vec[1] * b_market_raw[m]
          + sigma_b_market_vec[1] * c_s0_raw[m]) /
          market_s0_joint_scale;
        base_logit = mu_s0_pool_vec[1] + sigma_s0_pool_vec[1] * z;
        s0_raw_reconstructed[mky] = z;
      } else if (market_s0_rotation_active == 1) {
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
    } else base_logit = market_s0_rotation_active == 1
      ? c_s0[mky] - b_market[m] : s0_logit[mky];
    if (market_s0_coordinate_active == 0) {
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
  }

  if (is_single_market == 0) for (m in 1:N_market) {
    if (n_mky_per_market[m] > 1)
      gamma_market[m] = b_market[m] + mean_logit_s0_market[m];
    else gamma_market[m] = b_market[m];
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
  if (quality_parameterization == 1) for (block in 1:N_quality_block) {
    int G = quality_block_n[block];
    real suffix = 0.0;
    if (G > 1) {
      for (reverse_j in 1:(G - 1)) {
        int j = G + 1 - reverse_j;
        int k = j - 1;
        real v_k = b_group_centered[quality_coord_start[block] + k - 1]
          / sqrt(1.0 * k * (k + 1.0));
        b_quality_all[quality_block_start[block] + j - 1] = suffix - k * v_k;
        suffix += v_k;
      }
      b_quality_all[quality_block_start[block]] = suffix;
    }
  }
  for (i in 1:N) {
    int m = market[i];
    int mky = market_year_idx[i];
    real b_quality = (quality_parameterization == 1 ? b_quality_all[obs_to_quality[i]] : 0.0);
    share_mean[i] = logit(s0_prob[mky]) + b_market[m] + b_quality
      + market_sign * a_market[m] * price_std[i];
  }
}

model {
  real collapsed_share_ll = 0.0;
  vector[N_market] collapsed_share_ll_market = rep_vector(0.0, N_market);
  vector[N_quality] share_a = rep_vector(0.0, N_quality);
  vector[N_quality] share_c = rep_vector(0.0, N_quality);
  vector[N] share_resid = log_shareIn - share_mean;
  vector[N_market] share_within_market = rep_vector(0.0, N_market);
  vector[N_market] share_log_v_market = rep_vector(0.0, N_market);
  real share_within = 0.0;
  real share_log_v = 0.0;
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
  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);
  if (quality_scale_parameterization == 1) {
    target += normal_lpdf(sigma_b_strat | 1, prior_sigma_meanval_strat)
      - normal_lccdf(0 | 1, prior_sigma_meanval_strat);
    target += normal_lpdf(log_sigma_share | log(prior_sigma_share),
      prior_log_sigma_share_sd);
  } else {
    target += normal_lpdf(sigma_b_strat_economic | 1,
      prior_sigma_meanval_strat)
      - normal_lccdf(0 | 1, prior_sigma_meanval_strat);
    target += normal_lpdf(log_sigma_share_economic |
      log(prior_sigma_share), prior_log_sigma_share_sd);
    target += log(sigma_b_strat_economic) - log(2.0) - log(sigma_b_strat);
  }
  if (use_s0_market_pool == 1) {
    target += normal_lpdf(mu_s0_pool_vec | s0_prior_mean, s0_prior_sd);
    if (is_single_market == 0) {
      target += normal_lpdf(sigma_s0_pool_vec | 0, prior_sigma_s0_pool) - num_elements(sigma_s0_pool_vec) * normal_lccdf(0 | 0, prior_sigma_s0_pool);
      if (market_s0_joint_active == 1) {
        target += std_normal_lpdf(c_s0_raw);
      } else if (market_s0_rotation_active == 1) {
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

  // The share likelihood is common across all ordered-k states.
  if (likelihood_mode != 3) {
    if (quality_parameterization == 2 && share_likelihood == 2) {
      for (i in 1:N) {
        int q = obs_to_quality[i];
        int m = market[i];
        real z = log_shareIn[i] - share_mean[i];
        real v = square(quality_scope == 1 ? sigma_share : sigma_share_market[m]);
        if (share_a[q] > 0) {
          real residual = z - share_c[q] / share_a[q];
          real contribution = (share_a[q] / (share_a[q] + 1.0 / v))
            * square(residual) / v;
          share_within += contribution;
          share_within_market[m] += contribution;
        }
        share_a[q] += 1.0 / v;
        share_c[q] += z / v;
        share_log_v += log(v);
        share_log_v_market[m] += log(v);
      }
      if (quality_scope == 1) {
        collapsed_share_ll = centered_gaussian_from_stats(
          share_within, share_a, share_c, share_log_v, N, N_quality,
          sigma_b_strat_economic);
      } else for (m in 1:N_market) {
        int start = quality_block_start[m];
        int G = quality_block_n[m];
        vector[G] a;
        vector[G] c;
        for (j in 1:G) {
          a[j] = share_a[start + j - 1];
          c[j] = share_c[start + j - 1];
        }
        collapsed_share_ll_market[m] = centered_gaussian_from_stats(
          share_within_market[m], a, c, share_log_v_market[m],
          market_obs_n[m], G, sigma_b_strat_economic);
        collapsed_share_ll += collapsed_share_ll_market[m];
      }
      target += collapsed_share_ll;
    }
    else for (i in 1:N) {
      int m = market[i];
      if (share_likelihood == 2)
        target += normal_lpdf(log_shareIn[i] | share_mean[i], sigma_share_market[m]);
      else target += student_t_lpdf(log_shareIn[i] | nu, share_mean[i],
        sigma_share_market[m] * sqrt((nu - 2) / nu));
    }
  }

  if (likelihood_mode != 2) for (m in 1:N_market) {
    int mky0 = market_mky_order[market_mky_start[m]];
    int u = mky_to_partition_unit[mky0];
    int cstart = partition_unit_start[u];
    int nk = partition_unit_count[u];
    vector[nk] score = rep_vector(0, nk);
    for (kk in 1:nk) {
      int cc = cstart + kk - 1;
      int k = partition_component_k[cc];
      vector[market_mky_n[m]] response = rep_vector(0, market_mky_n[m]);
      score[kk] = partition_prior_probability[cc] > 0
        ? partition_prior_logprob[cc] : negative_infinity();
      // Compute Stack's follower aggregate once per component-period.
      if (structure == 2) for (q in 1:market_mky_n[m]) {
        int mky = market_mky_order[market_mky_start[m] + q - 1];
        int st = firm_mky_start[mky];
        if (supply_model == 1) {
          real B = 0;
          for (j in 1:firm_mky_n[mky]) {
            int fj = st + j - 1;
            real S = firm_mky_shareIn[fj] * (1 - s0_prob[mky]);
            if (firm_mky_rank[fj] > k)
              B += square(S) / (1 - S + square(S));
          }
          response[q] = 1 / (1 - B);
        } else {
          real S_F = 0;
          for (j in 1:firm_mky_n[mky]) {
            int fj = st + j - 1;
            if (firm_mky_rank[fj] > k)
              S_F += firm_mky_shareIn[fj] * (1 - s0_prob[mky]);
          }
          response[q] = S_F;
        }
      }
      if (quality_parameterization == 2 && share_likelihood == 2 &&
          margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1 &&
          (quality_scope == 2 || N_market == 1)) {
        if (quality_scope == 1) {
          vector[N_margin] margin_resid;
          for (q in 1:N_margin) {
            int i = margin_idx[q];
            int mky = market_year_idx[i];
            real invm = stack_logit_inv_margin(i, mky, k, s0_prob[mky],
              alpha_structural_market[m], supply_model, structure,
              response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
              firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
            margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
          }
          score[kk] += collapsed_joint_from_share_stats(
            share_within, share_a, share_c, share_log_v, share_resid,
            margin_resid, margin_idx, obs_to_quality, N_quality, sigma_share,
            sigma_margin, rho, sigma_b_strat_economic) - collapsed_share_ll;
        } else {
          int start = quality_block_start[m];
          int G = quality_block_n[m];
          int mn = market_margin_n[m];
          vector[G] a;
          vector[G] c;
          vector[mn] local_margin_resid;
          for (j in 1:G) {
            a[j] = share_a[start + j - 1];
            c[j] = share_c[start + j - 1];
          }
          for (q in 1:mn) {
            int pos = market_margin_order[market_margin_start[m] + q - 1];
            int i = margin_idx[pos];
            int mky = market_year_idx[i];
            real invm = stack_logit_inv_margin(i, mky, k, s0_prob[mky],
              alpha_structural_market[m], supply_model, structure,
              response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
              firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
            local_margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
          }
          score[kk] += collapsed_joint_from_local_share_stats(
            share_within_market[m], a, c, share_log_v_market[m],
            share_resid, local_margin_resid, market_margin_order,
            margin_idx, market_margin_start[m], market_obs_order,
          market_obs_position, market_obs_start[m], quality_local,
            market_obs_n[m], mn, G, sigma_share_market[m], sigma_margin_market[m],
            rho, sigma_b_strat_economic) - collapsed_share_ll_market[m];
        }
      } else for (pos in 1:market_obs_n[m]) {
        int i = market_obs_order[market_obs_start[m] + pos - 1];
        int mky = market_year_idx[i];
        if (margin_obs[i] == 1) {
          real invm = stack_logit_inv_margin(i, mky, k, s0_prob[mky],
            alpha_structural_market[m], supply_model, structure,
            response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
            firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
          real sharepred = share_mean[i];
          if (margin_likelihood == 2) {
            if (rho_active == 1)
              score[kk] += normal_lpdf(log_marginNormInv[i] |
                log(invm + 1e-9) + rho * (sigma_margin_market[m] / sigma_share_market[m])
                  * (log_shareIn[i] - sharepred),
                sigma_margin_market[m] * sqrt(1 - square(rho)));
            else score[kk] += normal_lpdf(log_marginNormInv[i] |
              log(invm + 1e-9), sigma_margin_market[m]);
          } else score[kk] += stack_margin_density(log_marginNormInv[i],
            log(invm + 1e-9), log_shareIn[i], sharepred,
            sigma_share_market[m], sigma_margin_market[m], rho, nu, rho_active);
        }
      }
    }
    target += log_sum_exp(score);
  }
}

generated quantities {
  vector[N] pred_shareIn = rep_vector(0, N);
  vector[N] pred_shareUncond = rep_vector(0, N);
  vector[N] pred_marginNormInv = rep_vector(0, N);
  vector[N] pred_margin = rep_vector(0, N);
  vector[N] log_lik_share_marginal = rep_vector(0, N);
  vector[N] log_lik_margin = rep_vector(0, N);
  vector[N] log_lik = rep_vector(0, N);
  vector[N_market_year] log_lik_period = rep_vector(0, N_market_year);
  vector[N_market] log_lik_market = rep_vector(0, N_market);
  vector[N_cluster] log_lik_cluster = rep_vector(0, N_cluster);
  vector[N_partition_component] posterior_component_prob = rep_vector(0, N_partition_component);
  vector[N_partition_component] prior_component_prob = rep_vector(0, N_partition_component);
  vector[N_partition_component] component_log_score = rep_vector(0, N_partition_component);
  array[N_partition_component] int component_k;
  vector[N_firm_mky] firm_leader_prob = rep_vector(0, N_firm_mky);
  vector[N_firm_mky] prior_firm_leader_prob = rep_vector(0, N_firm_mky);
  vector[N_firm_market] firm_market_core_prob = rep_vector(0, N_firm_market);
  vector[N_partition_unit] unit_margin_loglik = rep_vector(0, N_partition_unit);
  vector[N_partition_unit] posterior_k_mean = rep_vector(0, N_partition_unit);
  vector[N_partition_unit] posterior_k_median = rep_vector(0, N_partition_unit);
  vector[N_partition_unit] posterior_k_mode = rep_vector(0, N_partition_unit);
  vector[N_partition_unit] posterior_core_share_conditional = rep_vector(0, N_partition_unit);
  vector[N_market_year] posterior_core_share_conditional_mky = rep_vector(0, N_market_year);
  vector[N_market_year] posterior_core_share_unconditional_mky = rep_vector(0, N_market_year);
  vector[N_market] log_lik_share_market = rep_vector(0, N_market);
  vector[N_market] log_lik_margin_market = rep_vector(0, N_market);
  vector[N_market] log_lik_share_collapsed_market = rep_vector(0, N_market);
  vector[N_market] log_lik_joint_collapsed_market = rep_vector(0, N_market);
  int collapsed_market_lik_valid = quality_parameterization == 2
    && (quality_scope == 2 || N_market == 1)
    && likelihood_mode != 3 && share_likelihood == 2;
  real log_lik_share_collapsed = 0.0;
  real log_lik_joint_collapsed = 0.0;
  real share_ll_gq = 0.0;

  for (cc in 1:N_partition_component) component_k[cc] = partition_component_k[cc];

  {
    vector[N] sharepred;
    vector[N_quality] share_a_gq = rep_vector(0.0, N_quality);
    vector[N_quality] share_c_gq = rep_vector(0.0, N_quality);
    vector[N] share_resid_gq = rep_vector(0.0, N);
    vector[N_market] share_within_gq_market = rep_vector(0.0, N_market);
    vector[N_market] share_log_v_gq_market = rep_vector(0.0, N_market);
    real share_within_gq = 0.0;
    real share_log_v_gq = 0.0;
    vector[N_market] share_ll_gq_market = rep_vector(0.0, N_market);
    for (i in 1:N) {
      int m = market[i];
      int mky = market_year_idx[i];
      sharepred[i] = share_mean[i];
      share_resid_gq[i] = log_shareIn[i] - sharepred[i];
      if (quality_parameterization == 2 && likelihood_mode != 3) {
        int quality = obs_to_quality[i];
        real z = log_shareIn[i] - sharepred[i];
        real v = square((quality_scope == 1 || share_likelihood != 2)
          ? sigma_share : sigma_share_market[m]);
        if (share_a_gq[quality] > 0) {
          real residual = z - share_c_gq[quality] / share_a_gq[quality];
          real contribution = (share_a_gq[quality]
            / (share_a_gq[quality] + 1.0 / v))
            * square(residual) / v;
          share_within_gq += contribution;
          share_within_gq_market[m] += contribution;
        }
        share_a_gq[quality] += 1.0 / v;
        share_c_gq[quality] += z / v;
        share_log_v_gq += log(v);
        share_log_v_gq_market[m] += log(v);
      }
      pred_shareIn[i] = exp(sharepred[i]);
      pred_shareUncond[i] = (1 - s0_prob[mky]) * pred_shareIn[i];
      if (quality_parameterization != 2 && likelihood_mode != 3) {
        log_lik_share_marginal[i] = share_likelihood == 2
          ? normal_lpdf(log_shareIn[i] | sharepred[i], sigma_share_market[m])
          : student_t_lpdf(log_shareIn[i] | nu, sharepred[i],
              sigma_share_market[m] * sqrt((nu - 2) / nu));
        log_lik[i] = log_lik_share_marginal[i];
        log_lik_share_market[m] += log_lik_share_marginal[i];
      }
    }
    if (quality_parameterization == 2 && likelihood_mode != 3) {
      if (quality_scope == 1) {
        share_ll_gq = centered_gaussian_from_stats(
          share_within_gq, share_a_gq, share_c_gq, share_log_v_gq,
          N, N_quality, sigma_b_strat_economic);
      } else for (m in 1:N_market) {
        int start = quality_block_start[m];
        int G = quality_block_n[m];
        vector[G] a;
        vector[G] c;
        for (j in 1:G) {
          a[j] = share_a_gq[start + j - 1];
          c[j] = share_c_gq[start + j - 1];
        }
        share_ll_gq_market[m] = centered_gaussian_from_stats(
          share_within_gq_market[m], a, c, share_log_v_gq_market[m],
          market_obs_n[m], G, sigma_b_strat_economic);
        share_ll_gq += share_ll_gq_market[m];
      }
      if (quality_scope == 2 && share_likelihood == 2) for (m in 1:N_market) {
        log_lik_share_collapsed_market[m] = share_ll_gq_market[m];
        log_lik_share_market[m] = share_ll_gq_market[m];
      }
    }

    for (m in 1:N_market) {
      int mky0 = market_mky_order[market_mky_start[m]];
      int u = mky_to_partition_unit[mky0];
      int cstart = partition_unit_start[u];
      int nk = partition_unit_count[u];
      vector[nk] prior;
      vector[nk] priorprob;
      vector[nk] score;
      vector[nk] prefix;
      matrix[market_obs_n[m], nk] ell;
      matrix[market_obs_n[m], nk] invmargin;
      real previous;
      prior = rep_vector(0, nk);
      priorprob = rep_vector(0, nk);
      score = rep_vector(0, nk);
      ell = rep_matrix(0, market_obs_n[m], nk);
      invmargin = rep_matrix(0, market_obs_n[m], nk);

      for (kk in 1:nk) {
        int cc = cstart + kk - 1;
        int k = partition_component_k[cc];
        vector[market_mky_n[m]] response;
        response = rep_vector(0, market_mky_n[m]);
        prior[kk] = partition_prior_probability[cc] > 0
          ? partition_prior_logprob[cc] : negative_infinity();
        priorprob[kk] = partition_prior_probability[cc];
        score[kk] = partition_prior_probability[cc] > 0
          ? partition_prior_logprob[cc] : negative_infinity();
        if (structure == 2) for (q in 1:market_mky_n[m]) {
          int mky = market_mky_order[market_mky_start[m] + q - 1];
          int st = firm_mky_start[mky];
          if (supply_model == 1) {
            real B = 0;
            for (j in 1:firm_mky_n[mky]) {
              int fj = st + j - 1;
              real S = firm_mky_shareIn[fj] * (1 - s0_prob[mky]);
              if (firm_mky_rank[fj] > k)
                B += square(S) / (1 - S + square(S));
            }
            response[q] = 1 / (1 - B);
          } else {
            real S_F = 0;
            for (j in 1:firm_mky_n[mky]) {
              int fj = st + j - 1;
              if (firm_mky_rank[fj] > k)
                S_F += firm_mky_shareIn[fj] * (1 - s0_prob[mky]);
            }
            response[q] = S_F;
          }
        }
        if (quality_parameterization == 2 && share_likelihood == 2 &&
            margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1 &&
            (quality_scope == 2 || N_market == 1)) {
          if (quality_scope == 1) {
            vector[N_margin] margin_resid;
            for (q in 1:N_margin) {
              int i = margin_idx[q];
              int mm = market[i];
              int mky = market_year_idx[i];
              real invm = stack_logit_inv_margin(i, mky, k, s0_prob[mky],
                alpha_structural_market[mm], supply_model, structure,
                response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
                firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
              margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
            }
            score[kk] += collapsed_joint_from_share_stats(
              share_within_gq, share_a_gq, share_c_gq, share_log_v_gq,
              share_resid_gq, margin_resid, margin_idx, obs_to_quality,
              N_quality, sigma_share, sigma_margin, rho, sigma_b_strat_economic)
              - share_ll_gq;
          } else {
            int qstart = quality_block_start[m];
            int G = quality_block_n[m];
            int mn = market_margin_n[m];
            vector[G] a;
            vector[G] c;
            vector[mn] local_margin_resid;
            for (j in 1:G) {
              a[j] = share_a_gq[qstart + j - 1];
              c[j] = share_c_gq[qstart + j - 1];
            }
            for (q in 1:mn) {
              int pos = market_margin_order[market_margin_start[m] + q - 1];
              int i = margin_idx[pos];
              int mky = market_year_idx[i];
              real invm = stack_logit_inv_margin(i, mky, k, s0_prob[mky],
                alpha_structural_market[m], supply_model, structure,
                response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
                firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
              local_margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
            }
            score[kk] += collapsed_joint_from_local_share_stats(
              share_within_gq_market[m], a, c, share_log_v_gq_market[m],
              share_resid_gq, local_margin_resid, market_margin_order,
              margin_idx, market_margin_start[m], market_obs_order,
          market_obs_position, market_obs_start[m], quality_local,
              market_obs_n[m], mn, G, sigma_share_market[m], sigma_margin_market[m],
              rho, sigma_b_strat_economic) - share_ll_gq_market[m];
          }
        }
        for (pos in 1:market_obs_n[m]) {
          int i = market_obs_order[market_obs_start[m] + pos - 1];
          int mky = market_year_idx[i];
          invmargin[pos, kk] = stack_logit_inv_margin(i, mky, k,
            s0_prob[mky], alpha_structural_market[m], supply_model, structure,
            response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
            firm_mky_shareIn, firm_mky_rank, obs_to_firm_mky, price);
          if (likelihood_mode != 2 && margin_obs[i] == 1 &&
              !(quality_parameterization == 2 && share_likelihood == 2 &&
                margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1 &&
                (quality_scope == 2 || N_market == 1))) {
            if (margin_likelihood == 2) {
              ell[pos, kk] = rho_active == 1
                ? normal_lpdf(log_marginNormInv[i] |
                    log(invmargin[pos, kk] + 1e-9)
                      + rho * (sigma_margin_market[m] / sigma_share_market[m])
                        * (log_shareIn[i] - sharepred[i]),
                    sigma_margin_market[m] * sqrt(1 - square(rho)))
                : normal_lpdf(log_marginNormInv[i] |
                    log(invmargin[pos, kk] + 1e-9), sigma_margin_market[m]);
            } else ell[pos, kk] = stack_margin_density(log_marginNormInv[i],
              log(invmargin[pos, kk] + 1e-9), log_shareIn[i], sharepred[i],
              sigma_share_market[m], sigma_margin_market[m], rho, nu, rho_active);
            score[kk] += ell[pos, kk];
          }
        }
        component_log_score[cc] = score[kk];
      }

      {
        vector[nk] postprob = softmax(score);
        real mode_prob = -1;
        for (kk in 1:nk) {
          int cc = cstart + kk - 1;
          posterior_component_prob[cc] = postprob[kk];
          prior_component_prob[cc] = priorprob[kk];
          posterior_k_mean[u] += postprob[kk] * partition_component_k[cc];
          if (postprob[kk] > mode_prob) {
            mode_prob = postprob[kk];
            posterior_k_mode[u] = partition_component_k[cc];
          }
        }
        {
          int median_k = partition_component_k[cstart];
          // k_support is positional and need not be sorted.  Evaluate each
          // numeric support point's CDF and choose the smallest weighted
          // median, rather than treating component order as numeric order.
          for (kk in 1:nk) {
            int cc = cstart + kk - 1;
            if (partition_component_k[cc] > median_k)
              median_k = partition_component_k[cc];
          }
          for (kk in 1:nk) {
            int cc = cstart + kk - 1;
            real cum = 0;
            for (rr in 1:nk) {
              int cr = cstart + rr - 1;
              if (partition_component_k[cr] <= partition_component_k[cc])
                cum += postprob[rr];
            }
            if (cum >= 0.5 && partition_component_k[cc] < median_k)
              median_k = partition_component_k[cc];
          }
          posterior_k_median[u] = median_k;
        }
        unit_margin_loglik[u] = likelihood_mode == 2 ? 0 : log_sum_exp(score);
        log_lik_margin_market[m] = unit_margin_loglik[u];
        for (pos in 1:market_obs_n[m]) {
          int i = market_obs_order[market_obs_start[m] + pos - 1];
          for (kk in 1:nk) {
            pred_marginNormInv[i] += postprob[kk] * invmargin[pos, kk];
            pred_margin[i] += postprob[kk] / invmargin[pos, kk];
          }
        }
        // Ordered accounting factors telescope to the single market mixture.
        prefix = prior;
        previous = log_sum_exp(prefix);
        if (likelihood_mode != 2) for (pos in 1:market_obs_n[m]) {
          int i = market_obs_order[market_obs_start[m] + pos - 1];
          if (margin_obs[i] == 1) {
            for (kk in 1:nk) prefix[kk] += ell[pos, kk];
            {
              real current = log_sum_exp(prefix);
              log_lik_margin[i] = current - previous;
              log_lik[i] += log_lik_margin[i];
              previous = current;
            }
          }
        }
      }
      log_lik_market[m] = log_lik_share_market[m] + log_lik_margin_market[m];
    }
  }

  // Fixed ranking makes these role/core probabilities invariant to row order.
  for (fm in 1:N_firm_mky) {
    int u = mky_to_partition_unit[firm_mky_market_year[fm]];
    int cstart = partition_unit_start[u];
    for (kk in 1:partition_unit_count[u]) {
      int cc = cstart + kk - 1;
      if (firm_mky_rank[fm] <= partition_component_k[cc]) {
        firm_leader_prob[fm] += posterior_component_prob[cc];
        prior_firm_leader_prob[fm] += prior_component_prob[cc];
      }
    }
  }
  for (fm in 1:N_firm_market) {
    int m = firm_market_market[fm];
    int u = mky_to_partition_unit[market_mky_order[market_mky_start[m]]];
    int cstart = partition_unit_start[u];
    for (kk in 1:partition_unit_count[u]) {
      int cc = cstart + kk - 1;
      if (firm_market_rank[fm] <= partition_component_k[cc])
        firm_market_core_prob[fm] += posterior_component_prob[cc];
    }
  }
  for (mky in 1:N_market_year) {
    int u = mky_to_partition_unit[mky];
    int cstart = partition_unit_start[u];
    int st = firm_mky_start[mky];
    for (kk in 1:partition_unit_count[u]) {
      int cc = cstart + kk - 1;
      real core = 0;
      for (j in 1:firm_mky_n[mky]) {
        int fj = st + j - 1;
        if (firm_mky_rank[fj] <= partition_component_k[cc])
          core += firm_mky_shareIn[fj];
      }
      posterior_core_share_conditional_mky[mky] += posterior_component_prob[cc] * core;
      posterior_core_share_unconditional_mky[mky] += posterior_component_prob[cc]
        * core * (1 - s0_prob[mky]);
    }
  }
  for (m in 1:N_market) {
    int u = mky_to_partition_unit[market_mky_order[market_mky_start[m]]];
    int cstart = partition_unit_start[u];
    for (kk in 1:partition_unit_count[u]) {
      int cc = cstart + kk - 1;
      real core = 0;
      for (j in 1:firm_market_n[m]) {
        int fj = firm_market_start[m] + j - 1;
        if (firm_market_rank[fj] <= partition_component_k[cc])
          core += firm_market_shareIn[fj];
      }
      posterior_core_share_conditional[u] += posterior_component_prob[cc] * core;
    }
  }
  for (i in 1:N) {
    log_lik_period[market_year_idx[i]] += log_lik[i];
    log_lik_cluster[cluster[i]] += log_lik[i];
  }
  if (quality_parameterization == 2) {
    log_lik_share_collapsed = share_ll_gq;
    log_lik_joint_collapsed = share_ll_gq;
    if (quality_scope == 2 && share_likelihood == 2) for (m in 1:N_market) {
      int mky0 = market_mky_order[market_mky_start[m]];
      int u = mky_to_partition_unit[mky0];
      log_lik_joint_collapsed_market[m] = log_lik_share_collapsed_market[m]
        + unit_margin_loglik[u];
      log_lik_joint_collapsed += unit_margin_loglik[u];
    } else {
      if (N_market == 1 && likelihood_mode != 3 && share_likelihood == 2) {
        int mky0 = market_mky_order[market_mky_start[1]];
        int u = mky_to_partition_unit[mky0];
        log_lik_share_collapsed_market[1] = share_ll_gq;
        log_lik_joint_collapsed_market[1] = share_ll_gq
          + (likelihood_mode != 2 ? unit_margin_loglik[u] : 0);
      }
      if (likelihood_mode != 2)
        for (u in 1:N_partition_unit)
          log_lik_joint_collapsed += unit_margin_loglik[u];
    }
  }
}
