// Specification 9: centered orthonormal firm-quality coordinates.
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

  // Return the inverse Lerner (the observed CES margin response is
  // log(1 / Lerner)).  response is phi for Stackelberg Bertrand and C_F for
  // Stackelberg Cournot.  Fringe states have exact product-level behavior.
  real stack_ces_inv_margin(int i, int mky, int k, real s0,
                            real outside_price, real inside_price_index,
                            real gamma, int conduct, int structure,
                            real response, array[] int firm_mky_start,
                            array[] int firm_mky_n, vector firm_mky_revShare,
                            array[] int firm_mky_rank,
                            array[] int obs_to_firm_mky) {
    int me = obs_to_firm_mky[i];
    int st = firm_mky_start[mky];
    int nf = firm_mky_n[mky];
    real R0 = inv_logit(logit(s0) + log(outside_price)
      - log(inside_price_index));
    real R = firm_mky_revShare[me] * (1 - R0);
    real h = gamma - 1;
    real H = 1 + h * R0;
    real b = h / (gamma * H);
    real m;
    if (structure == 1) {
      if (firm_mky_rank[me] > k) m = 1 / gamma;
      else if (conduct == 1) m = 1 / (gamma - h * R);
      else m = 1 / gamma + b * R;
    } else if (conduct == 1) {
      m = firm_mky_rank[me] > k ? 1 / (gamma - h * R)
        : 1 / (gamma - h * response * R);
    } else {
      m = firm_mky_rank[me] > k ? 1 / gamma + b * R
        : 1 / gamma + b * R / (1 + response);
    }
    return 1 / m;
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
    vector[G] a = share_a;
    vector[G] c = share_c;
    real within = share_within;
    real logv = share_log_v;
    real slope = rho * sigma_margin / sigma_share;
    real v_margin = square(sigma_margin) * (1.0 - square(rho));
    for (j in 1:Q) {
      int i = margin_idx[j];
      real z = margin_resid[j] - slope * share_resid[i];
      real weight = square(slope) / v_margin;
      real residual = z + slope * c[g[i]] / a[g[i]];
      within += (a[g[i]] / (a[g[i]] + weight)) * square(residual) / v_margin;
      a[g[i]] += weight;
      c[g[i]] += -slope * z / v_margin;
      logv += log(v_margin);
    }
    return centered_gaussian_from_stats(within, a, c, logv, N0 + Q, G, tau);
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
  array[N_market_year] int<lower=1, upper=N_market_year> market_mky_order;
  array[N_market] int<lower=1> market_mky_start;
  array[N_market] int<lower=1> market_mky_n;
}

transformed data {
  int is_bargaining = supply_model == 5;
  int pool_alpha_by_market = is_single_market == 0 && use_alpha_market_pool == 1;
  array[N_market] int n_mky_per_market = rep_array(0, N_market);
  array[N_market_year] int mky_pos_in_market = rep_array(0, N_market_year);
  int any_paired = 0;
  int rho_active;
  int market_s0_rotation_active =
    (market_s0_parameterization == 2 && is_single_market == 0) ? 1 : 0;
  for (mky in 1:N_market_year) {
    n_mky_per_market[mky_to_market[mky]] += 1;
  }
  for (m in 1:N_market)
    for (pos in 1:market_mky_n[m])
      mky_pos_in_market[market_mky_order[market_mky_start[m] + pos - 1]] = pos;
  if (structure < 1 || structure > 2)
    reject("exact ordered-k models require structure = 1 (fringe) or 2 (Stack)");
  if (supply_model != 1 && supply_model != 3)
    reject("exact ordered-k CES model supports only Bertrand or Cournot");
  if (partition_prior_code != 1 || k_scope_code != 1)
    reject("exact ordered-k models require supplied market priors and market scope");
  if (structure == 2) for (cc in 1:N_partition_component)
    if (partition_component_k[cc] < 1)
      reject("Stackelberg ordered-k support must be strictly positive");
  for (i in 1:N) {
    if (margin_obs[i] == 1) any_paired = 1;
  }
  if (quality_parameterization == 2 && share_likelihood == 2 &&
      margin_likelihood == 2 && use_rho == 1 && likelihood_mode == 1 &&
      any_paired == 1 && N_market > 1)
    reject("collapsed Gaussian rho exact-k likelihood currently requires one market");
  rho_active = (use_rho == 1 && likelihood_mode == 1 && any_paired == 1) ? 1 : 0;
}

parameters {
  vector[is_single_market == 0 ? 1 : 0] mu_b_market_vec;
  vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;
  vector[is_single_market == 0 ? N_market : 0] b_market_raw;
  vector[quality_parameterization == 1 ? N_group - 1 : 0] b_group_centered;
  real<lower=0> sigma_b_strat;
  real mu_log_a;
  vector<lower=0>[pool_alpha_by_market == 1 ? 1 : 0] sigma_log_a_vec;
  vector[(pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
  vector[(pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1 ||
          market_s0_rotation_active == 1)
         ? 0 : N_market_year] s0_logit;
  vector[(use_s0_market_year_pool == 1 || use_s0_market_pool == 1)
         ? 1 : 0] mu_s0_pool_vec;
  vector<lower=0>[(use_s0_market_year_pool == 1 ||
                  (use_s0_market_pool == 1 && is_single_market == 0))
                  ? 1 : 0] sigma_s0_pool_vec;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 0 &&
          market_s0_rotation_active == 0) ? N_market_year : 0] s0_pool_raw;
  vector[(use_s0_market_year_pool == 1 && use_centered_s0_pool == 1 &&
          market_s0_rotation_active == 0) ? N_market_year : 0] eta_s0_pool;
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
  vector[N_market] b_market_rot;
  vector[N_market_year] s0_base_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] s0_raw_reconstructed = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_unbounded_mky = rep_vector(0.0, N_market_year);
  vector[N_market_year] c_effective_mky = rep_vector(0.0, N_market_year);
  vector[quality_parameterization == 1 ? N_group : 0] b_quality_all =
    rep_vector(0.0, quality_parameterization == 1 ? N_group : 0);
  vector[N] share_mean = rep_vector(0.0, N);
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
  if (is_bargaining == 1) {
    bargaining_power = bargaining_prior_lower
      + (bargaining_prior_upper - bargaining_prior_lower) * bargaining_unit[1];
  }
  if (rho_active == 1) rho = 0.95 * tanh(rho_raw[1]);

  // Reconstruct the unchanged market effect before applying the outside-share
  // rotation. Single-market fits retain the historical zero normalization.
  if (is_single_market == 1) {
    b_market = rep_vector(0.0, N_market);
    b_market_rot = rep_vector(0.0, N_market);
  } else {
    b_market = mu_b_market_vec[1]
      + sigma_b_market_vec[1] * b_market_raw;
  }

  for (mky in 1:N_market_year) {
    real base_logit;
    real x_beta = K_s0 > 0 ? X_s0[mky] * beta_s0 : 0.0;
    int m = mky_to_market[mky];
    if (use_s0_market_pool == 1) {
      // One quantity outside share per market, applied to every market-year.
      // The CES revenue share R0 is transformed below using period-specific
      // outside and inside price indices, so R0 may still vary by year.
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
    } else base_logit = market_s0_rotation_active == 1
      ? c_s0[mky] - b_market[m] : s0_logit[mky];
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
    // Physical outside revenue share; no second epsilon transform.
    R0_prob[mky] = inv_logit(logit(s0_prob[mky])
      + log(outside_price[mky]) - log(inside_price_index[mky]));
    if (R0_prob[mky] <= 0 || R0_prob[mky] >= 1)
      reject("CES outside revenue share must be strictly between zero and one");
    c_effective_mky[mky] = logit(s0_prob[mky]) + b_market[m];
    mean_logit_s0_market[mky_to_market[mky]] += logit(s0_prob[mky]);
  }
  for (m in 1:N_market) {
    mean_logit_s0_market[m] /= n_mky_per_market[m];
    if (is_single_market == 0) {
      if (n_mky_per_market[m] > 1) {
        b_market_rot[m] = b_market[m] + mean_logit_s0_market[m];
      } else b_market_rot[m] = b_market[m];
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
    real b_quality = (quality_parameterization == 1 ? b_quality_all[group[i]] : 0.0);
    share_mean[i] = logit(s0_prob[mky]) + b_market[m] + b_quality
      + market_sign * (a_market[m] + log_price_sd) * log_price_std[i]
      + gamma_structural_market[m] * log(outside_price[mky]);
  }
}

model {
  real collapsed_share_ll = 0.0;
  vector[N_group] share_a = rep_vector(0.0, N_group);
  vector[N_group] share_c = rep_vector(0.0, N_group);
  real share_within = 0.0;
  real share_log_v = 0.0;
  target += normal_lpdf(mu_log_a | prior_alpha_mean, prior_alpha_sd);
  if (pool_alpha_by_market == 1) {
    target += normal_lpdf(sigma_log_a_vec | 0, prior_sigma_alpha) - num_elements(sigma_log_a_vec) * normal_lccdf(0 | 0, prior_sigma_alpha);
    if (use_centered_alpha == 1) target += normal_lpdf(log_a_market_pool | mu_log_a, sigma_log_a);
    else target += std_normal_lpdf(r_market_a_raw);
  }
  if (is_single_market == 0) {
    target += normal_lpdf(mu_b_market_vec | 0, 1);
    target += normal_lpdf(sigma_b_market_vec | 0, 0.5) - num_elements(sigma_b_market_vec) * normal_lccdf(0 | 0, 0.5);
    target += std_normal_lpdf(b_market_raw);
  }
  target += normal_lpdf(b_group_centered | 0, sigma_b_strat);
  target += normal_lpdf(sigma_b_strat | 1, prior_sigma_meanval_strat) - 1 * normal_lccdf(0 | 1, prior_sigma_meanval_strat);
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
    target += (prior_lkj - 1.0) * log1m(square(rho))
      + log1m(square(tanh(rho_raw[1])))
      - (lbeta(0.5, prior_lkj) + beta_lcdf(square(0.95) | 0.5, prior_lkj) - log(0.95));
  }

  if (use_hmt == 1 && likelihood_mode != 2) {
    real gamma_check = 1.0 + exp(mu_log_a) / log_price_sd;
    real R0_check = mean(R0_prob);
    real agg_elasticity = 1.0 + (gamma_check - 1.0) * R0_check;
    real max_elasticity = 1.0 / (ssnip_hmt + avg_margin_hmt);
    target += soft_hmt_penalty(agg_elasticity, max_elasticity, 0.1);
  }

  // The share likelihood is common across all ordered-k states.
  if (likelihood_mode != 3) {
    if (quality_parameterization == 2 && share_likelihood == 2) {
      for (i in 1:N) {
        real z = log_shareIn[i] - share_mean[i];
        real v = square(sigma_share);
        if (share_a[group[i]] > 0) {
          real residual = z - share_c[group[i]] / share_a[group[i]];
          share_within += (share_a[group[i]] / (share_a[group[i]] + 1.0 / v))
            * square(residual) / v;
        }
        share_a[group[i]] += 1.0 / v;
        share_c[group[i]] += z / v;
        share_log_v += log(v);
      }
      collapsed_share_ll = centered_gaussian_from_stats(
        share_within, share_a, share_c, share_log_v, N, N_group, sigma_b_strat);
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
        real R0 = R0_prob[mky];
        if (supply_model == 1) {
          real C_F = 0;
          for (j in 1:firm_mky_n[mky]) {
            int fj = st + j - 1;
            real R = firm_mky_revShare[fj] * (1 - R0);
            real h = gamma_structural_market[m] - 1;
            real t = gamma_structural_market[m] - h * R;
            real af = square(h) * R / (t * (t - 1) + square(h) * R);
            if (firm_mky_rank[fj] > k) C_F += R * af;
          }
          response[q] = 1 / (1 - C_F);
        } else {
          real H = 1 + (gamma_structural_market[m] - 1) * R0;
          real h = gamma_structural_market[m] - 1;
          real C_F = 0;
          for (j in 1:firm_mky_n[mky]) {
            int fj = st + j - 1;
            real R = firm_mky_revShare[fj] * (1 - R0);
            real af = h * (H - 2 * R + gamma_structural_market[m] * h * R0 * R / H)
              / (H * (H + (gamma_structural_market[m] - 2) * R));
            if (firm_mky_rank[fj] > k) C_F += R * af;
          }
          response[q] = C_F;
        }
      }
      if (quality_parameterization == 2 && share_likelihood == 2 &&
          margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1 &&
          N_market == 1) {
        vector[N_margin] margin_resid;
        for (q in 1:N_margin) {
          int i = margin_idx[q];
          int mky = market_year_idx[i];
          real invm = stack_ces_inv_margin(i, mky, k, s0_prob[mky],
            outside_price[mky], inside_price_index[mky], gamma_structural_market[market[i]],
            supply_model, structure, response[mky_pos_in_market[mky]],
            firm_mky_start, firm_mky_n, firm_mky_revShare, firm_mky_rank,
            obs_to_firm_mky);
          margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
        }
        score[kk] += collapsed_joint_from_share_stats(
          share_within, share_a, share_c, share_log_v, log_shareIn - share_mean,
          margin_resid, margin_idx, group, N_group, sigma_share, sigma_margin,
          rho, sigma_b_strat) - collapsed_share_ll;
      } else for (pos in 1:market_obs_n[m]) {
        int i = market_obs_order[market_obs_start[m] + pos - 1];
        int mky = market_year_idx[i];
        if (margin_obs[i] == 1) {
          real invm = stack_ces_inv_margin(i, mky, k, s0_prob[mky],
            outside_price[mky], inside_price_index[mky],
            gamma_structural_market[m], supply_model, structure,
            response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
            firm_mky_revShare, firm_mky_rank, obs_to_firm_mky);
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
  real log_lik_share_collapsed = 0.0;
  real log_lik_joint_collapsed = 0.0;
  real share_ll_gq = 0.0;

  for (cc in 1:N_partition_component) component_k[cc] = partition_component_k[cc];

  {
    vector[N] sharepred;
    vector[N_group] share_a_gq = rep_vector(0.0, N_group);
    vector[N_group] share_c_gq = rep_vector(0.0, N_group);
    real share_within_gq = 0.0;
    real share_log_v_gq = 0.0;
    for (i in 1:N) {
      int m = market[i];
      int mky = market_year_idx[i];
      sharepred[i] = share_mean[i];
      if (quality_parameterization == 2 && likelihood_mode != 3) {
        real z = log_shareIn[i] - sharepred[i];
        real v = square(sigma_share);
        if (share_a_gq[group[i]] > 0) {
          real residual = z - share_c_gq[group[i]] / share_a_gq[group[i]];
          share_within_gq += (share_a_gq[group[i]] / (share_a_gq[group[i]] + 1.0 / v))
            * square(residual) / v;
        }
        share_a_gq[group[i]] += 1.0 / v;
        share_c_gq[group[i]] += z / v;
        share_log_v_gq += log(v);
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
    if (quality_parameterization == 2 && likelihood_mode != 3)
      share_ll_gq = centered_gaussian_from_stats(
        share_within_gq, share_a_gq, share_c_gq, share_log_v_gq,
        N, N_group, sigma_b_strat);

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
          real R0 = R0_prob[mky];
          if (supply_model == 1) {
            real C_F = 0;
            for (j in 1:firm_mky_n[mky]) {
              int fj = st + j - 1;
              real R = firm_mky_revShare[fj] * (1 - R0);
              real h = gamma_structural_market[m] - 1;
              real t = gamma_structural_market[m] - h * R;
              real af = square(h) * R / (t * (t - 1) + square(h) * R);
              if (firm_mky_rank[fj] > k) C_F += R * af;
            }
            response[q] = 1 / (1 - C_F);
          } else {
            real H = 1 + (gamma_structural_market[m] - 1) * R0;
            real h = gamma_structural_market[m] - 1;
            real C_F = 0;
            for (j in 1:firm_mky_n[mky]) {
              int fj = st + j - 1;
              real R = firm_mky_revShare[fj] * (1 - R0);
              real af = h * (H - 2 * R + gamma_structural_market[m] * h * R0 * R / H)
                / (H * (H + (gamma_structural_market[m] - 2) * R));
              if (firm_mky_rank[fj] > k) C_F += R * af;
            }
            response[q] = C_F;
          }
        }
        if (quality_parameterization == 2 && share_likelihood == 2 &&
            margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1) {
          vector[N_margin] margin_resid;
          for (q in 1:N_margin) {
            int i = margin_idx[q];
            int mm = market[i];
            int mky = market_year_idx[i];
            real invm = stack_ces_inv_margin(i, mky, k, s0_prob[mky],
              outside_price[mky], inside_price_index[mky], gamma_structural_market[mm],
              supply_model, structure, response[mky_pos_in_market[mky]],
              firm_mky_start, firm_mky_n, firm_mky_revShare, firm_mky_rank,
              obs_to_firm_mky);
            margin_resid[q] = log_marginNormInv[i] - log(invm + 1e-9);
          }
          score[kk] += collapsed_joint_from_share_stats(
            share_within_gq, share_a_gq, share_c_gq, share_log_v_gq,
            log_shareIn - sharepred, margin_resid, margin_idx, group, N_group,
            sigma_share, sigma_margin, rho, sigma_b_strat) - share_ll_gq;
        }
        for (pos in 1:market_obs_n[m]) {
          int i = market_obs_order[market_obs_start[m] + pos - 1];
          int mky = market_year_idx[i];
          invmargin[pos, kk] = stack_ces_inv_margin(i, mky, k, s0_prob[mky],
            outside_price[mky], inside_price_index[mky],
            gamma_structural_market[m], supply_model, structure,
            response[mky_pos_in_market[mky]], firm_mky_start, firm_mky_n,
            firm_mky_revShare, firm_mky_rank, obs_to_firm_mky);
          if (likelihood_mode != 2 && margin_obs[i] == 1 &&
              !(quality_parameterization == 2 && share_likelihood == 2 &&
                margin_likelihood == 2 && likelihood_mode != 3 && rho_active == 1 &&
                N_market == 1)) {
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
          core += firm_mky_revShare[fj];
      }
      posterior_core_share_conditional_mky[mky] += posterior_component_prob[cc] * core;
      posterior_core_share_unconditional_mky[mky] += posterior_component_prob[cc]
        * core * (1 - R0_prob[mky]);
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
          core += firm_market_revShare[fj];
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
    if (likelihood_mode != 2)
      for (u in 1:N_partition_unit)
        log_lik_joint_collapsed += unit_margin_loglik[u];
  }
}
