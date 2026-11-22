  // 0 = irregular incidence, 1 = rectangular shares, 2 = conditional
  // variances are also constant within each period (including masked periods).
  int panel_design(array[] int group, int G, array[] int cell, int T,
      array[] int observed, int rho_active) {
    int L = size(group);
    if (L % T != 0 || L %/% T != G) return 0;
    {
      array[G, T] int seen = rep_array(0, G, T);
      array[T] int paired = rep_array(0, T);
      for (i in 1:L) {
        if (seen[group[i], cell[i]]) return 0;
        seen[group[i], cell[i]] = 1;
        paired[cell[i]] += observed[i];
      }
      if (rho_active) for (t in 1:T)
        if (paired[t] > 0 && paired[t] < G) return 1;
    }
    return 2;
  }

  // Each orthonormal within-period product contrast has time covariance
  // diag(v) + tau^2 * 11'. Center in product space, then integrate its time mean.
  real panel_gaussian_logpdf(int eligible, vector y, vector variance,
      array[] int group, int G, array[] int cell, int T, real tau) {
    if (!eligible) return relative_gaussian_logpdf(y, variance, group, G, cell, T, tau);
    if (G == 1) return 0;
    {
      matrix[G, T] residual;
      vector[T] v;
      vector[T] precision;
      vector[G] average;
      real total_precision;
      real quadratic;
      for (i in 1:num_elements(y)) {
        residual[group[i], cell[i]] = y[i];
        v[cell[i]] = variance[i];
      }
      for (t in 1:T) residual[, t] -= mean(residual[, t]);
      precision = inv(v);
      total_precision = sum(precision);
      average = residual * (precision / total_precision);
      quadratic = dot_self(average) / (inv(total_precision) + square(tau));
      for (t in 1:T)
        quadratic += precision[t] * dot_self(residual[, t] - average);
      return -0.5 * (quadratic + (G - 1) * (T * log(2 * pi())
        + sum(log(v)) + log1p(square(tau) * total_precision)));
    }
  }

  real panel_local_gaussian_logpdf(int eligible, vector y, vector variance,
      array[] int obs_order, int start, int L, array[] int group,
      int G, array[] int cell, int T, real tau) {
    array[L] int rows = obs_order[start:(start + L - 1)];
    return panel_gaussian_logpdf(eligible, y[rows], variance[rows],
      group[rows], G, cell[rows], T, tau);
  }
