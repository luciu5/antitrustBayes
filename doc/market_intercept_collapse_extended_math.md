# Extended market-intercept collapse oracle

This is a reference-only extension to `market_intercept_collapse_reference.R`.
It uses dense Gaussian matrices only for small fixtures and does not implement
or change the production model. Run it with:

```sh
Rscript --vanilla doc/market_intercept_collapse_extended_reference.R
```

The script writes three compact CSVs beside itself: likelihood comparisons,
joint-posterior comparisons, and exact-k mixture comparisons.

## Conditional margins and normalized likelihood

Write the intercept as `b = mu_b + a`, with `a ~ Normal(0, sigma_b^2)`. For
each market, let `r_s = y_s - eta_s(k) - mu_b` and
`r_m = y_m - eta_m(k)` for observed margins only. If row `i` has share scale
`s_i`, observed margin `j` is paired with share row `i(j)`, its margin scale
is `t_j`, and their residual correlation is `rho`, conditioning on the
margins gives

```text
z_i = r_s,i - rho * s_i / t_j * r_m,j    (paired row i = i(j))
v_i = s_i^2 * (1 - rho^2)                (paired row)
v_i = s_i^2                              (missing-margin row)
```

The observed margin contribution is the normalized sum
`sum_j log Normal(r_m,j | 0, t_j^2)`. Conditional on those margins,
`z = 1*a + Z*q + epsilon`, where `epsilon ~ Normal(0, diag(v))`,
`q ~ Normal(0, tau^2 P)`, and `P = I - 11'/G`. Integrating `a` and `q`
therefore gives a normalized Gaussian with covariance

```text
V = diag(v) + tau^2 Z P Z' + sigma_b^2 11'.
```

The reference compares the structured density with the direct dense density
of the jointly observed shares and nonmissing margins. Both include all
Gaussian normalization constants. The scales differ by market in the exact-k
fixture; the partial-margin fixture also exercises conditional heteroskedastic
scales through its margin mask.

For a linear-cost calculation, group `g` contributes
`A_g = sum(i in g) 1/v_i` and `C_g = sum(i in g) z_i/v_i`. Define
`m_g = C_g/A_g`, `h_g = 1/(1 + tau^2 A_g)`,
`A = sum_g A_g h_g`, `m = sum_g A_g h_g m_g / A`, and `d = mean_g(h_g)`.
The centered intercept posterior is

```text
E[a | data] = sigma_b^2 A m / (d + sigma_b^2 A)
Var(a | data) = sigma_b^2 d / (d + sigma_b^2 A).
```

The script retains the within-group residual sum, contrast among group means,
determinant terms, and intercept direction in the normalized density.

## Joint recovery and exact-k factorization

Let `U` be an orthonormal basis for the centered product subspace, and write
`q = U xi`. With `H = [1, ZU]`, `x = (a, xi)`, prior covariance
`Lambda = diag(sigma_b^2, tau^2 I)`, and `D = diag(v)`, the posterior is

```text
S = (Lambda^-1 + H' D^-1 H)^-1
m = S H' D^-1 z.
```

Mapping `(a, xi)` to `(b, q)`, adding `mu_b` to the intercept mean, yields all
means, variances, and intercept-quality cross-covariances. The script checks
this precision form against the independent Gaussian conditioning form
`Lambda - Lambda H' (H Lambda H' + D)^-1 H Lambda`.

For coherent recovery, it draws `b` from its marginal posterior and then draws
`q` from its Gaussian conditional at that realized `b`. The conditional
moments are computed directly in contrast coordinates. A seeded simulation
checks all joint means and covariance entries, while a separate algebraic
check compares `q | b` with the conditional moments obtained from the joint
posterior covariance.

For market `m` and state `k`, let `ell_mk` be its fully normalized integrated
log likelihood and `pi_mk` its prior mass. The factorized exact-k target is

```text
sum_m log_sum_exp_k(log(pi_mk) + ell_mk).
```

The small oracle explicitly enumerates all `K^M = 3^3 = 27` market-state
assignments, compares the total log evidence, and recovers each market's state
probabilities from the enumeration. The same market state persists across its
fixture's multiple years.

## Results

All extension checks passed: 3/3 dense likelihood fixtures, 3/3 joint posterior
and sequential-recovery fixtures, and 1/1 exact-k factorization fixture. The
maximum dense-versus-structured likelihood error was `7.11e-15`; maximum
joint-posterior mean and covariance errors were `5.07e-16` and `2.09e-17`.
The maximum `q | b` formula discrepancy was `6.18e-16`. Across 30,000
sequential draws per fixture, the largest mean and covariance Monte Carlo
standardized errors were 1.40 and 2.32, below the prespecified threshold of 6.
Exact enumeration had zero log-evidence discrepancy and a maximum state
probability discrepancy of `3.34e-16`.

These checks establish the small-fixture Gaussian identities and market
factorization only. They do not validate a Stan implementation, production
gradients, policy reconstruction, or sampling performance.
