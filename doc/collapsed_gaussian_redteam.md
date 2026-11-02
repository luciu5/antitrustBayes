# Collapsed Gaussian red-team review

Review date: 2026-09-16. Scope: current uncommitted collapsed Gaussian,
Gaussian-rho, exact ordered-k, posterior-quality, and s0-pooling changes. This
was an independent source review. No production files, Git refs, or index were
changed by the review.

## Gate result

No release-blocking defect was found in the normalized collapsed target, the
stable sparse sufficient-statistic implementation, exact-k mixture accounting,
the s0 automatic ladder, or the benchmark FOC oracle. The primary benchmark can
continue under its frozen settings. The policy reconstruction limitation below
remains material and should stay explicit in claims and documentation.

## Findings and evidence

### Collapsed Gaussian algebra and numerical stability

The four Stan programs use the same centered Helmert model as the sampled path:

- `inst/stan/logit_supply.stan:63-118` and the corresponding helpers in
  `ces_supply.stan`, `logit_stackelberg.stan`, and `ces_stackelberg.stan` keep
  positive within-group weighted-regression contributions, group sufficient
  statistics `(a,c)`, and the centered rank-`G-1` determinant correction.
- `inst/stan/logit_supply.stan:121-150` adds Gaussian-rho rows with loading
  `-rho*sigma_margin/sigma_share` and variance
  `sigma_margin^2*(1-rho^2)`. This is the unit-Jacobian triangular
  transformation of the joint Gaussian residual block.
- The collapsed parameter declarations use zero-width quality coordinates and
  avoid a dense covariance or Helmert matrix in the collapsed path; the
  conditional Helmert matrix is allocated only for sampled quality.

The prior, Gaussian normalization constants, centered determinant term, and
rho transformation prior/Jacobian are present. The stable update was checked
against dense calculations and tiny-share-noise cases; the prior cancellation
failure is gone.

### Gaussian rho and missing margins

The direct data paths mask missing margins before transforming them and retain
only observed indices (`R/helpers.R:525-536`, `R/helpers.R:821-838`). Active
collapsed rho with multiple markets is rejected where global quality centering
would otherwise create unsupported cross-market coupling
(`R/helpers.R:474-485`, `R/helpers.R:877-890`). Exact ordered-k active rho is
also restricted to one market, matching the Stan target.

`posterior_quality()` now defaults to the full matched Gaussian conditional for
standard/smooth-fringe rho fits and offers an explicit share-only option
(`R/posterior_quality.R:272-310`, `R/posterior_quality.R:326-340`). It rejects
exact ordered-k rho recovery because generated margin predictions do not retain
component-specific means. This is the correct safe behavior; the share-only
conditional must not be described as the full rho posterior.

### Exact ordered-k mixtures

The stack programs compute one common collapsed share contribution and add each
component's margin contribution only once before the market-level
`log_sum_exp` (`inst/stan/logit_stackelberg.stan:572-675`, with the analogous
CES path). The resulting posterior component probabilities and generated
market margin likelihood use the same normalized mixture target. Ordered-k
support is validated as market-local (`R/estimate.R:343-363`,
`R/stackelberg_bayes.R:344-417`). No duplicate share likelihood or missing
mixture normalization was found.

### s0 ladder and economic equations

The automatic ladder resolves one direct outside share for one market-year and
one market-level outside share for panels/multiple markets
(`R/helpers.R:69-140`). The direct one-market Stan declarations omit unused
hierarchical coordinates and use `mu_s0_pool_vec[1]` consistently. The legacy
flags remain separately supported.

The production Logit and CES inverse-margin equations are unchanged in economic
meaning. I independently compared the benchmark's simultaneous, fringe, and
Stackelberg Bertrand/Cournot equations (`doc/collapsed_gaussian_benchmark.R:127-189`)
against the production formulas in the Stan programs over 1,000 random cases
per cell. Maximum absolute errors were `2.22e-16` for Logit and `8.88e-16` for
CES. The reproducible check is
`doc/collapsed_gaussian_oracle_validation.R`, with results in
`doc/collapsed_gaussian_oracle_validation.csv`.

### Material policy-reconstruction limitation

The deterministic policy path still anchors on observed conditional shares:
`R/simulate.R:231-281` extracts common draw parameters but no latent quality,
and `R/stackelberg_bridge.R:329-417,473-577` records
`observed_conditional_shares` as its anchoring method. Thus posterior quality
uncertainty does not enter those policy reconstructions. The architecture
document acknowledges this at `doc/collapsed_gaussian_architecture.md:301-310`.
The implementation should retain this limitation in user-facing policy claims;
addressing it would require a coordinated simulation redesign and was outside
this review.

### Predictive generated quantities

Collapsed Stan generated quantities intentionally report the common-parameter
mean prediction, `pred_shareIn = exp(share_mean)`, while sampled quality paths
include sampled `b_quality` (`inst/stan/logit_supply.stan:560-581` and the
corresponding three programs). The benchmark adapter now uses
`posterior_quality()` for comparable collapsed predictive summaries
(`doc/collapsed_gaussian_benchmark.R:791-798`). Raw generated-quantity values
should continue to be described as common-mean predictions.

## Verification limits

This review did not run HMC or a new compilation while the primary benchmark
was running. Independent reference checks, Stan source checks, and focused
posterior-quality checks supplied by the working session passed. Exact CF/CS
HMC cells are not part of the frozen 21-cell primary plan; source and algebra
coverage includes all four ordered-k conduct labels.
