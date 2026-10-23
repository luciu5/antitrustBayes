# Logit and CES inside-market utility intercept

Version `0.2.0.9007`, model specification 12, corrected the single-market
Logit utility normalization in both simultaneous and ordered-k kernels.
Version `0.2.0.9008`, model specification 13, applies the same correction to
both CES kernels; the sections below describe the Logit derivation and then
state the CES translation.

The outside good has utility zero. Product qualities are centered within their
quality block, and the price regressor retains its existing standardization.
Every market also has an inside-market intercept. For an output market:

\[
\log r_{jmt}=\operatorname{logit}(s_{0,mt})+a_m+b_{jm}
             -\lambda_m\widetilde p_{jmt}+\epsilon_{jmt}.
\]

Here `r` is the conditional inside share. The intercept measures the overall
utility of inside goods relative to the outside option. Centering product
qualities does not set this economic parameter to zero.
The standardized-price coefficient is `lambda_m = alpha_m * price_sd`, where
`alpha_m` is the physical price coefficient used in the markup equations below.
Stan reports these as `a_market` and `alpha_structural_market`, respectively.

## The removed restriction

For one product per firm in a noiseless single-market/single-period example,
centered qualities and centered prices imply

\[
\overline{\log r}=\operatorname{logit}(s_0)+a.
\]

The former single-market branch fixed `a = 0`. It therefore forced
`s0 = plogis(mean(log(r)))`. The geometric mean of seven positive conditional
shares summing to one is at most `1/7`, so the forced outside share was at
most `1/8`. With noisy observations this becomes an artificial likelihood
restriction rather than an exact equality.

The corrected intercept is free. At any chosen `0 < s0 < 1`,
`a = mean(log(r)) - qlogis(s0)` satisfies the mean demand equation without
changing the conditional shares. Centered product qualities explain the
remaining cross-product variation. No centered-log-ratio likelihood is used.

## Parameter and prior contract

For one market, `mu_b_market_vec[1]` is the intercept itself, with exactly one
normalized `Normal(0, 1)` prior. `sigma_b_market_vec` and `b_market_raw` have
zero length. `market_intercept_location` and sampled `b_market` contain that
intercept. A request for collapsed market-intercept coordinates still resolves
to the sampled single-market intercept; it does not fix the intercept.

Multiple markets retain the existing population mean, hierarchical market
deviations and selected SD prior. Their optional analytical integration is
unchanged. `is_single_market` keeps its original value and continues to govern
the existing alpha and outside-share pooling parameterizations. Outside-share
pooling and intercept estimation are separate decisions.

The existing `gamma_market` reporting convention now also applies to one
market: it reports the intercept with one period, and intercept plus mean
outside-share log odds with multiple periods. `c_unbounded_mky` and
`c_effective_mky` continue to include the intercept once. The positive
`a_market` field remains the price coefficient, despite its historical name.

## Identification and posterior use

Conditional demand identifies the sum of intercept and outside-share log odds.
Supply information and proper priors can separate them. Under monopolistic
competition, dollar markups are `1/alpha`, independent of `s0`; shares alone
therefore cannot precisely identify `s0` with a free intercept. Its posterior
can remain sensitive to the prior.

For single-product Cournot, at conditional share `r_j`, the dollar markup is

\[
p_j-c_j=\frac{1}{\alpha}\left\{1+\frac{1-s_0}{s_0}r_j\right\}.
\]

Heterogeneous shares and accurate margins identify a slope and intercept in
this relationship, providing information about both `alpha` and `s0`. This is
an identification result, not a promise of precise recovery from every small,
noisy sample.

Posterior quality, predictive and counterfactual code uses the fitted
`share_mean`, which already includes the intercept. It must not add an
independent intercept to that mean. Economic model construction preserves the
posterior structural outside share, normalizes the within-market-year quality
scores and retains existing physical-cost and role rules.

Both CmdStan and Nutpie compile the same corrected installed kernels. The
normalized target includes the additional scalar and its proper prior;
bridge sampling and conditional-k evidence retain their existing mathematics.
Source hashes and the specification identify the new target. Earlier fitted
objects require refitting before current comparison or simulation. There is
no legacy normalization switch.

## CES translation (specification 13)

The two CES kernels previously kept the zero single-market intercept that the
Logit correction removed. They now follow the Logit contract exactly:

- Outside-good utility is zero and product qualities keep their symmetric
  sum-to-zero coordinates. For one market, `mu_b_market_vec` has length one and
  *is* the inside-market intercept, with exactly one normalized `Normal(0, 1)`
  prior. `sigma_b_market_vec` and `b_market_raw` have zero length; no
  single-market hierarchy SD exists. Multiple markets keep the population mean,
  hierarchical deviations, SD prior and optional analytical integration.
- Conditional quantity-share demand is
  `log r = logit(s0) + a + b_j + (alpha + sd_log_p) z_j + gamma log(p_0)`, where
  `z` is the standardized log price, `gamma = 1 + alpha / sd_log_p`, and `a` is
  the free intercept (market-specific for multiple markets). Only the
  intercept's value changed for one market; the price and outside-price terms
  are as before.
- `s0` remains the *quantity* outside share, which the intercept and the
  outside-share log odds enter additively. The CES *revenue* outside share
  `R0 = plogis(qlogis(s0) + log(p_0) - log(inside price index))` still drives
  the supply first-order conditions, so margins depend on `R0`, not on `s0`
  directly. Outside prices other than one therefore shift `R0` and the
  `gamma log(p_0)` term but not the sum-to-zero qualities or the intercept prior.
- `b_market_rot` (the CES counterpart of Logit's `gamma_market`) now reports the
  intercept for one market and one period, and the intercept plus the mean
  outside-share log odds when a market has several periods.
  `c_unbounded_mky` and `c_effective_mky` include the intercept once.
- Posterior quality, predictive and counterfactual code already used the fitted
  `share_mean`, so no independent intercept is added anywhere in R.

Executed checks are in `tests/testthat/test-ces-intercept-normalization.R`:
exact recovery of the seven-firm conditional shares at 10%, 20% and 40% outside
shares and outside prices 0.8, 1 and 1.4; fixed-parameter Stan evaluation of
both CES kernels against `antitrust`/`coordination` margins; the normalized
prior; and posterior-to-counterfactual identities for B, C, CF and CS.

## Effect of the global specification bump (12 to 13)

`.antitrust_bayes_model_spec_version` is a single package-wide integer. The CES
correction therefore advances it for *every* fit, including Logit fits whose
kernels did not change between specifications 12 and 13. Consequences:

- Any fitted object, saved fit or cached comparison created at specification 12
  (Logit or CES) or earlier is rejected by `validate_fit_spec()` and by the
  bridge-evidence validation until it is refit with `bayes_est()`. Logit
  posteriors are not numerically invalidated by this change; the rejection is a
  consequence of the single global guard and was chosen over per-kernel version
  bookkeeping.
- Downstream consumers that pin or assert `model_spec_version == 12`
  (for example the `mergerBayes` and banking runners) must update the pin to 13
  and refit. No compatibility shim was added.
- Stan source hashes change only for the two CES kernels, so cached compiled
  Logit executables remain valid; CES executables are rebuilt on first use.

## Integration scope

The `mergerBayes` revised 20% outside-share baseline is economically compatible
with this specification. Its later-period equilibrium shocks can change the
true outside share, while market-level outside-share pooling imposes one value
across periods. That separate pooling assumption is unchanged here. No
production experiment or DGP tuning is part of this correction.

Executed verification is recorded in
[the correction validation report](validation/single_market_intercept_results.md).
