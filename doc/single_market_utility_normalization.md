# Logit inside-market utility intercept

Version `0.2.0.9007`, model specification 12, corrects the single-market
Logit utility normalization in both simultaneous and ordered-k kernels.

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
Source hashes and specification 12 identify the new target. Earlier fitted
objects require refitting before current comparison or simulation. There is
no legacy normalization switch. The CES kernels are unchanged.

## Integration scope

The `mergerBayes` revised 20% outside-share baseline is economically compatible
with this specification. Its later-period equilibrium shocks can change the
true outside share, while market-level outside-share pooling imposes one value
across periods. That separate pooling assumption is unchanged here. No
production experiment or DGP tuning is part of this correction.

Executed verification is recorded in
[the correction validation report](validation/single_market_intercept_results.md).
