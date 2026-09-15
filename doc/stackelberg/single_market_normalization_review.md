# Single-market Bayes normalization review

Date: 2026-09-13. No production files were changed. The independent check is
[`/tmp/single_market_normalization_check.R`](/tmp/single_market_normalization_check.R).

The proposed adjustment is algebraically correct for the existing
single-market Stan equations, provided the same centering conventions are
used: `b_market=0`, and `b_quality_all` is centered over unique firms by
`b_group_raw - mean(b_group_raw)`.

For Logit, with signed price coefficient `sign*a` (`sign=-1` output,
`sign=+1` input), define

```
deltaNorm_j = log(qc_j) - sign*a*(p_j-p0)
eta0 = mean_unique_firm(deltaNorm) - sign*a*(p0-mean_all(p))
bq_j = deltaNorm_j - mean_unique_firm(deltaNorm).
```

Then `eta0 + bq_j + sign*a*(p_j-mean_all(p)) = log(qc_j)` exactly, so the
Stan predictor `logit(s0)+delta` reproduces the conditional share logs when
`s0=plogis(eta0)`. The audit's `+a*(p-p0)` and `+a*(p0-mean(p))` are exactly
the output case. Input Logit requires the sign-coded minus version; using the
plus version there would be a sign error.

For CES output, with `gamma=g`, define

```
lnQNorm_j = log(qc_j) + g*log(p_j/p0)
eta0 = mean_unique_firm(lnQNorm) - g*mean_all(log(p))
bq_j = lnQNorm_j - mean_unique_firm(lnQNorm).
```

The existing predictor then satisfies
`eta0+bq_j-g*(log(p_j)-mean(log(p)))+g*log(p0)=log(qc_j)` exactly.
The R check obtains errors below `2e-15` for Logit output/input and CES.

These formulas require the normalized quality quantity (`deltaNorm` or
`lnQNorm`) to be constant for a firm across all products and periods that
share one `b_quality` draw. A fixed-quality DGP satisfies this. If qualities
vary by year or by product beyond the shared firm effect, a single-market
centered model cannot reproduce all rows exactly; an additional effect or
residual is then required.

The implied outside share is generally not an arbitrary input. In the check,
the output Logit data imply `s0=0.158322...`, so forcing `s0=.30` conflicts
with the centered no-intercept equation except for special data. For CES,
the outside-share formula changes when every price and the outside price are
rescaled by the same unit factor (from `0.092839` to `0.000407` under a
10-fold rescaling), while predicted share logs remain exact. This confirms
the stated unit-normalization caveat: a free market intercept could absorb
that common shift; `b_market=0` cannot.

## Ordered-k addendum

The `ordered_k` specification is mathematically coherent. A market-scoped
count must be summed over all years before `log_sum_exp`; multiplying
period-level mixtures would define a different model. Top-k role sets can
split exact size ties because this prior is a count prior rather than a
scalar boundary; no scalar cutoff should be inferred from such a split.
The stated restrictions on LOO units, shared simulation state, and
market-year/market aggregation follow from the shared latent `k`. Required
tests include two-period differing preferred counts, interleaved rows, ties,
zero/all leaders, absent margins, and equality of whole-unit target and GQ
mixture calculations.
