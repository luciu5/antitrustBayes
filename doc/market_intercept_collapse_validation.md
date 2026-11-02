# Gaussian market-intercept collapse validation

## Scope

Specification 11 adds an opt-in `market_intercept_parameterization =
"collapsed"` path. It analytically integrates the existing Gaussian market
intercept together with centered market-scoped Gaussian quality. It does not
change priors, outside-share pooling, demand, conduct, exact-k support, or
counterfactual definitions. The routine default remains `"sampled"` until the
frozen HMC comparison is complete.

The supported experimental path requires multiple markets, Gaussian share and
margin likelihoods, collapsed market-scoped quality, a joint likelihood, and
base (`legacy`) s0 coordinates. A single-market request is an explicit no-op
because the existing model fixes its market intercept to zero.

## Analytically established

For a market, margin-first conditioning gives a share residual vector `z` and
diagonal covariance `D`. Paired rows have variance
`sigma_share^2 * (1 - rho^2)`; rows with missing margins have variance
`sigma_share^2`. Integrating centered quality `q` and the zero-mean market
intercept offset `b - mu_b` gives

```text
V = D + tau^2 Z P Z' + sigma_b_market^2 11'.
```

The production density uses group sufficient statistics and determinant/
Woodbury identities in `O(N + G)` time and memory. It retains every Gaussian
normalizing constant. Conditional reconstruction first draws the market
intercept from its marginal conditional posterior and then draws centered
quality conditional on that same intercept. In ordered-k models, both draws
use one market-level k realization from the fully integrated component
probabilities.

## Numerically verified

`market_intercept_collapse_extended_reference.R` compares the structured
formula with dense covariance calculations and an explicit three-market,
three-state enumeration oracle.

| Check | Result |
|---|---:|
| Structured density cases | 3 / 3 passed |
| Joint posterior cases | 3 / 3 passed |
| Exact-k cases | 1 / 1 passed |
| Maximum density error | 7.105427e-15 |
| Maximum posterior mean error | 5.065393e-16 |
| Maximum posterior covariance error | 2.081668e-17 |
| Maximum q-given-b formula error | 6.175616e-16 |
| Enumerated joint k states | 27 |
| Exact-k log-evidence error | 0 |
| Maximum exact-k probability error | 3.330669e-16 |

`market_intercept_collapse_native_validation.R` extracts the exact production
Stan helpers into a small native model. It checks normalized targets, including
the production Gaussian intercept, quality-scale, residual-scale, and rho
hyperpriors, and Stan automatic derivatives against an independent dense R
oracle. Cases include small share noise, small quality/intercept scales, large
positive and negative rho, missing margins, and one market with no margins.

| Check | Result |
|---|---:|
| Native normalized target cases | 5 / 5 passed |
| Maximum target error | 3.860806e-10 |
| Native AD comparisons | 165 / 165 passed |
| Maximum scaled gradient error | 5.072867e-6 |

The larger gradient discrepancy occurs in the deliberately difficult
small-share-noise case and remains inside the prespecified absolute-plus-
relative tolerance. The validation script compiles a helper harness and does
not claim a whole-model HMC result. SHA-256 extraction checks confirm that the
three integrated-intercept helper bodies exercised by the native harness are
byte-identical across the Logit/CES simultaneous and exact-k kernels.

All four production Stan programs pass `rstan::stanc`. The generated bindings
compile and the package installs. The complete non-sampling test suite and the
new exact-k/reconstruction tests pass. The bridge migration regression suite
also passes with the new kernels registered as an explicit forward-only stage;
it injects sampled-intercept code 1 only for approved historical sources and
rejects reduced-dimension collapsed-intercept sources. Sampling tests remain
separately reported because they are opt-in and materially more expensive.

## Independent review

A Luna xhigh reviewer who did not author the implementation found no target or
economic-equation defect. The review confirmed:

- market likelihoods integrate the intercept and centered quality once;
- exact-k component scores retain the whole market contribution;
- missing margins retain unconditional share variance;
- one-group markets deterministically set centered quality to zero;
- sampled and collapsed coordinates do not become duplicate model-universe
  hypotheses;
- reconstruction uses one coherent k, intercept, and quality state;
- policy plans retain and fingerprint reconstructed intercept draws.

Inside-share score normalization in latent-quality policy simulation removes a
common market intercept while the fitted outside share fixes the extensive
margin. This is the existing policy-anchoring convention, now documented and
tested rather than a new economic assumption.

## Sampler evidence still required

One ad hoc two-market Logit smoke fit executed the new code path through four
default 2,000-iteration chains and reported four total divergences. The command
incorrectly supplied sampler arguments through `...` instead of the public
`stan_options` list, so the intended short budget and control settings were not
used; the process did not retain a complete diagnostic table. This failed smoke
is retained here as execution evidence only. It is neither a frozen comparison
nor support for promotion. The benchmark runner passes every sampler setting
through `stan_options` and records all chains.

The required two-seed Supreme Court comparison has not been executed. Until it
passes the frozen divergence, treedepth, BFMI, Rhat, ESS, posterior-equivalence,
and policy-equivalence gates, the integrated intercept remains experimental
and the package default must not change.

The runner validates the canonical 176-row, six-market fixture without
sampling by default. Run mode verifies the SHA-256 hashes of all four installed
Stan kernels before starting, records those hashes and the dirty source diff,
and initializes every case at the common zero-unconstrained economic center.
