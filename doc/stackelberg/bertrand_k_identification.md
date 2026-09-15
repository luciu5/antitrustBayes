# Logit Bertrand-Stackelberg ordered-k identification audit

Status: conditional numerical audit, 2026-09-15. This note addresses whether
the seed-1103 single-product recovery result reflects a scale, normalization,
or ordering error. It does not change the prespecified DGP, prior, likelihood,
or sampler gate.

## Analytically established

The Stan margin response is `log_marginNormInv = log(1/margin)`. The DGP uses

    measured_margin = true_margin * exp(-epsilon),

so its observed response is exactly

    log(1/measured_margin) = log(1/true_margin) + epsilon.

The DGP draws `epsilon` from a Student-t with five degrees of freedom and scale
`0.04*sqrt(3/5)`. Its standard deviation is therefore 0.04. The Stan likelihood
uses the identical scale conversion. There is no response-scale conversion
between the DGP and likelihood.

For each period and candidate k, the audit aggregates observed conditional
shares to firms, converts them to unconditional shares using the true outside
share, assigns the fixed whole-market rank, and evaluates the audited exact
Bertrand-Stackelberg equations. Product rows and all covariates remain in their
original order. Only the firm-level role partition changes.

## Numerically verified at the true continuous parameters

The reproducible calculation is in `bertrand_k_identification.R`; its tables are
`bertrand_k_oracle_scores.csv`, `bertrand_k_oracle_power.csv`,
`bertrand_k_zero_noise.csv`, and `bertrand_k_layout_scores.csv`.

For the realized true-k=2 sample, the exact conditional log likelihoods are:

| k | noiseless RMSE from k=2 | RMSE / margin-noise SD | log likelihood relative to k=2 | conditional probability |
|---:|---:|---:|---:|---:|
| 1 | 0.012140 | 0.304 | 0.8624 | 0.4767 |
| 2 | 0 | 0 | 0 | 0.2012 |
| 3 | 0.006916 | 0.173 | 0.0372 | 0.2088 |
| 4 | 0.009992 | 0.250 | -0.8575 | 0.0854 |
| 5 | 0.012711 | 0.318 | -1.9752 | 0.0279 |

Thus adjacent partitions change the log inverse-margin prediction by only
0.17--0.32 of one allowed noise standard deviation in RMS. The realized noise
happens to favor k=1 over true k=2 by 0.86 log units. This conclusion is
available before HMC and does not depend on sampler behavior.

With the realized noise removed but the likelihood SD retained at 0.04, the
conditional probability of k=2 is 0.689. Holding the same design fixed and
using SD 0.02 raises it to 0.989. This is a power comparison, not a proposal to
retune the fitted model after seeing its result.

A prespecified 5,000-replication noise experiment (seed 20260915) makes k=2 the
modal state in 59.16 percent of replications. Mean posterior probability on k=2
is 0.4257 and median probability is 0.4180. The design is informative about k,
but it is not expected to identify the true state decisively in every eight-
period realization.

## Empirically estimated by the current HMC fit

The completed matching fit gives component probabilities
`(0.4832, 0.1967, 0.2110, 0.0823, 0.0269)`, close to the true-parameter oracle
above. It recovers alpha at 2.99993 versus truth 3, the outside share near
0.153451 versus truth 0.153451, and the margin-noise SD at 0.0389 versus truth
0.04. This agreement rejects a material alpha, outside-share, unconditional-
share, margin-response, or role-order scaling mismatch in this case.

That fit is not yet a passing recovery draw: margin-noise R-hat is about 1.025
and effective sample size is about 97 under the prespecified two-chain run.
Those diagnostics limit precision, but they do not explain why k=1 is modal;
the exact oracle already has the same ordering.

The two other prespecified ownership layouts select the true state.  At the
true continuous parameters, the multiproduct-leader design assigns probability
0.6840 to k=2; HMC assigns 0.6184 and has posterior mode and median 2.  The
multiproduct-follower design assigns oracle probability 0.9001 to k=2; HMC
assigns 0.9050 and again has mode and median 2.  The HMC and oracle results are
therefore close in all three layouts:

| layout | oracle P(k=2) | HMC P(k=2) | HMC mode |
|:--|--:|--:|--:|
| single product | 0.2012 | 0.1967 | 1 |
| multiproduct leader | 0.6840 | 0.6184 | 2 |
| multiproduct follower | 0.9001 | 0.9050 | 2 |

These compact HMC fits all fail the prespecified sampler gate, so their
posterior probabilities remain exploratory.  The oracle calculations do not
use HMC and establish the same cross-layout pattern.

## Sampler geometry of the exact-share experiment

The primary recovery DGP deliberately has zero share-equation residual.  That
choice places the data-generating share-error standard deviation on the
boundary of the fitted model, whose parameter is strictly positive.  The
result is visible in every candidate structure because the share likelihood is
common to every k state.

The scale collapse has a quantitative check.  With x = log(sigma_share), a
Normal prior centered at log(.5) with SD .5, N exact observations, and p locally
fitted regression coordinates, the zero-residual scale kernel has approximate
mode

    x = log(.5) - (N - p) * .5^2.

There are p=13 effective conditional-demand coordinates here: one alpha, eight
market-year outside-share intercepts, and four centered firm-quality
coordinates.  The single-product layout has N=40, giving -7.443; its sampled
mean is -7.497.  Each multiproduct layout has N=48, giving -9.443; their sampled
means are -9.430 and -9.410.  This agreement diagnoses the small share scale as
the expected consequence of exact-fit data rather than a conduct-equation
scale error.

The collapsed scale creates difficult posterior geometry and low effective
sample size.  It does not determine k directly: the common share likelihood
cancels from the conditional k odds when rho is fixed at zero.  A separately
prespecified compositional-share-noise stress experiment is the appropriate
way to test sampling away from this boundary.  It cannot replace or relabel
the failed primary sampler gate.

## Still unresolved

Cournot and CES recovery remain to be established.  A useful future DGP-design
exercise can vary the number and precision of margin observations
prospectively. It must remain separate from the current prespecified recovery
gate and from the PNB application.
