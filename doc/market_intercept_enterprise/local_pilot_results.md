# Local market-intercept pilot results

## Scope

These runs compare two mathematically equivalent representations of the
market intercept while holding the economic model, priors, likelihood,
market-specific collapsed quality, active `rho`, data, seed, and HMC budget
fixed:

- `market_intercept_parameterization = "collapsed"` analytically integrates
  the market intercept;
- `market_intercept_parameterization = "sampled"` explicitly samples it and
  uses the joint market-intercept/outside-share coordinates.

Each fit used two chains, 1,200 iterations per chain, 600 warmup iterations,
`adapt_delta = 0.99`, and `max_treedepth = 13`. The budget is a diagnostic
pilot, not a production posterior. The compact machine-readable results are
in `local_pilot_hmc.csv`; incremental Stan CSVs remain in the ignored
`doc/statistical_toolkit_runs/20260922` directories.

## Executed evidence

### Ten-market cross-sections and panels

| Demand/conduct | Periods | Intercept | Wall seconds | Divergences | Treedepth | Min BFMI | Max R-hat | Min bulk ESS | Gate |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| Logit B | 1 | collapsed | 55 | 0 | 0 | 0.653 | 1.008 | 228 | pass |
| Logit B | 1 | sampled | 86 | 0 | 0 | 0.704 | 1.004 | 460 | pass |
| Logit B | 4 | collapsed | 251 | 8 | 0 | 0.263 | 1.058 | 31 | fail |
| Logit B | 4 | sampled | 1,048 | 0 | 0 | 0.881 | 1.004 | 501 | pass |
| Logit BS | 4 | collapsed | 670 | 56 | 0 | 0.238 | 1.007 | 56 | fail |
| Logit BS | 4 | sampled | 2,568 | 0 | 0 | 0.844 | 1.005 | 499 | pass |
| CES C | 4 | collapsed | 641 | 0 | 0 | 0.279 | 1.0101 | 63 | near-pass |
| CES CS | 4 | collapsed | 1,764 | 0 | 0 | 0.396 | 1.003 | 93 | near-pass |

The CES C fit missed the frozen gate on BFMI, R-hat, and bulk ESS by small
amounts under the short budget. The CES CS fit missed only the bulk-ESS gate
of 100. Neither had divergences or treedepth hits.

### Two-market, two-period sparse panels

| Demand/conduct | Intercept | Wall seconds | Divergences | Treedepth | Min BFMI | Max R-hat | Min bulk ESS | Gate |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Logit B | collapsed | 25 | 0 | 0 | 0.773 | 1.004 | 384 | pass |
| Logit B | sampled | 72 | 0 | 0 | 0.883 | 1.007 | 398 | pass |
| Logit BS | collapsed | 53 | 0 | 0 | 0.856 | 1.005 | 327 | pass |
| Logit BS | sampled | 157 | 0 | 0 | 0.705 | 1.004 | 462 | pass |
| CES C | collapsed | 25 | 0 | 0 | 0.757 | 1.009 | 328 | pass |
| CES C | sampled | 104 | 3 | 0 | 0.771 | 1.011 | 327 | fail |
| CES CS | collapsed | 53 | 18 | 0 | 0.851 | 1.011 | 212 | fail |
| CES CS | sampled | 287 | 0 | 0 | 0.807 | 1.014 | 415 | fail on R-hat only |

The exact-k posterior probabilities summed to one to machine precision in all
reported BF/CF/BS/CS fits.

## Posterior equivalence

For each collapsed/sampled pair, reported economic posterior means were
compared using the combined Monte Carlo standard error. The maximum absolute
difference was below 1.8 combined MCSE for the usable pairs. The CES C pair
reached 2.8 MCSE in `mu_log_a`; its sampled fit had three divergences and was
not usable. The 10-market Logit B pair differed by at most 1.5 combined MCSE
despite the collapsed fit's poor diagnostics. These checks support a geometry
difference rather than a changed statistical target, but fixed-state target
equality remains the stronger validation.

The independent code/target audit found that the collapsed kernel contains
the full Gaussian quadratic form, log determinant, residual-scale
normalization, correlation adjustment, normalized priors, and exact-k
mixture. Existing reference checks report a dense-density maximum error of
`7.1e-15`, posterior-mean error of `5.1e-16`, native normalized-target error
of `3.9e-10`, and 165/165 passing automatic-differentiation checks. The
sampled joint s0 rotation is orthogonal with determinant one. These results
support exact target equivalence and do not support an implementation error.

## Divergence localization

Post-warmup divergent draws in the 10-market Logit panels concentrate where
outside-share dispersion dominates market-intercept dispersion. Define

\[
S=\sqrt{\sigma_{s0}^2+\sigma_b^2},\qquad
\phi=\frac{\sigma_{s0}^2}{\sigma_{s0}^2+\sigma_b^2}.
\]

For ordinary Logit B, median `phi` was 0.949 in divergent draws and 0.670 in
nondivergent draws, while median `S` was nearly unchanged (0.127 versus
0.125). For Logit BS the contrast was 0.996 versus 0.726; the divergent median
market-intercept scale was 0.0086 versus 0.0639 otherwise. This is direct
evidence of a variance-fraction boundary. It does not by itself establish
whether the cause is inherent marginal variance-component geometry or an
implementation defect.

## Current assessment

Analytically established:

- the compared representations are intended to preserve the same economic
  model and priors;
- the market-intercept and market-level outside-share effects enter the share
  mean through a common market direction;
- collapsing a Gaussian intercept yields a low-rank covariance contribution.

Numerically verified:

- sparse 2x2 Logit panels are usable under both representations;
- collapsed intercepts are roughly three to four times faster in the clean
  paired cases;
- sampled intercepts remove the observed large-panel Logit divergences;
- the large-panel collapsed Logit divergences concentrate near `phi = 1`;
- CES behavior is conduct-specific rather than uniformly favoring one
  representation.

Still unresolved:

- whether total-scale/variance-fraction coordinates eliminate the collapsed
  boundary without changing the target;
- whether a prespecified automatic selection rule is preferable to making one
  representation universal;
- confirmation across additional seeds after the geometry is repaired.

No Supreme Court fit is used as a validity gate in this pilot.

## Candidate exact repair

If the speed advantage of collapse is to be retained, sample the identified
total scale and weakly identified allocation directly:

\[
H=\sqrt{\sigma_b^2+\sigma_{s0}^2},\qquad
\lambda=\frac{\sigma_{s0}^2}{H^2},
\]

and reconstruct

\[
\sigma_b=H\sqrt{1-\lambda},\qquad
\sigma_{s0}=H\sqrt{\lambda}.
\]

For constrained `(H, lambda)` coordinates, the required Jacobian is

\[
\log H-\log 2-\tfrac12\log\lambda-\tfrac12\log(1-\lambda).
\]

The existing independent half-normal priors must be evaluated on the two
reconstructed scales; placing new independent priors on `H` and `lambda`
would change the model. This repair remains experimental until fixed-state
target/Jacobian equality, boundary gradients, and paired 2x2/10x4 HMC gates
pass. Until then, the sampled-intercept joint-s0 representation is the
validated robust option for the difficult large Logit panels.
