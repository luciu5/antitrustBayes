# Market-intercept scale-coordinate experiment

## Question

The Gaussian panel model contains a market-intercept scale
\(\sigma_b\) and, under market-level outside-share pooling, an outside-share
scale \(\sigma_{s0}\).  This experiment tested whether replacing those two
positive sampling coordinates by

\[
H=\sqrt{\sigma_b^2+\sigma_{s0}^2},\qquad
\lambda=\frac{\sigma_{s0}^2}{\sigma_b^2+\sigma_{s0}^2}
\]

would materially improve HMC without changing the statistical model.  The
inverse map is

\[
\sigma_b=H\sqrt{1-\lambda},\qquad
\sigma_{s0}=H\sqrt{\lambda},
\]

with constrained-coordinate Jacobian

\[
\left|\frac{\partial(\sigma_b,\sigma_{s0})}
{\partial(H,\lambda)}\right|
=\frac{H}{2\sqrt{\lambda(1-\lambda)}}.
\]

The experiment retained the existing priors by evaluating them on the
reconstructed scales and adding this Jacobian.  It did not alter demand,
supply, pooling, or likelihood equations.

## Analytically and numerically established

The transformation is exact.  Fixed-state checks used the native compiled
Logit/CES supply and Stackelberg kernels, including missing margins, active
\(\rho\), collapsed market quality, and exact ordered-k likelihoods.  All
checks passed:

| check | passed | maximum error |
|---|---:|---:|
| normalized Stan target | 20/20 | 4.73e-12 |
| reconstructed/generated state | 150/150 | see CSV |
| automatic gradient vs finite difference | 12/12 | 6.52e-09 scaled |

The mapped unconstrained targets differ by the expected state-independent
\(-\log 2\), because
\(|\partial(\log\sigma_b,\log\sigma_{s0})/
\partial(\log H,\operatorname{logit}\lambda)|=1/2\).  This constant was
included in the bridge-normalized comparison.  Detailed results are in:

* `market_intercept_scale_target_checks.csv`
* `market_intercept_scale_state_checks.csv`
* `market_intercept_scale_gradient_checks.csv`

## Controlled HMC evidence

The local gate used one prespecified Logit-Bertrand fixture with 6 markets,
3 periods, 5 firms per market, active rho, Gaussian margins and shares,
collapsed market quality, market-level outside-share pooling, two chains,
600 iterations (300 warmup), `adapt_delta=.99`, and a common seed.  It compared
the historical collapsed coordinates, the total/fraction coordinates, and the
existing sampled-intercept joint rotation.  Two DGP regimes were run: a true
zero market-intercept dispersion and a true interior dispersion of 0.18.

| DGP | coordinates | seconds | divergences | max Rhat | min bulk ESS | min BFMI |
|---|---|---:|---:|---:|---:|---:|
| boundary | collapsed legacy | 46.2 | 9 | 1.028 | 64 | .506 |
| boundary | collapsed total/fraction | 59.4 | 6 | 1.045 | 58 | .610 |
| boundary | sampled joint | 265.7 | 0 | 1.021 | 143 | .854 |
| interior | collapsed legacy | 72.0 | 0 | 1.018 | 78 | .518 |
| interior | collapsed total/fraction | 69.7 | 0 | 1.033 | 91 | .548 |
| interior | sampled joint | 173.7 | 0 | 1.019 | 112 | .753 |

These short runs are geometry diagnostics, not a recovery study.  They show
that the severe divergence pattern is concentrated at the variance-component
boundary.  They do not support a claim that collapsed intercepts generally
fail as panels grow.

## Decision

The total/fraction coordinates are rejected as a production API.  They were
slower in the boundary comparison, did not eliminate divergences, and did not
improve Rhat or ESS consistently.  The small reduction from 9 to 6
divergences is not worth another permanent coordinate mode.

The existing architecture is retained:

* sampled market intercepts remain the default and robust fallback;
* collapsed market intercepts remain an opt-in speed path whose suitability
  must be diagnosed for the application;
* the sampled joint market-s0 rotation remains the leading fallback to test
  when the intercept variance may be near zero; this short run removed
  divergences but did not pass every convergence gate;
* the enterprise synthetic harness now covers both boundary and interior true
  market-intercept dispersions, because prior fixtures inadvertently tested
  only the zero-boundary case.

The HMC rows and failed gates are retained in
`market_intercept_scale_hmc.csv` and the per-case directory.  No sampler gate
was weakened, and no prior, likelihood, or economic equation was changed.

## Remaining limits

The runs use one seed and a short local HMC budget.  The interior fixture does
not strongly identify the variance decomposition in this sample; its posterior
scale median should not be read as a recovery result.  Future enterprise tests
should cross boundary/interior truth with Logit/CES, atomic/exact-k conduct,
small panels, and larger panels.  The evidence here is sufficient to reject
the new coordinate API, not to choose one market-intercept representation for
every dataset.

The rejected coordinate implementation was removed from the production Stan
and R interfaces.  The model-specification version therefore remains 11; old
serialized fits are not invalidated by this experiment.
