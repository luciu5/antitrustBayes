# Exact synthetic shares and share-noise concentration

This note concerns the prespecified synthetic recovery DGP, whose prices and
shares are exact and whose measured log inverse margins have independent
Student-t noise. It is separate from the Stackelberg FOC audit. It changes no
economic equation, prior, data in the running pilots, or diagnostic gate.

## Analytically established

For the pilot flags (fixed df=5, rho/HMT/pooling/market residual scales off),
each share observation has a marginal Student-t density with scale
kappa*sigma_s, kappa=sqrt(3/5). The baseline prior is
log(sigma_s) ~ Normal(mu,sqrt(v)), mu=log(.5), variance v=.25. Margin noise also has a proper
lognormal scale prior. A Student-t density is bounded by C/sigma independently
of its location. Consequently the full finite hard-game mixture likelihood is
bounded by C_s^N C_m^M sigma_s^(-N) sigma_m^(-M). Lognormal priors have finite
inverse moments of every finite order; the other active priors are normalized.
The posterior normalizing constant is therefore finite under these pilot
flags. Exact shares do not make this posterior improper.

They can nevertheless produce narrow concentration. For simultaneous or
Stackelberg supply, the Logit and CES log-share equations have an exact linear
representation y=X theta at fixed observed prices. The outside-price CES term
is absorbed by the period intercepts. There are T independent market-year
intercepts, G-1 firm quality contrasts and one price coefficient. The audited
within-period design has rank G=5, so the full rank is r=T+G.

Near an exact-fit point theta_0 with continuous positive prior density, write
theta=theta_0+sigma_s*u. The share likelihood and integration volume contribute
sigma_s^(-N) and sigma_s^r respectively. The integral over u is finite: select
r independent rows of X, integrate their Student-t densities under an invertible
linear map, and bound the other factors by their finite maxima. Thus the local
small-scale contribution is proportional to sigma_s^(r-N), up to continuous
prior/margin factors. In the log-scale coordinate l=log(sigma_s), its leading
log density is

    -(N-r)*l - (l-mu)^2/(2*v) + constant.

The corresponding local asymptotic center is mu-(N-r)*v. This is neither an
exact full-posterior mode nor a causal diagnosis of the observed HMC failures.
Other regions, supply-margin tradeoffs and weak directions can matter. The
fringe quality-shift column requires its own rank analysis; the table below
does not predict fringe-model scales.

| products | periods | N | full rank r | N-r | local asymptotic scale |
|---:|---:|---:|---:|---:|---:|
| 5 | 8 | 40 | 13 | 27 | 5.85e-4 |
| 5 | 32 | 160 | 37 | 123 | 2.21e-14 |
| 5 | 128 | 640 | 133 | 507 | 4.49e-56 |
| 6 | 8 | 48 | 13 | 35 | 7.92e-5 |
| 6 | 32 | 192 | 37 | 155 | 7.41e-18 |
| 6 | 128 | 768 | 133 | 635 | 5.68e-70 |

Repeating exact-share periods therefore need not ease computation. These
limiting scales are analytic illustrations; numerical floating-point residuals
can intervene well before the scales in the longer-panel rows.

## Verification and experiment decision

The independent source, full-rank, Student-t profile and curvature audit
passes. All four source equations and all six demand/layout designs are
checked. The T=8/32/128 dimension sweep confirms rank T+5. All 148 local profile
points have finite targets, improved local modes and positive-definite Hessians;
share-block log-scale slopes agree with -(N-r). This numerical mechanism test
uses simple orthonormal Normal nuisance priors, rather than the full native
prior. The full local tables and independent report remain local in
`exact_share_audit`. The initial minimization error is retained as failed
exploratory evidence and excluded from the passing checks.

Saved initial and dense cohorts each contain ten completed diagnostic failures.
Their scale summaries are exploratory observations from failed samplers, not
valid posterior estimates. The simultaneous median scales are approximately
6.07e-4 and 6.03e-4 respectively, consistent in order with the local calculation;
that comparison establishes no causal explanation or statistical recovery. Its final executed tables will be linked here.
The larger-data helper passed independent default/prefix/truth checks, but its
sampling screen remains on hold pending this diagnosis. The unchanged
eight-period orthonormal-quality pilot remains the next computational
experiment after its full package check passes. Neither mathematical
propriety nor unit-test success establishes mixing or statistical recovery.
