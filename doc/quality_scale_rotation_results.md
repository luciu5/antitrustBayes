# Collapsed joint-Gaussian geometry results

## Question

The exact joint market/outside-share rotation removed the original
market-\(s_0\) cancellation but left 12 and 9 divergences in two frozen
Supreme Court year-cross-section fits. This experiment asks whether the
remaining geometry comes from the quality/share-noise variance decomposition
and its interaction with active \(\rho\).

## Analytically established

The `total_fraction` coordinates

\[
T=\sqrt{\tau^2+\sigma_s^2},\qquad
\lambda=\sigma_s^2/T^2
\]

are an exact reparameterization. With sampler coordinates
\(x=\log T\) and \(y=\operatorname{logit}\lambda\), the Jacobian under the
historical measure \(d\log\sigma_s\,d\tau\) is \(\tau/2\). The implementation
evaluates the existing normalized priors at reconstructed economic scales and
includes this Jacobian exactly.

The margin-first collapsed Gaussian factorization is also exact. It replaces
the historical conditional-margin coefficient
\(\rho\sigma_m/\sigma_s\) with the bounded-at-small-\(\sigma_s\) conditional
share coefficient \(\rho\sigma_s/\sigma_m\). It does not alter the joint
covariance, likelihood, prior, conduct equation, or economic output.

## Numerically verified

* Historical conditional collapsed kernel versus dense joint covariance:
  120/120 likelihood checks and 72/72 gradient checks passed.
* Total/fraction R reference: 5/5 target/reconstruction checks and 10/10
  chain-rule gradient checks passed.
* Native Stan coordinate check before the margin-first change: 170/170 target
  and generated-quantity comparisons and 8/8 autodiff checks passed across
  Logit/CES and simultaneous/exact-\(k\) kernels.
* Margin-first collapsed kernel versus dense joint covariance: 120/120
  likelihood checks and 72/72 gradient checks passed.
* Conditional versus margin-first native Stan generated quantities: 1,760/1,760
  elementwise checks passed across Logit/CES and simultaneous/exact-\(k\)
  kernels; the maximum absolute difference was \(1.02\times10^{-10}\).
* The 25 shared economic posterior summaries in the paired sparse synthetic
  fits agreed within the prespecified Monte Carlo tolerance.

## HMC evidence before the margin-first change

All runs used the joint market-\(s_0\) coordinate and active \(\rho\).

| Fixture | Scale coordinates | Seed | Draws | Divergences | Treedepth | Max R-hat | Min bulk ESS | Min tail ESS | Seconds |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 6-market sparse synthetic | total/fraction | 19091 | 1,600 | 0 | 0 | 1.0081 | 605 | 433 | 47 |
| 6-market sparse synthetic | total/fraction | 29091 | 1,600 | 0 | 0 | 1.0064 | 514 | 405 | 45 |
| 6-market, 3-period synthetic panel | total/fraction | 19091 | 1,600 | 0 | 0 | 1.0056 | 329 | 419 | 550 |
| Supreme Court cross-section | legacy | 19091 | 4,000 | 12 | 0 | 1.0041 | 559 | 492 | 1,954 |
| Supreme Court cross-section | total/fraction | 19091 | 4,000 | 6 | 0 | 1.0091 | 582 | 829 | 2,187 |

In this one paired seed, total/fraction coordinates reduced the exact-data
divergence count from 12 to six, missed the zero-divergence gate, and were about
12 percent slower. This does not establish a consistent improvement across
seeds or a causal explanation for the remaining divergences.

## Margin-first HMC evidence

| Fixture | Seed | Draws | Divergences | Treedepth | Max R-hat | Min bulk ESS | Min tail ESS | Seconds | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 6-market sparse synthetic | 19091 | 1,600 | 0 | 0 | 1.0060 | 677 | 364 | 47 | pass |
| 6-market sparse synthetic | 29091 | 1,600 | 0 | 0 | 1.0068 | 563 | 531 | 47 | pass |
| 6-market, 3-period synthetic panel | 19091 | 1,600 | 0 | 0 | 1.0125 | 383 | 474 | 560 | fail: R-hat |
| Supreme Court cross-section | 19091 | 4,000 | 10 | 0 | 1.0021 | 713 | 1,164 | 1,742 | fail: divergences |
| Supreme Court cross-section | 29091 | 4,000 | 6 | 0 | 1.0037 | 598 | 1,046 | 2,678 | fail: divergences |

The panel result is retained as a failure under the prespecified maximum
\(\hat R\le1.01\) rule even though it had no divergences or treedepth hits.
All 27 shared economic posterior summaries in the paired sparse conditional
and margin-first fits agreed within the Monte Carlo tolerance.

On the frozen 176-observation Supreme Court fixture, margin-first evaluation
had 10 and six divergences in two prespecified seeds. The matched
conditional-margin total/fraction run had six divergences in the first seed.
All runs had zero treedepth hits and otherwise healthy diagnostics. The
different factorization therefore does not resolve the remaining sampler
geometry. The second margin-first seed selected a smaller mean step size
(0.0050) and took 2,678 seconds, so the first seed's 20 percent speed advantage
does not establish a general runtime improvement.

The original resumable runner preserved the primary HMC row but overwrote the
seed-19091 auxiliary posterior/divergence-detail CSV when seed 29091 resumed.
The headline seed-19091 metrics and divergence medians had already been
recorded above. The runner now resumes all auxiliary CSVs by appending existing
rows; a no-sampling resume check over both completed seeds exits cleanly.

The 10 divergent draws were associated with a lower share-variance fraction
(median 0.030 versus 0.081), a lower residual share scale (0.263 versus 0.436),
and a somewhat higher quality scale (1.544 versus 1.443). These are descriptive
associations from rare events, not a causal localization.

For collapsed joint exact-\(k\) fits, normalized market-level component scores
and posterior component probabilities are produced and were included in the
native equivalence checks. Per-observation margin log-likelihood entries are
not a valid decomposition after collapsing persistent quality and remain
unset on that path; pointwise LOO is consequently disabled for collapsed fits.

## Decision

The old conditional kernel is algebraically correct. Total/fraction coordinates
reduced divergences in one paired seed but have not established a general
improvement. Reversing the
conditional factorization is numerically exact and removes an explicit
division by the residual share scale, but it does not fix the posterior
geometry. At identical coordinates, the exact joint density and its mathematical
curvature are unchanged by the choice of factorization; numerical evaluation
and resulting sampler adaptation can still differ. Neither coordinate change becomes a package default on this
evidence, and the upper-ladder production run remains unauthorized.
