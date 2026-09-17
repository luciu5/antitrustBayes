# Market/outside-share rotation: evidence to date

Updated 2026-09-16. The experimental API retains `legacy` as its default.
This report distinguishes exact target validation from sampler performance.

## Analytically established

The approved rotation preserves the existing demand, supply, likelihood,
priors, and pooling assumptions. Priors are evaluated on the reconstructed
legacy coordinates, with the Jacobian for the actual sampling coordinates.
See `market_s0_rotation_math.md` for bounded outside-share, noncentered
hierarchy, and Cournot details. Economic posterior correlations need not
decrease under an exact coordinate transformation.

## Numerically verified

The independent native gate passed 468 supported fixed-parameter cases across
Logit/CES, Bertrand/Cournot, supported outside-share structures, sampled and
collapsed quality, exact ordered-k structures, and supported rho modes.
The largest normalized target difference in the main 176-case comparison was
2.85e-14; the largest automatic-gradient versus finite-difference difference
was 6.87e-9. These checks retain density constants and transformation
Jacobians. See `market_s0_rotation_native_validation.md`.

## Sampler evidence

The frozen primary fixture has two markets, two periods, 96 observations,
and 48 firm groups. Both pairs use collapsed Gaussian shares, Gaussian
margins, simultaneous Logit Bertrand, unpooled market-year outside shares,
common alpha, rho off, and identical economic initial states, priors, and
sampler budgets. Four chains run sequentially, with 1,000 iterations each,
500 warmup, adapt_delta 0.95, and maximum treedepth 13.
The two seeds change sampler randomness on the **same dataset**; they are
not independent data-generating-process replications.

| Sampler seed | Legacy wall seconds | Rotated wall seconds | Wall speedup | Structural ESS/second gain | Divergences, legacy / rotated |
| --- | ---: | ---: | ---: | ---: | ---: |
| 271828 | 548.3 | 43.0 | 12.76x | 9.37x | 3 / 0 |
| 314159 | 528.3 | 42.6 | 12.40x | 11.10x | 0 / 0 |

Both rotated fits pass the unchanged prespecified gate: zero divergences and
treedepth hits, BFMI at least 0.3, maximum rank R-hat at most 1.01, and minimum
bulk/tail ESS at least 100. Legacy seed 314159 also passes. Legacy seed 271828
fails because of its three divergences and remains in the results.
ESS/second uses the runner's minimum monitored structural bulk ESS divided
by native Stan time, including warmup. Wall and native time are recorded
separately; summed native chain times can differ from wall time.

The reconstructed economic market-intercept/outside-logit correlation remains
approximately -0.997 in both parameterizations. Maximum absolute correlation
between the new raw market-intercept and rotated outside-share coordinates
is 0.048 and 0.086 in the two rotated fits. This is consistent with improved
sampling geometry while retaining the economic ridge.

Compact results are in `market_s0_rotation_hmc_benchmark.csv` and
`market_s0_rotation_hmc_efficiency.csv`. Full fits, warnings, starts, and
failed attempts remain under the ignored Windows-drive archive
`market_s0_rotation_runs/20260916/paired-hmc-attempt003/`.

## Scaling evidence

Both rotated primaries passed, permitting the frozen scaling experiment.
All five initial scaling attempts (legacy M=10, rotated M=10/100/500/2000)
stopped **before HMC**: a length-one parameter declared as a vector was
supplied without its dimension attribute (`mu_b_market_vec`). Those attempts
establish neither scaling success nor scaling sampler failure. Their logs and
exit statuses are preserved. A correction limited to initialization shapes
passed private reference validation; no numeric starts, priors, support, or
HMC settings changed. The same audit exposed nonfinite values in the frozen
rotated scaling starts. It preserved those values and rejected the starts.
Scaling remains deferred pending a separately recorded fixture-generation
correction; no result from those failed starts is called a sampler result.

## Posterior, policy, and bridge equivalence

Native target equality has passed. Empirical posterior, repeated bridge,
and merger/cost-shock equivalence for the usable seed-314159 pair remain
pending. Seed 271828's legacy fit cannot establish empirical equivalence
under the prespecified gate. No completed policy comparison from another
experiment is substituted for this coordinate experiment.

## Recommendation and unresolved issues

The measured speedup is substantial for this tested configuration. Keep
rotation experimental until scaling and usable-pair posterior, bridge, and
policy comparisons are complete. The primary benchmark does not establish
performance for the automatic noncentered market-level outside-share pool,
Student-t likelihoods, Cournot, CES, or estimated rho; fixed-target validation
of these paths alone is not an HMC performance result.
