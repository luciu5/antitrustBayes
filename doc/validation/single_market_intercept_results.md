# Single-market Logit utility correction, 2026-10-21

Starting commit: `55c6fdc53739d29b8be98f2e943a822c8f3204e3` on `main`.
Correction branch: `fix/single-market-intercept`. Package `0.2.0.9007`, fitted
model specification 12. The final correction commit is recorded in Git history.

## Cause and correction

Both Logit kernels formerly omitted `mu_b_market_vec` for one market and fixed
the inside-market intercept to zero. With centered qualities and prices, the
noiseless seven-firm M1T1 share equation then imposed
`s0 = plogis(mean(log(r))) <= 1/8`. Positive measurement noise could absorb the
incorrect common level; it did not justify that economic restriction.

Both kernels now sample one `mu_b_market_vec[1]` directly for a single market,
with exactly one normalized `Normal(0,1)` prior. The single-market hierarchy SD
and raw deviations remain absent. Multiple-market hierarchy, outside utility
zero, quality centering, price standardization, likelihoods, conduct equations,
outside-share pooling, priors on other parameters and counterfactual rules are
unchanged. `is_single_market` retains its original meaning. Reporting fields
`b_market`, `market_intercept_location`, `gamma_market` and the outside-share
coordinate sums now carry the intercept consistently.

`R/helpers.R` resolves a single-market collapsed-coordinate request to the
sampled intercept and advances the specification guard to 12. `R/estimate.R`
records the active population/intercept prior. `R/priors.R`, help and the
[normalization contract](../single_market_utility_normalization.md) explain
the prior and identification. Existing fits must be refitted; no legacy
normalization option was added. The two CES kernels are unchanged.

`bayes_calibrated_market(..., outside_share = .20)` accepts an independent
baseline outside share, default .20. `R/calibrated_market.R` retains centered
log-share contrasts separately from the baseline utility index, includes its
intercept once and permits the outside share to respond to subsequent shocks.

## Deterministic and package checks

The ordinary tests independently construct valid seven-firm Logit demand at
outside shares .10, .20 and .40 with identical conditional shares. Nine explicit
positive-cost B/C/CF equilibria also verify all three outside-share levels with
unchanged primitive costs and conditional shares, independent demand identities,
physical-cost recovery and fresh equilibrium checks. A free
intercept satisfies all three equilibria to floating-point precision. MC
margin identities and bounded-link prior-sensitivity quadrature demonstrate
its weak outside-share identification. A heterogeneous Cournot markup line
independently identifies alpha and s0 in noiseless data; no requirement that
every model tightly recover s0 was imposed.

All eight conduct families are checked across M1T1, M1T4, M4T1 and M4T4, with
independent, market and market-year outside-share pooling. Dedicated downstream
tests pass 175 assertions covering nonzero-intercept quality/predictive means,
an independent dense Gaussian conditioning oracle, B/C/CF/CS demand, margins,
costs and equilibrium, ordered-k roles, merger and tariff equivalence. The
existing R consumers already use `share_mean`, so no additional intercept or
economics branch was introduced into posterior simulation.

| Tier | Passed assertions | Failures/errors | Test warnings | Skips |
| --- | ---: | ---: | ---: | ---: |
| Complete current development-source ordinary suite | 6,454 | 0 / 0 | 22 | 12 opt-in |
| Existing opt-in CmdStan filter | 590 | 0 / 0 | 9 | 0 |
| Installed source-package tests | 6,219 | 0 / 0 | 22 | 29 |

`R CMD build` and `R CMD check --no-manual` succeed with zero errors, zero
warnings and one existing NOTE for unqualified `tail` in `.cmdstan_assert_fit`.
Four additional solver `termcd: 2` warnings in the nine new equilibrium checks
are retained: all independent demand, cost and fresh-equilibrium assertions pass.
Installed-package load, all four byte-identical Stan sources, new-intercept
extraction, posterior policy plans and help parsing also pass. The source
snapshot preserves raw source bytes instead of applying Git line-ending
conversion. Native compilation took 72.92 seconds for Logit supply and 75.42
seconds for the ordered-k kernel through CmdStan; BridgeStan compilation is
recorded separately in the native outputs.

Current-specification test fixtures were advanced to 12 and an explicit
specification-11 refit rejection added. A pre-existing source-hash corruption
test assumed the first hash character was nonzero; the new source starts with
zero. Its mutation now always changes that character, retaining all mismatch
assertions. No scientific assertion or production CSV precision was weakened.

## Fresh native target and evidence protocol

The [bounded runner](single_market_intercept_validation.R) seals six new
fixtures: Cournot in all four market/time layouts and fixed-k3 CF/CS in M1T1.
Each backend uses four chains, 500 warmup plus 500 retained states per chain,
two cores, acceptance .95 and depth13, with fixed seeds and no retuning. The
panel fixtures deliberately use per-market-year outside shares; they do not
validate a production restriction of constant outside share across time.

The protocol comprises twelve fresh fits, six 23-state target comparisons,
six public three-repetition bridge calculations and an auxiliary MC prior
check. No old posterior draws are reused. A sealed old compiled target is used
only as a mathematical reference: with identical MC margin-only data and all
other coordinates matched, the absolute new-minus-old normalized log density
must equal `dnorm(a, log = TRUE)`, including its constant. The reference does
not provide a legacy estimation option.

All **25/25 tasks completed**, with no failed execution or replacement fits.
All 27 fresh equilibrium gates pass (maximum price/share/cost/FOC errors
`5.45e-9`/`3.69e-9`/`1.81e-14`/`1.38e-8`). Ten construction solver warnings
are preserved alongside those independent checks.
The 138 paired target states agree with maximum normalized log-density error
`8.19e-12`. All 425 posterior mean/median/SD/interval comparisons with defined
MCSE agree within four combined MCSE; the other 90 fixed-k role/component
constants agree exactly and keep their unavailable MCSE labels. The largest
standardized posterior difference is `2.60643` combined MCSE. Four-draw
posterior quality reconstruction passes for every pair.

The five-point auxiliary absolute density difference matches the normalized
`Normal(0,1)` log density to `2.00e-15`, including `-.5*log(2*pi)`. This directly
checks the new factor once, rather than canceling its constant in a density
ratio. Both libraries and the old mathematical reference are hash-sealed.

All twelve saved fits also pass a separate generated-quantity audit of
`mu_b_market_vec`, `market_intercept_location`, `b_market`, `gamma_market` and
`c_effective_mky`. No single-market hierarchy SD/raw columns are present.
The maximum CmdStan serialization discrepancy is `1.20e-7`, within the
explicit `1e-6` tolerance for its unchanged eight-significant-digit CSV;
Nutpie identities pass at `1e-10`. See [generated quantities](single_market_intercept_generated_quantities.csv)
and [source/data/prior provenance](single_market_intercept_provenance.csv).

Actual diagnostics from the fixed short budget are:

| Fixture | Backend | Divergences | Max R-hat | Min bulk ESS | Min tail ESS |
| --- | --- | ---: | ---: | ---: | ---: |
| C_M1T1 | CmdStan | 12 | 1.0109 | 321.8 | 492.6 |
| C_M1T1 | Nutpie diagonal | 1 | 1.0129 | 214.5 | 239.5 |
| C_M1T4 | CmdStan | 0 | 1.0063 | 379.2 | 409.6 |
| C_M1T4 | Nutpie diagonal | 0 | 1.0171 | 113.4 | 85.9 |
| C_M4T1 | CmdStan | 10 | 1.0065 | 509.6 | 654.6 |
| C_M4T1 | Nutpie diagonal | 1 | 1.0361 | 201.8 | 67.4 |
| C_M4T4 | CmdStan | 2 | 1.0109 | 666.1 | 846.1 |
| C_M4T4 | Nutpie diagonal | 1 | 1.0142 | 378.0 | 483.7 |
| CF_M1T1 | CmdStan | 17 | 1.0133 | 399.6 | 256.7 |
| CF_M1T1 | Nutpie diagonal | 2 | 1.0138 | 332.8 | 200.3 |
| CS_M1T1 | CmdStan | 4 | 1.0095 | 296.6 | 331.6 |
| CS_M1T1 | Nutpie diagonal | 0 | 1.0458 | 114.1 | 126.2 |

All fits have zero depth hits. There are 45 CmdStan and five Nutpie divergences
in total. Maximum R-hat `1.0458` and minimum tail ESS `67.4` retain important
mixing limitations. Early rejected proposals with underflowed scale parameters
and repeated loading of the same hash-verified BridgeStan library are retained
in the native logs/status records. None triggered retuning or scientific
candidate removal.

All **18/18 requested bridge repetitions are finite**:

| Fixture | CmdStan log Z | Nutpie log Z | Nutpie minus CmdStan | Four combined errors |
| --- | ---: | ---: | ---: | ---: |
| C M1T1 | 3.125000 | 3.191474 | +.066475 | .259282 |
| CF M1T1, fixed k=3 | .693413 | .722021 | +.028608 | .287106 |
| CS M1T1, fixed k=3 | .232991 | .246425 | +.013434 | .276369 |

All three evidence differences agree within their bridge Monte Carlo bounds.
All three remain **precision-inconclusive** against the prespecified `.25`
four-error cap. The short-run evidence checks therefore support target/prior
consistency and finite normalized calculation, rather than tight Bayes-factor
precision or broad sampler qualification. Conditional-k arithmetic, full
positive-prior support and no-survivor-renormalization tests remain unchanged;
this task does not add another full conditional-k sampling grid.

Independent review verified all 25 result/manifest bindings, all twelve saved
fit seals, 39 frozen runtime/source hashes and five historical-reference seals,
and cleared the scoped correction for merge. This does not certify precise
evidence, broad recovery, or a full new-target merger/tariff distribution
benchmark. Deterministic economic/policy equivalence is the counterfactual
validation executed for this correction.

Exact estimates, numerical uncertainties and statuses are in
[the generated report](single_market_intercept_validation_results.md),
[evidence](single_market_intercept_evidence.csv),
[paired evidence comparisons](single_market_intercept_evidence_comparisons.csv),
[prior check](single_market_intercept_prior.csv), and the adjacent posterior,
target and diagnostic CSVs. Sampler diagnostics describe reliability and never
suppress available evidence. Historical specification-11 benchmark claims are
not used as statistical recovery or speed evidence for this corrected target.

## Revised mergerBayes baseline integration

The unchanged archived revised DGP has true s0 .20, while its conditional
shares imply .0908512008087206 under the old zero-intercept restriction. The
required intercept is -.916991402898777. Corrected noiseless demand errors are
at most 1.33e-15. The public calibrated-market API reproduces its C/CF/CS
baseline prices, shares, margins, costs, alpha and mean values with exactly
zero differences. See [baseline identities](single_market_intercept_integration.csv).

A fresh Nutpie C M1T1 fit uses exactly the archived seven measured rows and
all seven observed margins, with four chains and 500 warmup/500 retained
states per chain. Outside-share mean is .180977, median .175548, 95% interval
[.125923,.264893], and mean MCSE .001904. It no longer faces the old geometric
restriction. See [posterior summaries](single_market_intercept_integration_posterior.csv).
There are zero divergences/depth hits, BFMI .830--.978, maximum R-hat 1.02206,
minimum bulk ESS 381.92 and tail ESS 470.92. This is a short integration smoke,
with visible mixing limitations, rather than precise recovery certification.

No estimator workaround is required in mergerBayes after installing the updated
package and rerunning fits. Its existing private DGP adapter can later use the
new public argument; this correction did not edit downstream code. Later-period
equilibrium shocks vary true outside shares (the revised report records
18.39--21.26%), whereas some estimation calls pool s0 at market level. That
separate time-pooling assumption remains unchanged and requires its own review.

All five inspected downstream source/master/manifest seals are unchanged. No
mergerBayes production analysis, DGP parameter, concentration, margin, cost
distribution, sampler profile, selected run or scientific artifact was changed.
No 64-cell experiment or Slurm job was launched.

## Durable raw evidence

Scripts, runtime and source seals, exact fixtures, fit checkpoints, full native
diagnostics and warnings are under ignored `.scratch/single-market-intercept/`.
The directory also contains `prechange-diagnostic.json`, `ordinary-final-candidate.*`,
`cmdstan-regressions.*`, `package-check-final.log`, the exact source package/check,
`installed-smoke.*`, and `merger-integration/`. Raw fits and native libraries are
not committed. Historical Nutpie benchmarks describe specification 11 and
must not be read as performance measurements of this corrected target.
