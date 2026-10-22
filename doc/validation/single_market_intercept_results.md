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

## 2026-10-21 follow-up: revised-paper Bertrand M1T1 cases

Both authorized revised-paper Bertrand-truth cases, `share04_B_rep101_M1T1`
and `share08_B_rep101_M1T1`, completed under specification 12 at commit
`112c4480ff087e595dba10170a996b8d7ce28fe1` (version `0.2.0.9007`). These
are **CmdStanR-only** results for one fixed rep101 economy at two measurement
noise levels. They do not establish repeated-sample recovery or Nutpie speed.
The complete archived master/economy hashes and policy-target objects match
the revised `B_share04`/`B_share08` masters exactly, including the seven
observations, margin masks, fixed merger, and 25% tariff definitions.

The unchanged production profile fits MC/B/C/A and all six fixed-k components
of each BF/CF/BS/CS family: **28 fits per case, 56 total**. MC and A each have
prior probability 1/4; the other six candidates each have 1/12. Every rich
family retains its full original k=1:6 support and uniform k prior, restored
exactly once through the existing conditional-k evidence calculation.
All fits use four chains and at most four sampling cores. MC/B/A retain
2000 warmup and 2000 posterior draws per chain with `.99` acceptance,
depth 14 and diagonal adaptation; C retains 4000/4000 with `.999`, depth 14
and dense adaptation. Each rich-family component retains 2000/2000 with
`.999`, depth 14 and dense adaptation. Seeds, raw inputs, exact calls,
continuous priors, source/data hashes, and complete native chains were
independently checked against the sealed manifest. All 56 resolved
continuous-prior metadata hashes agree.

All **168/168 requested normal bridge repetitions are finite** (84 per case;
three per fitted atomic/component target, `maxiter=5000`). Independent
log-sum-exp calculations reproduce every family evidence, posterior k weight,
and eight-model probability. Exact values and repetition SDs are in
[model comparison](single_market_intercept_bertrand_models.csv), with
[family probabilities](single_market_intercept_bertrand_families.csv).

| Scenario | Bertrand-family prior | B/BF/BS posterior | Atomic B posterior | MAP candidate |
| --- | ---: | ---: | ---: | --- |
| share04 | 25% | 60.067% | 20.631% | A, 21.379% |
| share08 | 25% | 59.527% | 20.460% | A, 21.584% |

The evidence increases total probability on the Bertrand family; it does not
select atomic Bertrand as the most probable individual candidate.

| Scenario | Divergences | Depth hits | Max R-hat | Min bulk ESS | Min tail ESS | End-to-end elapsed |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| share04 | 2 | 0 | 1.003512 | 1649.98 | 824.70 | 39m24s |
| share08 | 3 | 0 | 1.005432 | 1881.54 | 850.06 | 36m58s |

Elapsed times are controller start-to-finish times including worker startup,
compilation/cache use, evidence and policies, not sampling-only timings.
[Fit diagnostics and provenance](single_market_intercept_bertrand_fits.csv)
retain all 56 source/data/prior/checkpoint hashes, actual controls, BFMI and
per-fit elapsed times. Divergences remain disclosed. Raw execution logs also
retain repeated unchanged-ownership/Logit-share adapter notices. Raw CmdStan
CSVs were kept when the historical cleanup comparator could not confirm
serialized function/source-reference object identity; independent exact
source, data, prior, call and native-completion checks passed.

Merger and tariff use the same posterior draw, k, quality and plan identities
within every model. Each atomic candidate requests 200 states; each rich
family requests 200 per k component. Both cases produced all 130 existing
policy-summary rows. Selected market-basket price effects, in percent, are:

| Scenario | Policy | Oracle | BMA mean | BMA SD | BMA 95% interval | B-only mean | MAP mean |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: |
| share04 | Merger | 5.520 | 4.371 | 2.131 | [0, 8.087] | 5.315 | 3.100 |
| share08 | Merger | 5.520 | 4.480 | 2.143 | [0, 8.386] | 5.392 | 3.347 |
| share04 | Tariff 25% | 8.332 | 7.975 | .639 | [6.222, 8.940] | 8.272 | 7.867 |
| share08 | Tariff 25% | 8.332 | 7.989 | .627 | [6.270, 8.933] | 8.278 | 7.872 |

[Policy results](single_market_intercept_bertrand_policy.csv) retain all 20
`price_pct` rows: five existing procedures, both policies and both cases,
including threshold probabilities and explicit support/attrition fields.
The equal-model-prior and joint-support sensitivity procedures remain
separately labeled.

CF retains three finite nonpositive-cost states in share04 (k1 draws
2019/2495 and k3 draw 3227) and four in share08 (k1 draws 3424/3026/1413 and
k2 draw 5514). They are the same states in both policies. CF has 1197/1200
and 1196/1200 valid states respectively; all other candidates retain every
requested state. Primary BMA summaries condition on successful draws within
each model while retaining the original model probabilities, with original
weighted negative-cost mass `9.853273e-6` and `6.286872e-6`. These omissions
and the separate joint-support treatment remain visible; no candidate or k
state was removed from evidence or model weighting.

Raw fits, policies, warnings, manifests and controller receipts remain under
ignored `.scratch/single-market-intercept/paper-bertrand-m1t1-20261021-v2/`;
the earlier failed preflight-only attempt is preserved separately. No
mergerBayes production source, selected analysis, DGP, or scientific/sampler
setting was changed. Hosted ordinary CI passed at `112c448`
([run 39710993394](https://github.com/luciu5/antitrustBayes/actions/runs/39710993394)).

## 2026-10-21 follow-up: PNB 1960 M1T1

The requested PNB cross section completed locally with fresh specification-12
CmdStanR fits at `112c4480ff087e595dba10170a996b8d7ce28fe1`. The unchanged
audited banking inputs contain 39 banks and 32 observed margins; the seven
banks with missing margins remain in demand. The fixed merger joins target
`889018` and acquirer `298816`. No Slurm job, banking-data rebuild, DGP change,
or selected production-result replacement was performed.

The seven-model universe is MC/B/C/BF/CF/BS/CS, with MC prior 1/3 and each
other candidate 1/9. Thus MC, the Bertrand family and the Cournot family each
start at 1/3. Auction is excluded by the existing banking design. Every rich
fit uses the original native marginalized **uniform k=1:5** support, rather
than the synthetic experiment's conditional-k decomposition over k=1:6.

All seven fits retain the maintained banking controls: four chains, 1000
warmup and 1000 retained draws per chain, seed 19091, initialization zero,
diagonal metric, acceptance .99 and maximum depth 13. The likelihood is joint
Gaussian with active correlation, input-market equations, collapsed market
quality and sampled intercept. `s0_pooling="auto"` resolves to `"none"` in
this cross section. No sampler control was tuned after inspecting outcomes.
The existing paper script still pins specification 11 in its loader and
sidecar contracts. The local research wrapper extracts its unchanged function
definitions and fitting-argument expressions, records specification 12
explicitly, and avoids those obsolete guards without impersonating version
11. Updating the maintained downstream runner remains an integration task
before production use; no mergerBayes source was edited here.

All **21/21 requested normal bridge repetitions are finite** (three per
model, `maxiter=5000`). The finalizer requires all seven finite evidences;
independent probability reconstruction applies the original priors once.
Complete values, bridge SDs, diagnostics and source/data/prior hashes are in
[model comparison](single_market_intercept_pnb_models.csv).

| Model | Prior | Posterior | Log marginal likelihood | Bridge repetition SD | Divergences |
| --- | ---: | ---: | ---: | ---: | ---: |
| MC | 33.333% | 40.621% | -143.9828 | .004724 | 0 |
| B | 11.111% | 12.527% | -144.0606 | .005029 | 0 |
| C | 11.111% | 6.393% | -144.7332 | .006172 | 2 |
| BF | 11.111% | 12.775% | -144.0410 | .003791 | 1 |
| CF | 11.111% | 6.706% | -144.6855 | .006368 | 3 |
| BS | 11.111% | 12.644% | -144.0513 | .007646 | 2 |
| CS | 11.111% | 8.334% | -144.4682 | .011890 | 1 |

[Family probabilities](single_market_intercept_pnb_families.csv) are
**40.621% MC, 37.946% Bertrand and 21.433% Cournot**, from equal one-third
priors. These are modest evidence differences, rather than decisive conduct
identification. Diagnostics retain nine divergences, zero depth hits, maximum
R-hat 1.004457, minimum bulk ESS 1159.10, minimum tail ESS 884.42 and minimum
BFMI .84631. All seven fits and their evidence remain available for review.
The controller ran 23:05:56–23:11:55 EDT (about six minutes), sequentially after
the synthetic runs, using at most four sampling cores and their original
shared absolute deadline. This elapsed time includes startup, evidence and
policy calculations, not just HMC.

The saved posterior directly retains the new nonzero inside-market intercept.
For atomic B, its mean is -2.9631 with 95% interval [-3.8677, -2.0663]; the
outside-share mean is .12039 with interval [.05334, .22918]. The complete
[structural summaries](single_market_intercept_pnb_structural.csv) include
the intercept, outside share and physical price sensitivity for every model;
the original banking [parameter summaries](single_market_intercept_pnb_parameters.csv)
retain noise and correlation estimates. These empirical data have no known
true outside share, and freeing the intercept does not guarantee precise
identification of outside share.

All seven models retain 200/200 valid merger states, with no reported policy
attrition. The maintained banking summary pools fractional deposit-rate
changes across **banks and posterior draws, with equal bank weights**. Its
mean is -1.285% for BMA, -1.893% for B alone and zero for MC; the corresponding
pooled bank/draw SDs are 9.305, 7.535 and zero percentage points. These SDs
are not posterior SDs of a market-basket aggregate. The maintained helper
does not retain an aggregate posterior interval, so none is inferred here.
[Policy summaries](single_market_intercept_pnb_policy.csv) preserve the exact
fractional-unit results and explicitly labeled percent conversions. Tariffs
are not part of this banking route.

The active continuous priors match the synthetic M1T1 fits: in particular,
`logit(s0) ~ Normal(-1.446786403764666, .456234014550216)` (outside-share prior
mean .20, central 95% interval [.08778, .36527]), one `Normal(0,1)` inside
intercept, and identical standardized price, quality, noise and correlation
priors. Comparing the full resolved manifests finds only one differing
entry: the market-intercept hierarchy SD family is Half-Normal(.5) in PNB
and Gamma(shape=2, rate=5) in synthetic. That parameter is absent with one
market. Identical standardized price priors imply different physical-unit
priors when the empirical price SD differs. Model weights, k support and
input/output-market economics therefore need explicit labels when comparing
the two experiments.

Sealed inputs, attempts, raw fits, native CSVs, bridge objects, diagnostic
arrays, policy summaries and the reviewed wrappers remain under ignored
`.scratch/single-market-intercept/pnb-m1t1-20261021/`. `summarize.R` performs
only postprocessing and verifies the source/wrapper/data seals, all seven
model probabilities, all saved fit seals and identical active prior settings.
The pre-HMC contract tests passed 111 assertions in six blocks. No model was
dropped, no prior mass was renormalized over survivors and no further HMC was
launched during summarization.

## 2026-10-22 follow-up: matched Nutpie PNB M1T1 rerun

The authorized Nutpie rerun completed all seven PNB candidates under the same
corrected specification 12 at `112c448`. It reused the sealed observed input
from the preceding CmdStan run. Every fit has **identical full Stan data,
Stan-source/data hashes, resolved priors, model specification and model prior**
to its CmdStan counterpart. The original native uniform k=1:5 support remains
in all rich models. The only fitting changes are the sampler backend, explicit
Nutpie diagonal adaptation replacing CmdStan's `diag_e`, and removal of the
CmdStan-specific CSV output option. Common controls remain four chains, 1000
warmup plus 1000 retained draws per chain, seed 19091, initialization zero,
acceptance .99, maximum depth 13 and at most four sampling cores. NutpieR
1.9.1 and the BridgeStan 2.8.0 R client were selected from the existing library;
both matching kernel libraries were already cached. No dependency was installed.

All **21/21 normal bridge repetitions are finite**, using the same three-repeat
and `maxiter=5000` budgets. No diagnostic threshold gates model evidence or
probabilities. Every positive-prior candidate is retained. The maximum
absolute log-evidence difference is .02133; all seven differences are within
.720 of the conservative combined bridge uncertainty used here. That
uncertainty is the square root of the sum of the two maximum squared
within-repetition log-evidence MCSEs and the two squared repetition SDs; shared
posterior draws within a backend are not treated as independent repetitions.
The maximum individual model-probability difference is .576 percentage points.

| Model | Prior | CmdStan posterior | Nutpie posterior | CmdStan divergences | Nutpie divergences |
| --- | ---: | ---: | ---: | ---: | ---: |
| MC | 33.333% | 40.621% | 41.197% | 0 | 0 |
| B | 11.111% | 12.527% | 12.428% | 0 | 0 |
| C | 11.111% | 6.393% | 6.301% | 2 | 0 |
| BF | 11.111% | 12.775% | 12.723% | 1 | 1 |
| CF | 11.111% | 6.706% | 6.646% | 3 | 0 |
| BS | 11.111% | 12.644% | 12.292% | 2 | 1 |
| CS | 11.111% | 8.334% | 8.413% | 1 | 1 |

[Full comparison](single_market_intercept_pnb_nutpie_comparison.csv) includes
log evidence, MCSE, diagnostics, timings and exact source/data hashes;
[bridge repetitions](single_market_intercept_pnb_nutpie_bridge_repetitions.csv)
preserve all 42 estimates from the paired runs. Family posterior probabilities
are **41.197% MC, 37.443% Bertrand and 21.360% Cournot** under Nutpie, compared
with 40.621%, 37.946% and 21.433% under CmdStan, from equal one-third priors.
See the [family comparison](single_market_intercept_pnb_nutpie_families.csv).

The 42 paired means for outside share, intercept, physical price sensitivity,
share/margin error scales and correlation differ by at most 2.902 combined
mean MCSEs. [Posterior comparisons](single_market_intercept_pnb_nutpie_parameters.csv)
also retain both medians, SDs and central 95% intervals. This bounded result
supports agreement for this empirical cross section, rather than general
recovery or universal sampler superiority.

| Diagnostic | CmdStan | Nutpie diagonal |
| --- | ---: | ---: |
| Divergences, seven models combined | 9 | 3 |
| Depth hits | 0 | 0 |
| Maximum R-hat | 1.004457 | 1.016147 |
| Minimum bulk ESS | 1159.10 | 469.09 |
| Minimum tail ESS | 884.42 | 798.82 |
| Minimum BFMI | .84631 | .82919 |

Nutpie's weaker maximum R-hat/bulk ESS occurs in C's outside-share coordinate;
its intercept has R-hat 1.01355 and bulk ESS 499.84. These mixing warnings
remain visible despite fewer divergences. No tuning or additional sampling
was performed after seeing them.

The sum of CmdStan process elapsed times over seven fits is 115.27 seconds
(111.68 seconds summing each fit's longest chain warmup-plus-sampling time).
Nutpie reports 73.23 seconds of native warmup and sampling, including callbacks
and native extraction; its complete seven `bayes_est()` calls total 78.54
seconds. Nutpie's cache/compile checks total .268 seconds. These are sums of
parallel-chain wall times across models, not CPU totals; backend timing scopes
differ. The whole Nutpie run took 222.46 seconds including R preflight,
evidence and policies. Its single R process avoids the eight worker startups
in the earlier CmdStan orchestration, so the whole-run elapsed difference is
not a pure sampler benchmark. Lower elapsed sampling time and fewer
divergences here do not imply better ESS per second for every parameter.

The unchanged PNB merger calculation retained **200/200 states per model**,
with no reported attrition. BMA mean bank deposit-rate change is -1.27198%
under Nutpie versus -1.28521% under CmdStan, a difference of .01324 percentage
points. The pooled bank/draw SDs are 9.40463 and 9.30545 percentage points.
[Policy comparison](single_market_intercept_pnb_nutpie_policy.csv) preserves
the existing equal-bank-weight statistic and its units; these SDs are not
posterior SDs of an aggregate policy effect, and no aggregate interval or
policy MCSE is inferred from them. The samplers produce different posterior
states; their 200-state policy calculations do not compare identical states.

The frozen runner, manifest, exact paired-identity checks, raw fits, native
diagnostics, all bridge attempts, warnings and execution receipt remain in
ignored `.scratch/single-market-intercept/pnb-m1t1-nutpie-20261022/`. The
outer four-hour process-group timeout was not approached. Cached-library
reload warnings were retained and no binary/source identity changed. No
mergerBayes source, scientific specification, selected result or production
inference was modified. CmdStan remains the package default.

## 2026-10-22 follow-up: Nutpie PNB acceptance .999

The requested acceptance sensitivity rerun completed separately in
`.scratch/single-market-intercept/pnb-m1t1-nutpie-delta999-20261022/`.
This is **Nutpie diagonal**, not a new CmdStan run. A sealed comparison of all
seven argument lists confirms that the only fitting change from the preceding
Nutpie run is `stan_options$control$adapt_delta = .999` (translated to native
`target_accept=.999`). Full Stan data, source/data hashes, model specification
12, priors, uniform k=1:5 support, initialization, seed, chain/draw budgets,
depth limit, libraries, adaptation mode and policy definitions are unchanged.
The earlier runs and their evidence remain intact. No production sampler
default or scientific specification was changed.

All seven fits contain the requested four complete chains, each with 1000
retained draws after 1000 warmup iterations. All **21/21 requested normal
bridge repetitions are finite**. Every positive-prior candidate remains in
the comparison. The [three-run comparison](single_market_intercept_pnb_nutpie_delta999_comparison.csv)
retains the corrected CmdStan .99 reference, Nutpie .99 and Nutpie .999.
Columns labeled `Nutpie` in these new tables refer to **.999**; `Nutpie99`
refers to the preceding run. The [repetition table](single_market_intercept_pnb_nutpie_delta999_bridge_repetitions.csv)
preserves all 63 bridge estimates and explicitly records target acceptance.

| Diagnostic | CmdStan .99 | Nutpie .99 | Nutpie .999 |
| --- | ---: | ---: | ---: |
| Divergences, all seven models | 9 | 3 | 0 |
| Depth hits | 0 | 0 | 0 |
| Maximum R-hat | 1.004457 | 1.016147 | 1.007238 |
| Minimum bulk ESS | 1159.10 | 469.09 | 679.72 |
| Minimum tail ESS | 884.42 | 798.82 | 805.26 |
| Minimum BFMI | .84631 | .82919 | .83954 |

The C outside-share coordinate improves from R-hat 1.01615/bulk ESS 469.09
to 1.00598/859.85. Improvements are not uniform: CF's minimum bulk ESS is
679.72 versus 754.26 at .99, and BF's maximum R-hat is 1.00724 versus 1.00510.
Zero observed divergences in this one run does not establish their absence
in future runs or universal efficiency superiority.

| Family | Prior, every run | CmdStan .99 | Nutpie .99 | Nutpie .999 |
| --- | ---: | ---: | ---: | ---: |
| MC | 33.333% | 40.621% | 41.197% | 40.473% |
| Bertrand | 33.333% | 37.946% | 37.443% | 38.035% |
| Cournot | 33.333% | 21.433% | 21.360% | 21.492% |

[Family probabilities](single_market_intercept_pnb_nutpie_delta999_families.csv)
and the full seven-model table retain the original priors. Maximum individual
model-probability change from Nutpie .99 is .7245 percentage points. The
maximum absolute log-evidence change is .02016 versus .99 and .02660 versus
CmdStan, within .643 and .758 respectively of the same conservative combined
bridge uncertainty defined in the preceding section. The 42 structural means
change by at most 2.2784 combined mean MCSEs versus Nutpie .99 and 1.6902
versus CmdStan. [Parameter comparisons](single_market_intercept_pnb_nutpie_delta999_parameters.csv)
also retain all three means, SDs, medians and central intervals.

All seven merger calculations retain **200/200 valid states**. The unchanged
BMA mean bank deposit-rate change is -1.28198%, versus -1.27198% at Nutpie
.99 and -1.28521% under CmdStan. The .999-minus-.99 change is -.01000
percentage points. The .999 pooled bank/draw SD is 9.17446 percentage points;
it is not an aggregate posterior SD. The [policy comparison](single_market_intercept_pnb_nutpie_delta999_policy.csv)
retains the same equal-bank-weight estimand and all successful-state counts.
No aggregate policy interval or MCSE is inferred from those pooled moments.

Nutpie native warmup-plus-sampling time, summed over seven fits, increases
from 73.226 to **111.878 seconds**, or 52.8%. Complete fitting calls total
117.248 seconds; the whole run, including preflight, evidence and policies,
takes 258.125 seconds. Backend timing scopes remain as described above;
the earlier CmdStan reference was not rerun at .999. No additional tuning or
sampling followed inspection of these outcomes.

Independent review checked all 31 new RDS seals, exact paired targets and
control changes, complete native draws, bridge medians/SDs, all model/family
weights and BMA policy moments. No material issue was found. The recorded
zero divergences and improved worst mixing come with greater computation,
rather than establishing a general backend or acceptance recommendation.

## 2026-10-22 follow-up: CmdStan PNB acceptance .999

The requested **CmdStanR .999** rerun completed separately in
`.scratch/single-market-intercept/pnb-m1t1-cmdstan-delta999-20261022/`.
All seven call lists match the original corrected CmdStan .99 calls after
changing only `adapt_delta` and the separate CSV output directory. The
diagonal `diag_e` metric, four chains, 1000 warmup/1000 retained draws per
chain, initialization zero, seed 19091, maximum depth 13, all scientific
arguments, full native k=1:5 support, model priors and policy definitions
remain unchanged. Direct checks against the saved Nutpie .999 fits confirm
identical full Stan data, model specification, source/data hashes and resolved
priors. This provides a matched .999 backend comparison under specification 12.

All seven fits contain every requested native draw; **21/21 normal bridge
repetitions are finite**, and all models retain **200/200 valid merger states**.
No diagnostics threshold suppresses evidence, and no positive-prior model is
dropped. The [model table](single_market_intercept_pnb_cmdstan_delta999_models.csv)
retains all 28 model/configuration rows from CmdStan/Nutpie at .99/.999,
including controls, diagnostics, elapsed times and source/data hashes. The
[bridge table](single_market_intercept_pnb_cmdstan_delta999_bridge_repetitions.csv)
retains all 84 repetition estimates across the four configurations.

| Diagnostic | CmdStan .99 | Nutpie .99 | CmdStan .999 | Nutpie .999 |
| --- | ---: | ---: | ---: | ---: |
| Divergences, seven models combined | 9 | 3 | 0 | 0 |
| Depth hits | 0 | 0 | 0 | 0 |
| Maximum R-hat | 1.004457 | 1.016147 | 1.006629 | 1.007238 |
| Minimum bulk ESS | 1159.10 | 469.09 | 976.57 | 679.72 |
| Minimum tail ESS | 884.42 | 798.82 | 839.88 | 805.26 |
| Summed warmup/sampling elapsed, seconds | 111.68 | 73.23 | 150.18 | 111.88 |

CmdStan elapsed in the last row sums each fit's longest native chain
warmup-plus-sampling time; Nutpie reports its native warmup/sampling total,
including callbacks/extraction. These are parallel-chain wall times summed
over models, not CPU totals, and backend timing scopes differ. CmdStan .999
process elapsed totals 153.11 seconds; its whole run, including R preflight,
evidence and policies, takes 331.45 seconds, compared with 258.13 seconds for
the single-process Nutpie .999 run. Existing compiled models/caches were reused;
CmdStan fits do not retain a separate compile-time measurement.

Raising CmdStan acceptance removes the observed divergences in this run but
does not improve every mixing metric. Its worst R-hat/bulk/tail summaries are
slightly weaker than in its .99 run. Under matched .999 acceptance, Nutpie
takes less measured sampling time while CmdStan retains more effective draws
for some coordinates, particularly outside share. The 42 prespecified scalar
structural/noise comparisons record direct bulk/tail ESS per native elapsed
second: Nutpie has higher bulk ESS/sec in 30/42 and higher tail ESS/sec in 28/42,
with median Nutpie/CmdStan ratios 1.296 and 1.182. Bulk ESS/sec ratios range
from .426 to 2.330. These are descriptive single-run results with the stated
timing scopes, not universal speed or efficiency guarantees.

| Model | Prior | CmdStan .999 posterior | Nutpie .999 posterior |
| --- | ---: | ---: | ---: |
| MC | 33.333% | 41.030% | 40.473% |
| B | 11.111% | 12.500% | 12.654% |
| C | 11.111% | 6.245% | 6.307% |
| BF | 11.111% | 12.800% | 12.812% |
| CF | 11.111% | 6.603% | 6.667% |
| BS | 11.111% | 12.431% | 12.570% |
| CS | 11.111% | 8.390% | 8.518% |

Family posteriors at .999 are **41.030% MC, 37.732% Bertrand, 21.238% Cournot**
under CmdStan and 40.473%, 38.035%, 21.492% under Nutpie, from one-third priors.
[Family probabilities](single_market_intercept_pnb_cmdstan_delta999_families.csv)
retain every configuration. The [matched evidence comparison](single_market_intercept_pnb_cmdstan_delta999_comparison.csv)
shows a maximum model-probability difference of .5576 percentage points and
maximum absolute log-evidence difference of .02298. All seven log-evidence
differences are within .729 of the conservative combined bridge uncertainty
defined above. [Parameter comparisons](single_market_intercept_pnb_cmdstan_delta999_parameters.csv)
retain both means, SDs, medians, intervals, mean MCSEs and direct R-hat/ESS
values, including the ESS/sec measures. The largest of 42 paired mean
differences is 2.2135 combined mean MCSEs.

The .999 BMA mean bank deposit-rate change is **-1.29404% CmdStan versus
-1.28198% Nutpie**, a difference of -.01206 percentage points. Their pooled
bank/draw SDs are 9.48925 and 9.17446 percentage points. The
[policy table](single_market_intercept_pnb_cmdstan_delta999_policy.csv)
preserves all four configurations, the existing equal-bank-weight statistic
and all successful-state counts. These pooled SDs are not aggregate posterior
SDs; no aggregate interval or policy MCSE is inferred from them.

All raw CmdStan CSVs and materialized posterior/diagnostic arrays are retained
with sealed fits, manifests, bridge attempts and execution receipts. The
original CmdStan .99 and both Nutpie runs remain intact. No mergerBayes
production source, data, selected result or scientific setting was changed;
no package production code or default sampler control was changed.

Independent review verified all 38 new RDS seals, complete native chains and
controls, all 21 bridge medians/SDs, model/family probabilities, BMA policy
moments and the six compact comparison tables. No material issue was found;
no additional fitting, bridge evaluation or policy calculation was performed
for that review.
