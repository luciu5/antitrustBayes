# Market-quality sampler and application evidence

Status: implementation, native validation and independent package tests passed.
The first six frozen Bertrand fits have completed and all miss the prespecified
usability gate. The remaining thirty primary cases have not run. All eight
reconstructed-input PNB fits have completed; seven pass the full gate. Bertrand
has two divergences and therefore misses that gate. Historical
collapsed-Gaussian sampler evidence is recorded separately and is not a pass
for these new multi-market, rho-active targets.

## Frozen comparison

The comparison distinguishes three specifications:

| label | shares / margins | quality parameterization | scope | rho |
|---|---|---|---|---|
| legacy | Student-t / Student-t, df=5 | sampled | global | active |
| Gaussian expanded market | Gaussian / Gaussian | sampled | market | active |
| proposed routine model | Gaussian / Gaussian | collapsed | market | active |

The two Gaussian market models have exactly the same statistical prior and
likelihood after quality integration. The legacy model differs in both
residual family and multi-market quality prior; its comparison cannot isolate
Student-t robustness from the scope change. Single-market Gaussian scope
comparisons isolate the normalization/regression implementation exactly.

Primary fits use frozen data, common prior settings, seeds and sampler budget.
No fit-specific tuning follows diagnostic failures. Any existing coordinate
rotation sensitivity remains separately labeled and preserves the statistical
model. All failed and interrupted fits remain in the benchmark ledger.

## Scaling and PNB

Scaling includes target/gradient timings for M=10,100,500,2000, active rho and
an exact family; HMC runs are recorded where feasible. No dense global
covariance or production Cartesian product of states is permitted.

PNB is one market. Global/market Gaussian scope is therefore mathematically
identical; Gaussian versus Student-t results need not agree. Canonical B, MC,
C, A fits precede richer-family fits and normalized bridge evidence. The
non-overlapping model universe represents shared endpoints once. Model prior
weights must retain their documented scientific provenance; they are never
chosen after seeing evidence.

Complete-universe policy comparisons use a common supported demand anchoring.
Existing auction/bargaining latent-quality anchoring restrictions are explicit;
latent-quality comparisons for supported conduct families are reported
separately rather than combined under inconsistent anchoring.

## Current recommendation

Do not change routine defaults yet. Independent integration/factorization,
native target/gradient and package gates have passed. Multi-seed usable HMC
and PNB evidence remain mandatory promotion gates.

## First active-rho paired result (one seed; provisional)

The multi-market Bertrand fixture at seed 271828 fails the frozen usability
gate in both Gaussian parameterizations. Collapsed market quality has 17
post-warmup divergences, minimum BFMI .761, maximum R-hat 1.0118 and minimum
bulk ESS 599.5; sampled market quality has 20 divergences, BFMI .501, maximum
R-hat 1.0150 and minimum bulk ESS 393.5. Wall times are 125.25 and 67.15 seconds,
giving minimum bulk ESS/wall-second 4.79 and 5.86 respectively. This first result
does not establish a collapse speed advantage or acceptable geometry with
active rho. Original fits/warnings are retained; the second seed is reported
below. Other primary fixtures remain frozen and unrun. No primary settings
were tuned after this result.

## Located PNB comparator

The user-updated mergerBayes checkout at commit `3dc20ce` contains the recent
`rich_main` run; its recorded antitrustBayes source is `ea29e96`. It used
Student-t shares/margins, sampled global quality, active rho, fixed df=5, the
baseline parameter prior, and legacy market/outside-share coordinates. All
eight old fits have zero divergences and treedepth hits, but every model fails
minimum BFMI >0.30 (range approximately 0.041–0.107). Auction additionally has
maximum R-hat approximately 1.036. These failures remain part of the comparator.

The frozen PNB comparison will retain its four chains, 8,000 iterations, 2,000
warmup, adapt_delta=0.99, max_treedepth=12, and seed 1960. Exact families use
k=1:5 with uniform within-family probabilities, matching that run. Its
non-overlapping universe has prior weights MC=A=1/4 and
B=C=BF=CF=BS=CS=1/12. Equal-1/8 sensitivity, if reported, reuses the same
evidence and remains a separately labeled prior sensitivity. Single-market
market/global Gaussian equality is a separate gate from Gaussian/Student-t
posterior similarity. New Gaussian PNB fits are reported separately below.

The original rich-main input file is unavailable in this checkout. The user
approved the bundled 39-bank PNB data for a separately labeled diagnostic.
The bundled mean margin is 0.70985 versus 0.47369 in rich-main, despite matching
price summary statistics. Its data/settings are frozen before fitting. Any
comparison with rich-main must retain this input confounding; it is not a
likelihood-only robustness comparison.

## Both initial Bertrand seeds

The complete initial comparison is retained in
[the benchmark ledger](market_quality_hmc_benchmark.csv). All fits use the
same frozen sampler budget; no settings were changed after these results.

| seed | specification | wall seconds | divergences | minimum BFMI | maximum R-hat | minimum bulk ESS / wall second |
|---|---|---:|---:|---:|---:|---:|
| 271828 | Gaussian collapsed market | 125.25 | 17 | .761 | 1.0118 | 4.79 |
| 271828 | Gaussian sampled market | 67.15 | 20 | .501 | 1.0150 | 5.86 |
| 271828 | Student-t sampled global | 66.19 | 1 | .641 | 1.0334 | 3.58 |
| 314159 | Gaussian collapsed market | 137.25 | 10 | .695 | 1.0047 | 5.28 |
| 314159 | Gaussian sampled market | 81.67 | 13 | .423 | 1.0063 | 5.90 |
| 314159 | Student-t sampled global | 47.13 | 0 | .665 | 1.0113 | 7.45 |

None reaches maximum treedepth. Collapsing improves minimum bulk ESS and BFMI
in both Gaussian comparisons but costs more wall time, so minimum bulk ESS per
second is lower. These minima cover sampled coordinates and reconstructed
economic outputs; they are not a direct policy ESS calculation. The second
Student-t fit fails the R-hat gate despite zero divergences. Statistical-model
differences prevent attributing the Student-t/global comparison solely to
residual robustness.

Descriptive first-seed analysis finds correlations around -.992 between
reconstructed logit outside share and the market intercept in both Gaussian
fits. Their sum has only .002--.065 correlation with the intercept. This
confirms that a separate ridge remains in these legacy sampling coordinates;
it does not identify the cause of every divergent transition or establish that
a rotation sensitivity would succeed.

## Executed native scaling

[Native timings](market_quality_scaling.csv) use active rho=.31, market-level
outside-share/alpha pooling, and exact BF states 0:3 independently in each
market. Each size runs in its own R process, with 25 target and gradient
evaluations after fixed-state construction. Targets, gradients and generated
component probabilities are finite and probabilities sum to one by market.

| markets | observations | target ms | gradient ms | peak RSS MiB |
|---:|---:|---:|---:|---:|
| 10 | 80 | .32 | .64 | 519.5 |
| 100 | 800 | 2.88 | 5.96 | 522.8 |
| 500 | 4,000 | 15.24 | 31.76 | 579.1 |
| 2,000 | 16,000 | 63.20 | 132.28 | 789.2 |

Runtime growth is approximately linear in this fixed-K design. RSS includes
the R/Stan/native-library baseline; it is not kernel-only memory. No dense
global covariance or joint K^M configurations were constructed. HMC time and
ESS/sec at these sizes remain unmeasured: the small primary sampler gate has
not passed, so this table cannot establish a practical large-market posterior
sampler.

## Recovered historical PNB input

Read-only Git recovery located mergerBayes commit `4be1c47`, blob
`a26a38e07c2a9127a77a4a43fcb9d158a13be6da`, containing the historical
`supreme_data.RData` (216,541 bytes). Applying the existing rich-main preparation
matches its recorded PNB summary, including mean markdown .4736931142. The
reconstructed data and unchanged application settings were frozen before any
PNB fit, and are now used for the primary application diagnostic. The original
process's serialized input remains unavailable, so this is an archived-input
reconstruction. The bundled alternative is preserved separately without fits.

The first PNB Bertrand fit completed in 53.001 wall seconds (171.452 summed
native seconds over four parallel chains), with two post-warmup divergences.
It therefore fails the frozen usability gate. No sampler or prior setting was
changed after this outcome. The global-quality regression produces identical
checked posterior draws and the same divergences, as single-market equivalence
requires.

The initial BF driver attempt
failed before HMC because an explicit legacy `use_cutoff=FALSE` contradicted
`structure="fringe"`. That attempt is retained. Removing the unnecessary flag
leaves the frozen data, priors, exact-k support and HMC settings unchanged.
Fixed-state preflight checks confirm that the inherited extra R data flag is
unused by the exact-k Stan program, including native density, gradient and
generated-quantity parity. The remaining frozen fits use the repaired driver.

## Complete PNB sampler summary

[The compact PNB ledger](market_quality_pnb_hmc.csv) retains all eight fits.
All use the archived-input reconstruction, Gaussian shares/margins, collapsed
market quality, active rho, and legacy market/outside-share coordinates.

| model | wall seconds | divergences | minimum BFMI | maximum R-hat | minimum bulk ESS | full gate |
|---|---:|---:|---:|---:|---:|---|
| B | 53.00 | 2 | .863 | 1.00058 | 7,474 | fail |
| MC | 56.56 | 0 | .917 | 1.00044 | 6,945 | pass |
| C | 60.26 | 0 | .940 | 1.00061 | 7,949 | pass |
| A | 67.50 | 0 | .917 | 1.00069 | 7,684 | pass |
| BF | 148.52 | 0 | .912 | 1.00053 | 6,792 | pass |
| BS | 135.25 | 0 | .916 | 1.00039 | 6,876 | pass |
| CF | 161.15 | 0 | .924 | 1.00086 | 7,823 | pass |
| CS | 164.96 | 0 | .946 | 1.00048 | 8,983 | pass |

There are no treedepth hits or constant estimated coordinates. Minimum tail
ESS exceeds 5,294 in every model. These are executed diagnostic results, not
claims based on the absence of warnings. Bridge calculations, policy draws
and complete-universe BMA are separate gates. Complete-universe BMA remains
withheld while any contributing fit is unusable.

The two B divergences are consecutive saved draws in one chain. At those
draws, sigma_share is below its first posterior percentile, tau above its
98th percentile, and rho unusually negative (-.808 and -.881). This rare-tail
association is a diagnostic clue, not evidence identifying the cause. The
[CSV diagnostic](market_quality_pnb_divergence_diagnostic.py) retains the
divergent states and all chain counts. Native CSV precision limits these
descriptive comparisons; they are not fixed-target or gradient checks.

[The parameter comparison](market_quality_pnb_parameters.csv) records both
new and historical posterior summaries with sampler eligibility. For example,
mean alpha is 1.190 versus 1.061 in B and 1.847 versus 1.861 in C; mean outside
share is .0191 versus .0141 in B and .0350 versus .0359 in C. Historical
Student-t summaries are diagnostic comparisons from BFMI-failing fits, not a
validated reference posterior. The archived-input provenance qualification
also remains in force.

[The exact-k posterior](market_quality_pnb_k.csv) is available for every richer
family. BS probabilities remain approximately .20 for each k: these data do
not substantially update its uniform prior. CF puts .483 on k=1; CS puts
.379 on k=5. BF is close to uniform (.191--.208). This is posterior evidence
conditional on each family and its frozen specification, not proof of a
particular strategic-core size or a scaling defect. Model probabilities and
policy conclusions require the separate bridge/simulation calculations.

## Active-rho rotation sensitivity: pre-HMC gate

The separately labeled sensitivity copies the two original collapsed
Bertrand fixtures exactly, changes only market/outside-share sampling
coordinates, and maps economic initial states. The priors, data, seeds,
sampler budget and native model library remain fixed. Original failed primary
fits remain in the ledger. A sampled-quality rotated pair is prespecified,
conditional on both collapsed rotated fits passing the full gate.

Eight initial-state checks and two nonzero economic probes pass before HMC.
The probes include rho=.31, nonzero market intercepts and unequal outside
shares. Maximum normalized target difference is 2.85e-14; automatic-gradient
chain-rule and finite-difference relative errors are below 7.51e-10. Recovered
economic parameters and checked likelihood/prediction outputs differ by at
most 8.89e-16. This establishes coordinate equivalence for the tested states;
sampler improvement is a separate executed comparison, pending below.

## Active-rho rotation: executed short-budget results

[The rotation ledger](market_quality_rotation_hmc.csv) includes the original
no-sample launch, its paired not-run case, and the two valid retries. The
failed launch used a missing seed column, producing integer(0); the operational
retry restores exactly the original primary seeds. No posterior was sampled
in that first launch. The immutable failed artifacts remain separately saved.

| seed | legacy / rotated wall seconds | legacy / rotated divergences | rotated max R-hat | rotated minimum bulk ESS/sec | full gate |
|---|---:|---:|---:|---:|---|
| 271828 | 125.25 / 26.09 | 17 / 0 | 1.01045 | 26.27 | fail: R-hat |
| 314159 | 137.25 / 27.63 | 10 / 0 | 1.00784 | 28.05 | pass |

The sampling-coordinate c/intercept correlations are at most .0764/.0677;
reconstructed outside-logit/intercept correlations remain about -.993, as
expected for the same economic posterior. Total leapfrog counts fall from
1,087,392/1,082,660 to 176,366/187,190. An independent native-CSV calculation
reproduces divergences, BFMI, trajectory counts, R-hat and ESS to CSV precision.
The normalized target and gradient checks establish the same statistical model.
The first rotated run's maximum R-hat of 1.0104476, on the share residual scale,
is a marginal Monte Carlo precision miss of the 1.01 threshold. It is not
evidence that the economic model or coordinate transformation failed. The
strict gate label is retained separately from that interpretation. Fully usable
paired posterior equivalence and default promotion remain unestablished.

The original conditional sampled-quality short pair is not run because the
two collapsed gates did not both pass. One separately labeled budget-length
sensitivity is specified before new posterior sampling: four chains, 2,000
iterations, 1,000 warmup, the same seeds/fixtures/priors/economic starts and
unchanged .95/depth-13/diagonal-metric controls. Its sampled-quality pair is
conditional on both longer collapsed runs passing. It does not replace any
original failure, and speed comparisons between quality modes must use equal
budgets. Longer chains can improve Monte Carlo precision; they do not by
themselves prove that divergent exploration has been repaired.

## Executed doubled-budget sensitivity

Both collapsed runs pass, allowing the frozen sampled-quality pair to run.
These four fits use exactly the same Gaussian market-quality statistical
model, active rho, rotated coordinates and 2,000/1,000 iteration/warmup budget.
No priors or sampler controls changed.

| seed | quality coordinates | wall seconds | divergences | BFMI | max R-hat | min bulk ESS/sec | min tail ESS/sec |
|---|---|---:|---:|---:|---:|---:|---:|
| 271828 | collapsed | 38.29 | 0 | .777 | 1.00219 | 40.45 | 47.53 |
| 271828 | sampled | 27.33 | 2 | .549 | 1.00459 | 32.29 | 56.28 |
| 314159 | collapsed | 37.99 | 0 | .721 | 1.00272 | 48.26 | 38.58 |
| 314159 | sampled | 27.72 | 0 | .458 | 1.00385 | 31.73 | 53.92 |

There are no treedepth hits, constant estimated coordinates, or low-ESS
coordinates. Collapsing improves minimum bulk ESS/sec by 1.25 and 1.52 times
and BFMI at both seeds, despite costing more raw wall time. Minimum tail
ESS/sec is lower, and one sampled run has two divergences. These are diagnostic
minima, not measured policy ESS/sec; there is no uniformly faster winner.

Independent CSV calculations reproduce the trajectory and posterior
diagnostics to saved precision. The
[economic posterior comparison](market_quality_rotation_posteriors.csv)
records means, intervals and estimated
Monte Carlo uncertainty for both modes. Alpha and outside-share mean
differences are below 1.7 combined MCSE. At seed 314159 the margin scale
differs by about .00076, or 3.42 estimated combined MCSE; this discrepancy is
retained rather than described as established posterior equality. The
combined-MCSE calculation is descriptive and treats chains from the separate
fits as independent. A sampled fit with divergences cannot establish fully
usable paired posterior equivalence.

The initial marginal R-hat miss disappears at the longer budget. This provides
two-seed usable collapsed evidence for this Bertrand fixture, while the
unrun conduct fixtures, PNB Bertrand divergences and application evidence
remain separate promotion gates. It does not replace the original failures.
