> Specification 14 (`0.2.0.9009`) defaults to relative shares. The current
> focused route is [the native target checks](relative_share_target_checks.R)
> followed by [paired backend validation](relative_share_backend_validation.R),
> including Logit/CES M4T4. [Executed implementation results](relative_share_implementation/README.md)
> retain the incomplete noiseless CES panel stress attempt separately from the
> completed finite-error panel. The earlier recovery gate and numeric results below
> concern older absolute-share specifications; their identification and
> intercept assumptions do not certify the new likelihood.

# Specification-14 coordinated release gate

Use the extended-tier [bounded runner](coordinated_release_gate.R) with current
development source and a compatible explicitly selected R library:

```sh
R_LIBS=/path/to/compatible/library \
  Rscript --vanilla doc/validation/coordinated_release_gate.R --sample
R_LIBS=/path/to/final/compatible/library \
  Rscript --vanilla doc/validation/coordinated_release_gate.R --analyze --policy
```

The default durable ignored output is
`.scratch/coordinated-release-gate` in the authorized antitrustBayes workspace;
set `AB_RELEASE_GATE_OUTPUT` to select a fresh output root. An existing manifest
blocks resampling. Fits and failures are preserved before analysis. Exactly two
current relative-share fits use `relative_validation_fixtures()`'s Logit C M4T4
and the prescribed finite-noise CES C M4T4 transformation. Both have four chains,
1,000 total iterations, 500 warmup, `adapt_delta=.99`, depth 13, and recorded seeds.
No bridge sampling or historical absolute-share fixture enters this gate.

Acceptance is frozen before sampling: each market's structural truth lies in
its 95% interval with median relative error at most 8%; at least three of the
four t1 reference outside-share truths lie in 95% intervals with mean absolute
median error at most .06; active sampled parameters have maximum R-hat at most
1.05, minimum bulk ESS at least 100, and zero divergences. All 16 outside-share
cells are reported. The t1 reference is necessary because the current Logit
fixture leaves outside shares unpooled across years. Both error-scale medians
must be between .5 and 2 times their fixture's actual log-error RMS, including
the finite CES errors; no historical Logit noise band is reused for CES.

Raw demand, ownership and cost primitives independently check baseline and
merger FOCs at the true state and one posterior draw fixed before fitting.
The maximum allowed residual is `1e-6`. CES follows the legacy Cournot contract,
which holds its calibrated outside-revenue alpha fixed over policy changes.
Source, kernel, data, fitted specification, sampler settings and runtime library
provenance are sealed. [Dated results](coordinated_release_gate_results.md)
record execution and any failed criterion without changing these thresholds.

## Sampling-free follow-up

Use the [saved-fit/current-policy replay](saved_gate_policy_replay.md) when
packaging/dependency metadata changes while the sealed R/Stan inference source
is unchanged. Historical seals remain intact; current runtime and policy FOCs
are recorded separately. This does not certify new statistical recovery.

# Historical absolute-share CmdStanR validation

The following runner and reported numeric results predate specification 14.
They remain historical evidence and do not qualify the current default target.


Run the known-truth recovery gate after installing the current package and
CmdStan:

```sh
AB_RECOVERY_OUTPUT=/tmp/antitrustBayes-recovery \
  Rscript --vanilla doc/validation/known_truth_recovery.R
```

Set `AB_RECOVERY_CMDSTAN` to a **writable** CmdStan installation for the
default full gate, and `AB_RECOVERY_CACHE` to a writable model cache if the
default temporary cache is unsuitable. The script uses fixed seeds, two HMC
chains per fit, 250 warmup and 250 retained iterations per chain. It writes
recovery, evidence, and diagnostic CSVs to `AB_RECOVERY_OUTPUT` and exits
nonzero if any required criterion fails. Use
`AB_RECOVERY_SCOPE=recovery` for the shorter parameter-recovery check.

The fixture is generated without antitrustBayes helpers. It has six markets,
three years per market, and five single-product firms per market. Conditional
Logit shares come from a softmax with known structural price coefficient
`alpha = 2.5`. The outside share is constant across years in each market.
Output margins follow the single-product Bertrand FOC
`1 / Lerner = price * alpha * (1 - (1 - s0) * inside_share)`
with independently generated margin noise. The resulting 90 observations
are fitted with the current public `bayes_est()` CmdStanR path.

The gate requires the true alpha inside its 95% posterior interval and its
posterior median within 8% of truth; at least five of six true outside
shares inside their 95% intervals with mean absolute median error at most
six percentage points; a plausible recovered margin-noise scale; key
rank-normalized R-hat at most 1.05, minimum bulk ESS at least 100, and zero
divergences. It also checks the economic invariants that a common alpha
does not vary across markets and a market-level outside share does not vary
across years in any posterior draw. These thresholds were chosen for a
bounded known-truth diagnostic, not as a claim of repeated-sample coverage.

The full gate also fits Bertrand and Cournot to a second, independently
generated **one-market**, eight-year Bertrand market. Its outside share is
constant across years and a free market intercept is inactive, so Cournot
cannot imitate Bertrand by raising an unconstrained market outside share.
The gate uses the public `bayes_compare()` bridge path with two repetitions.
It requires both fits to have sound chain diagnostics, both bridge
repetitions to be finite and consistent, the current Stan source hash to
match the fit, a Bertrand/Cournot Bayes factor of at least 10, and posterior
model probability for the true Bertrand model of at least 0.9 under equal
model priors. The completed full gate gave log marginal likelihoods 135.627
(Bertrand) and 82.129 (Cournot), with bridge repetition SDs 0.023 and 0.003
log units.

A first multi-market evidence pilot did **not** meet that standard: Cournot
could inflate its market outside shares, and its bridge repetitions differed
by about five log-evidence units. That more general model-selection problem
remains open; the successful one-market check establishes a bounded current
CmdStanR evidence path, not universal conduct identification.

The historical 124-task enterprise benchmark in
`doc/market_intercept_enterprise/` still requires an RStan-to-CmdStanR port.
These bounded checks establish current-backend parameter recovery and
model-evidence preference in two identified designs; they do not establish
broad repeated-sample coverage or replace the enterprise study.

## Paired hierarchy-SD prior validation

`hierarchical_prior_recovery.R` compares the package's Gamma(2, rate=5)
alpha/outside-share SD defaults with their historical HalfNormal(0,.5)
priors. This is a separate bounded diagnostic; it is not a production
prerequisite or a change to the paper's explicit alpha-SD Gamma rate of 10.

```sh
AB_HIERARCHY_OUTPUT=/tmp/antitrustBayes-hierarchy-recovery \
AB_HIERARCHY_CMDSTAN=/path/to/writable/cmdstan \
  timeout --signal=TERM --kill-after=30s 7200 \
  Rscript --vanilla doc/validation/hierarchical_prior_recovery.R
```

Dependencies and the current package must already be installed in the R
library. `AB_HIERARCHY_SCOPE=prepare` validates the fixtures without HMC.
The default output is `antitrustBayes-hierarchy-recovery` in the system's
temporary root, outside R's session directory so results survive R's exit.
An output directory containing an inference manifest cannot be used to reset
the deadline. After a reporting interruption, `AB_HIERARCHY_RESUME=true`
reuses saved posteriors and identical fit arguments under the original
deadline. An attempted fit without a saved posterior is not silently refitted.
The full run has exactly four fits: small and
moderate heterogeneity under each prior, with four chains, 4,000 iterations,
2,000 warmup, four cores, dense mass matrix, adapt_delta .995 and depth 14.

Each fixture contains eight markets, four periods and five single-product
Bertrand firms. Normal innovations generate log-alpha SD .05/.30 and
logit-outside-share SD .08/.40; market-intercept SD is .25 in both fixtures.
All primitive and measurement innovations are shared across fixture scales,
and each prior pair receives identical observations and sampler seeds.
Independent algebra sets prices, persistent qualities and positive physical
costs satisfying Bertrand FOCs. A global intercept level enforces the fitted
price-centering convention, without changing the hierarchy SDs. All 64
market-period equilibria are freshly solved and checked before sampling.
Compositional share error (.025) and log-margin error (.035) are applied
afterward. The structural hierarchy is matched; the correlated compositional
measurement process remains an approximation to the independent Gaussian
log-share working likelihood. This is not simulation-based calibration.

Both tested hierarchies are active. Other priors, including the market-
intercept Half-Normal(.5) prior, explicitly fixed for this comparison, stay
fixed. Recovery tables distinguish the
known generating SD from its realized finite-market sample SD. Diagnostics
include maxima/minima across all reported quantities (covering all active
parameters), selected recovery-variable diagnostics, divergences, BFMI and
depth hits. Evidence reports three bridge repeats for the same true Bertrand
model under each prior; it is not an eight-conduct recovery comparison.

Policy checks use 200 shared posterior parameter/quality states for a fixed
market-1/period-1 HH merger with zero efficiencies and a statutory 25% tariff
on the other firms (internal consumer-price fraction .2). Both predictions
use the same paired admissible states and actual known-economy oracles.
Failures/warnings are retained. Outputs include `manifest.rds`, `fixtures.rds`,
`equilibrium_validation.csv`, `recovery.csv`, `diagnostics.csv`, `evidence.csv`
and `policy.csv`; completed fits are saved before downstream work so a timeout
does not erase them. Two fixed fixtures assess sensitivity, not general
coverage or a universally preferable prior.

The [completed local comparison](hierarchical_prior_recovery_results.md)
records the four fits, independent audit, recovery, policy and evidence
results, and execution warnings. Both families recover these fixtures;
Gamma is not a uniform sampling improvement. This validation adds no
production prerequisite.
