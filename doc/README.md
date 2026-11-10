# Development and evidence index

For public package use, start with the [supplied-state PNB guide](../vignettes/pnb-merger.Rmd),
then the optional [posterior workflow](../vignettes/posterior-workflow.Rmd) and
[model reference](../vignettes/supply-and-demand.Rmd). These guides do not run
inference or bridge sampling during documentation builds.

The 2026-11-08 total/fraction extraction closeout preserves the economic quality-SD alias
through nested conditional-quality recovery. Focused `posterior-quality` and
`selective-extraction` fixtures verify distinct economic/total scales and
nonzero-quality no-op policies without changing targets or specification 14.
The downstream installation includes additional uncommitted source changes;
its recorded Git SHA alone does not reproduce that build. See the mergerBayes
execution ledger for this closeout's source/build provenance and scope.

Opt-in [anchored Logit price scaling](validation/anchored_price_standardization.md)
has deterministic source checks; installed integration and empirical refits
remain pending the shared-library dependency refresh.

Current package version is `0.3.0.9002`, model specification 15. New Gaussian
relative-share Logit/CES fits default to persistent quality plus period demand
and independent measurement error. Counterfactuals retain the jointly recovered
quality, period demand and role states. See the
[demand-state contract](share_residual_demand_state.md).
Specification-14 fits retain their original model and verified frozen kernels;
new-model inference/evidence is required before adopting updated model weights.
Estimated-fit counterfactuals default to fitted posterior demand and role states;
observed-share anchoring is an explicit alternative. Specified fits retain
their supplied-state/calibration route. See the
[counterfactual default contract and local verification](validation/posterior_policy_default_20261108.md).
The historical simulation-default change does not invalidate compatible inference. Relative-share
Option B for Logit/CES was validated 2026-10-23 at version `0.2.0.9009`.
See the bounded
[implementation results](validation/relative_share_implementation/README.md),
including the preserved incomplete noiseless CES panel stress attempt.
The 2026-11-01 [coordinated release gate](validation/coordinated_release_gate_results.md)
passed 20/20 frozen inference criteria and 32/32 independent baseline/merger
FOCs on fresh current relative-share Logit/CES panel fits; broad recovery
remains separate. Package checks now have no errors, warnings, or NOTEs.
The previous committed baseline is `01e5a22` (version `0.2.0.9008`, specification 13).
This is the current operational agenda and navigation entrypoint. Topic contracts
and results remain the detailed sources; historical checklists do not select work.
Update this index when a completed task changes a status or verification route.

## Current development

The specification-15 package implementation includes the estimated measurement
prior, four production kernels, version-4 complete-demand policy plans and an
explicit original-model route. Six bounded known-truth fits completed with zero
divergences/depth hits, maximum R-hat 1.00934 and minimum bulk ESS 573; all
2,200 draw/market-year policies passed independent FOCs and no-op checks.
The ordinary non-HMC suite passed, with optional sampling skips and retained
warnings. See the [integration evidence](validation/residual_model_integration_results.md)
for native likelihood checks, sampler diagnostics, initial reporting failures,
saved-fit compatibility and package-check results. No full paper rerun or new
model probabilities are implied by this validation.

The completed [PNB residual and variance experiments](validation/share_residual_experiment_record_20261109.md)
are recorded for comparison; none selects an unfinished runner. Inverse
institution-assets variance substantially improved snapshot baseline shares;
hierarchical year effects improved historical temporal fit. These corrections
already enter the pilots' conditional quality recovery and both equilibria.
Production adoption remains separate, including revision of the recorded
hierarchy prior's heavy tail.

The [bounded measurement-prior diagnostic](validation/pnb_measurement_prior_diagnostic_20261109.md)
completed both new Bertrand snapshot/history fits with a prior on the previously
fixed measurement SD. Both had zero divergences/depth hits and adequate bounded
mixing; large-bank mean errors were at most .30/.59 percentage points, with
little margin change. All 200 new policies and all 800 comparison states were
valid. Snapshot precision remains prior-like; history shifts toward .079 but
the variance split remains prior-dependent. No production defaults, model
probabilities or manuscript estimates were changed. See the report for source
seals, independent checks, the reporting-only correction and replay commands.

The [public installation and supplied-state quick start](../README.md#installation)
separates package installation from optional estimation toolchains. Saved-fit
evidence accepts only proven LF/CRLF source portability; raw seals, substantive
source identity, data, priors and model-specification checks remain enforced.

Start with the requested behavior, current source/help, and its focused tests.
Use [package usage](../README.md) for the public interface and
[default-model history](default_model_history.md) for approved defaults and
frozen-fit provenance. Current defaults are Gaussian shares/margins, collapsed
market quality, and active rho where identified. Gamma hierarchy-SD defaults
do not retroactively change saved fits or establish universal HMC improvement.

| Task | Start here | Smallest relevant test filter |
| --- | --- | --- |
| Relative-share likelihood, reconstruction and evidence | [mathematical contract](relative_share_likelihood.md), [executed results](validation/relative_share_implementation/README.md), [native target runner](validation/relative_share_target_checks.R), [paired sampling runner](validation/relative_share_backend_validation.R) | `relative-share` |
| Absolute-share Logit and CES utility normalization | [normalization contract](single_market_utility_normalization.md), [correction results](validation/single_market_intercept_results.md), [fresh native checks](validation/single_market_intercept_validation_results.md) | `single-market-intercept|market-intercept|s0-pooling|conditional-k` |
| Estimation, controls, priors | [estimation](../R/estimate.R), [priors](../R/priors.R) | `stan-control-policy`, `hierarchical-sd-priors`, `market-intercept-sd-prior`, or `residual-sd-prior` |
| Sampler backends and paired validation | [backend adapters](../R/sampler_backend.R), [closeout summary](validation/nutpie_closeout_summary.md), [closeout protocol](validation/nutpie_closeout_protocol.md), [original paired results](validation/nutpie_backend_results.md), [evidence contract](validation/nutpie_evidence_design.md) | `sampler-backend|nutpie|backend-policy-equivalence|evidence-availability` |
| Quality/intercept reconstruction | [quality](../R/posterior_quality.R), [collapse contract](market_intercept_collapse_plan.md) | `posterior-quality|market-intercept-collapse` |
| Evidence and exact-k comparison | [comparison](../R/bma.R), [model universe](../R/model_universe.R), [conditional k](../R/conditional_k.R) | `bma|model-universe|restriction-comparison|conditional-k` |
| Policy states and counterfactuals | [draw plans](../R/posterior_draw_contract.R), [usage](statistical_toolkit_usage.md), [tariff contract](tariff_posterior_reuse_plan.md) | `posterior-policy-models|latent-quality-simulation|simulate-dispatch` |
| Conduct and cross-package models | [conduct](../R/conduct.R), [Stackelberg map](stackelberg_implementation_map.md) | `conduct|stackelberg-r-contract` |

`antitrust` owns common structural/counterfactual machinery; `coordination` owns
fringe/Stackelberg models. Read [DESCRIPTION](../DESCRIPTION) for dependency
minima and [CI](../.github/workflows/qa.yml) for dependency selection and the
canonical R version. Current minima are `antitrust >= 0.99.43.9000`,
`coordination >= 0.1.0`, and `cmdstanr >= 0.9.0`; optional policy integration
uses `trade >= 0.9.0`. Ordinary development CI resolves antitrust/trade `master`
and coordination `main`, records their SHAs, and pins CmdStanR. Release runs
use the full dependency SHAs described below. Historical executed reports
retain their earlier compatible stacks. The ordinary tier explicitly excludes Nutpie/BridgeStan;
the opt-in CmdStan tier installs its frozen toolchain before testing. Hosted
ordinary CI passed at baseline `55c6fdc` after its parser prerequisite guard
fix ([executed run](https://github.com/luciu5/antitrustBayes/actions/runs/39570316189)).
The current correction has separate local source/package/native checks; its
hosted run follows publication. A clean native Nutpie CI route remains future
work. See [CI readiness](validation/nutpie_ci_readiness.md).
Check installed versions and library paths before interpreting integration skips
or failures. The inspected default user library lacked the compatible stack;
choose an existing compatible library explicitly rather than reinstalling blindly.
An installed version alone does not verify its source revision.
Automatic `release/**` push and pull-request CI reads exact structural-package
SHAs from [`ci/release-dependencies.json`](../ci/release-dependencies.json);
manual release certification still requires and checks all three full-SHA
workflow inputs. Other development branches retain their existing bootstrap
refs. Release-pin CI also executes both embedded dependency metadata verifiers
against their selected source checkouts and installed packages; see the
[regression route](../ci/README.md#structural-dependency-metadata-regression).

## Verification commands

Run from the repository root with a compatible selected R library. This preflight
prints package provenance without fitting or installing anything:

```sh
Rscript --vanilla -e 'for (p in c("antitrust", "coordination", "cmdstanr", "testthat", "pkgload")) { cat(p, ": "); if (requireNamespace(p, quietly=TRUE)) cat(as.character(packageVersion(p)), find.package(p), "\n") else cat("missing\n") }'
```

Focused checks load this development source; replace the filter with the task's
entry above. The ordinary suite uses the same command without `filter`:

```sh
ANTITRUSTBAYES_RUN_STAN_TESTS=false Rscript --vanilla -e 'testthat::test_local(filter="posterior-quality|market-intercept-collapse", reporter="summary")'
ANTITRUSTBAYES_RUN_STAN_TESTS=false Rscript --vanilla -e 'testthat::test_local(reporter="summary")'
```

Record integration skips and the resolved dependency versions; the current
minimum is `coordination >= 0.1.0`. Routine checks can parse Stan
when CmdStan exists. Sampling is separately opt-in; short chains are smoke tests.
See the [test audit](test_suite_audit.md) for tier definitions and dated evidence.
Extended deterministic enterprise-fixture checks must also use current source:

```sh
Rscript --vanilla -e 'pkgload::load_all("."); testthat::test_file("doc/market_intercept_enterprise/fixture_checks.R")'
```

For statistical recovery, use the explicit prerequisites, output configuration,
and frozen specification-14 gates in [current CmdStanR validation](validation/README.md).
The [coordinated release runner](validation/coordinated_release_gate.R) uses
exactly two current relative-share Logit/CES panel fits; its
[dated report](validation/coordinated_release_gate_results.md) preserves every
failed criterion and fit provenance. Older absolute-share recovery results do
not qualify specification 14. The archived
124-task [enterprise runner](market_intercept_enterprise/README.md) needs a
reviewed backend port; its passing fixture checks do not validate current HMC.
For documentation-only work, resolve changed local links and run
`git diff --check`; do not rerun sampling or package installation.

For sampling-free release follow-up, use the
[sealed-fit/current-policy replay](validation/saved_gate_policy_replay.md).
It audits the historical inference seals and checks current structural-policy
FOCs in a new output directory; it does not relax the original runner's strict
source guard or relabel old fits as fresh inference. Its provenance regressions
run with `Rscript --vanilla doc/validation/saved_gate_replay_checks.R`.

## Consolidated agenda and topic evidence

Select follow-up work from the user's task. These statuses are evidence bounds,
not authorization to launch unfinished experiments. They replace the operational
checklists in the historical umbrella and recovery notes below.

| Topic | Recorded status and next decision | Detailed source |
| --- | --- | --- |
| Relative-share demand (Option B) | Implemented as the Gaussian default for Logit/CES, including rho, sampled/collapsed persistent quality, derived market-year intercepts, policy reconstruction and normalized evidence. Bounded validation: 5,903 native checks; 1,092 actual M4T4 structural/policy checks; 26 complete paired-backend fits and 78 finite bridge repetitions across the original complete cases and a separately labeled finite-error CES panel. Original noiseless CES M4T4 stress remains incomplete. Full opt-in suite: 7,619 passing expectations, no failures/errors/skips, 30 retained warnings. Package check: no errors/warnings, one existing NOTE. Broad recovery and downstream paper adoption remain separate. | [contract](relative_share_likelihood.md), [executed results and limitations](validation/relative_share_implementation/README.md), [production validation runner](validation/relative_share_backend_validation.R), [earlier prototype](validation/option_b_relative_shares/README.md) |
| Market/outside-share rotation | Implemented experimental coordinates with bounded sampler evidence; scaling and broader equivalence remain incomplete. | [rotation results](market_s0_rotation_results.md) |
| Posterior-state policies | Default for estimated fits; explicit observed-share calibration remains available. Fixed plans align parameters, recovered quality and roles across equilibria and scenarios. The frozen historical observed/latent benchmark and broader statistical coverage remain distinct. | [current default contract](validation/posterior_policy_default_20261108.md), [historical policy results](latent_quality_policy_results.md), [usage](statistical_toolkit_usage.md) |
| Hybrid margins | Frozen seven-fixture benchmark completed; outlier/predictive evidence remains pending. | [hybrid results](collapsed_gaussian_hybrid_results.md) |
| Model universe | Public descriptors and evidence aggregation implemented; restricted-family policy BMA is outside that interface. | [usage](statistical_toolkit_usage.md), [design contract](model_universe_specification.md) |
| Hierarchy priors and recovery | Approved defaults and bounded comparisons recorded; existing fits retain frozen priors. | [default history](default_model_history.md), [paired results](validation/hierarchical_prior_recovery_results.md) |
| Market-intercept collapse | Implemented opt-in; fixed-target validation differs from broad sampler certification. Enterprise port remains required. | [validation](market_intercept_collapse_validation.md), [enterprise status](market_intercept_enterprise/README.md) |
| Tariff posterior reuse | Implemented initial route subject to provenance/equivalence conditions; heterogeneous post tariffs remain a separate extension. | [results](tariff_posterior_reuse_results.md), [contract](tariff_posterior_reuse_plan.md) |
| Sampler presets and diagnostic verdicts | Proposed, unimplemented; current controls remain explicit. | [proposal](sampler_usability_followup_plan.md) |
| Prior-predictive diagnostics | Historical requested follow-up; no general public prior-predictive interface was found in this review. Scope separately before implementation. | [historical roadmap](statistical_toolkit_followup_plan.md), [inference contract](market_prediction_inference_design.md) |
| Recovery and model discrimination | Current bounded CmdStanR gates exist; broad repeated-sample coverage and general conduct identification remain open. | [validation entrypoint](validation/README.md) |
| Single-market Logit normalization | Corrected in both Logit kernels: one direct `Normal(0,1)` inside intercept, no single-market hierarchy SD; multi-market hierarchy, outside-share pooling unchanged; the CES kernels followed at specification 13 (next row). Deterministic 10%/20%/40% support and B/C/CF/CS posterior-policy identities pass. Full ordinary suite: 6,454 assertions; CmdStan: 590; package check: zero errors/warnings, one existing NOTE; hosted ordinary CI passed at `112c448`. Fresh validation: 25/25 tasks, 138 paired target states, 18/18 finite bridge repetitions; evidence differences agree within MC uncertainty, with all three precision caps inconclusive and actual mixing warnings retained. Local revised-paper B M1T1 follow-up: 56/56 fits, 168/168 finite bridge repetitions, 60.067%/59.527% Bertrand-family posterior from 25% prior; merger/tariff results and small explicit cost attrition retained. PNB M1T1: 7/7 fits, 21/21 finite bridge repetitions, complete 200-state policies per model; family probabilities 40.621% MC/37.946% Bertrand/21.433% Cournot from equal thirds. Matched Nutpie PNB rerun: 7/7 fits, 21/21 finite bridges, exact source/data/prior identity, maximum posterior model-probability difference .576 percentage points, 3 divergences versus 9 CmdStan; Nutpie C has weaker mixing (maximum R-hat 1.016, minimum bulk ESS 469). All merger states retained. Nutpie PNB acceptance .999 sensitivity: 7/7 fits, 21/21 finite bridges, zero divergences/depth hits, maximum R-hat 1.00724 and minimum bulk ESS 680; native time rises 52.8%, with stable probabilities/policies. Matched CmdStan .999 rerun: 7/7 fits, 21/21 finite bridges, zero divergences/depth hits, maximum R-hat 1.00663 and minimum bulk ESS 977; matched .999 backend model probabilities differ by at most .558 percentage points and BMA merger means by .0121 percentage points. These bounded runs do not establish repeated-sample recovery or universal sampler superiority. Existing fits require refitting under specification 13 (the global guard also rejects specification-12 Logit fits); downstream mergerBayes received the specification-12 integration and .999 banking acceptance floor (Supreme cross-section MC retains .9999) at `d7520d2`, with targeted cache, banking, reporting and scientific-freeze checks. Its newer runner at `dac8a95` accepts specification 12 or later and records package builds as provenance; this does not bypass the package's specification-13 saved-fit guard. | [contract](single_market_utility_normalization.md), [results and local follow-ups](validation/single_market_intercept_results.md), [fresh native checks](validation/single_market_intercept_validation_results.md) |
| Single-market CES normalization | Applied to both CES kernels at specification 13: one direct `Normal(0,1)` inside intercept, no single-market hierarchy SD; multi-market hierarchy, CES revenue/quantity outside-share and outside-price handling and supply FOCs unchanged. The global specification bump also invalidates specification-12 Logit fits (kernels unchanged); downstream runners pinning 12 must update and refit. Verification is the deterministic and fixed-parameter Stan tests in the linked test file; no sampling-based recovery or bridge-evidence validation was run for CES. | [contract](single_market_utility_normalization.md), [tests](../tests/testthat/test-ces-intercept-normalization.R) |
| Conduct-schedule feasibility | `conduct_schedule()` now requires the aggregate condition `sum(q) <= 1 - s0` (tolerance 1e-8; `q` may list a subset of firms), and `conduct_distance()` bounds the profiled outside share by `1 - sum(q)` rather than `1 - max(q)`. | [theory](../R/theory.R), [tests](../tests/testthat/test-theory.R) |
| Optional Nutpie sampler | Historical specification-11 validation: public normalized evidence, mixed-backend BMA and exact-source cache reconstruction supported; 42/42 closeout tasks completed, 90/90 bridge repetitions finite, all four empty-cache rebuilds validated. Full conditional-k evidence/weights and policy means agree. Precision remains inconclusive for 1/21 evidence and 13/28 conservative policy comparisons; attrition/diagnostics stay visible. Original 34 fits and corrected CES replay: 1,385 scored posterior/policy checks within MCSE tolerances. CmdStanR remains default; Nutpie diagonal is a validated alternative for the tested targets, low-rank an expert option. Ordinary suite: 5,166 assertions with native dependencies, 5,159 without; CmdStan: 590 passed. Package check: zero errors/warnings, one existing NOTE. Hosted ordinary CI passed after the parser prerequisite fix. Original timings and posterior checks describe the earlier target; see the normalization correction row for fresh specification-12 checks. | [closeout summary](validation/nutpie_closeout_summary.md), [uncertainty audit](validation/nutpie_closeout_uncertainty_results.md), [regressions](validation/nutpie_closeout_regressions.md), [original results](validation/nutpie_backend_results.md), [promotion pilot](validation/nutpie_certification_pilot_results.md), [saved-fit results](validation/nutpie_serialization_results.md), [CI readiness](validation/nutpie_ci_readiness.md) |

## Historical checkpoints

These retain original plans, failures, source identities, and execution records.
Session IDs, temporary libraries, worktree constraints, and statements about
HEAD describe their recorded dates. Verify present source and durable artifacts
before reusing anything; this review does not establish later experiment outcomes.

- [Original toolkit roadmap](statistical_toolkit_followup_plan.md): overlapping
  topic work is consolidated in the table above; unique scientific caveats remain.
- Performance audit: [frozen plan](performance_audit/plan.md),
  [dated results](performance_audit/results.md),
  [environment recovery](performance_audit/recovery_status.md), and
  [old resumption checklist](performance_audit/recovery/next_steps.md).
  Incomplete vectorization evidence establishes no HMC speedup; the optional
  bridge candidate remained held in the recorded audit.
- Stackelberg: [execution log](stackelberg/lead_execution_state.md) and
  [September evidence index](stackelberg/results_index.md). Historical pilot
  progress is not a current process status or broad recovery result. Use current
  source/tests for supported APIs and retain the reports' economic limitations.

Keep math, contracts, executed results, and open proposals distinct. Add links
from this index when a topic becomes relevant; do not recursively load or rewrite
the documentation tree to answer a bounded task.
