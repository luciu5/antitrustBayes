# Development and evidence index

Reviewed 2026-10-21 for the single-market Logit normalization correction,
starting at `55c6fdc`, package version `0.2.0.9007`, model specification 12.
This is the current operational agenda and navigation entrypoint. Topic contracts
and results remain the detailed sources; historical checklists do not select work.
Update this index when a completed task changes a status or verification route.

## Current development

Start with the requested behavior, current source/help, and its focused tests.
Use [package usage](../README.md) for the public interface and
[default-model history](default_model_history.md) for approved defaults and
frozen-fit provenance. Current defaults are Gaussian shares/margins, collapsed
market quality, and active rho where identified. Gamma hierarchy-SD defaults
do not retroactively change saved fits or establish universal HMC improvement.

| Task | Start here | Smallest relevant test filter |
| --- | --- | --- |
| Single-market Logit utility normalization | [normalization contract](single_market_utility_normalization.md), [correction results](validation/single_market_intercept_results.md), [fresh native checks](validation/single_market_intercept_validation_results.md) | `single-market-intercept|market-intercept|s0-pooling|conditional-k` |
| Estimation, controls, priors | [estimation](../R/estimate.R), [priors](../R/priors.R) | `stan-control-policy`, `hierarchical-sd-priors`, `market-intercept-sd-prior`, or `residual-sd-prior` |
| Sampler backends and paired validation | [backend adapters](../R/sampler_backend.R), [closeout summary](validation/nutpie_closeout_summary.md), [closeout protocol](validation/nutpie_closeout_protocol.md), [original paired results](validation/nutpie_backend_results.md), [evidence contract](validation/nutpie_evidence_design.md) | `sampler-backend|nutpie|backend-policy-equivalence|evidence-availability` |
| Quality/intercept reconstruction | [quality](../R/posterior_quality.R), [collapse contract](market_intercept_collapse_plan.md) | `posterior-quality|market-intercept-collapse` |
| Evidence and exact-k comparison | [comparison](../R/bma.R), [model universe](../R/model_universe.R), [conditional k](../R/conditional_k.R) | `bma|model-universe|restriction-comparison|conditional-k` |
| Policy states and counterfactuals | [draw plans](../R/posterior_draw_contract.R), [usage](statistical_toolkit_usage.md), [tariff contract](tariff_posterior_reuse_plan.md) | `posterior-policy-models|latent-quality-simulation|simulate-dispatch` |
| Conduct and cross-package models | [conduct](../R/conduct.R), [Stackelberg map](stackelberg_implementation_map.md) | `conduct|stackelberg-r-contract` |

`antitrust` owns common structural/counterfactual machinery; `coordination` owns
fringe/Stackelberg models. Read [DESCRIPTION](../DESCRIPTION) for dependency
minima and [CI](../.github/workflows/qa.yml) for the current development pin and
canonical R version. CI now pins a verified compatible `antitrust` `0.99.42`,
`coordination` `0.0.0.9006` and `trade` `0.8.9.9002` stack, with the validated
CmdStanR revision. The former `antitrust` `0.99.39` pin was incompatible with
current coordination. The ordinary tier explicitly excludes Nutpie/BridgeStan;
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

Some integrations require a newer `coordination` than the declared minimum;
record skips and the resolved dependency versions. Routine checks can parse Stan
when CmdStan exists. Sampling is separately opt-in; short chains are smoke tests.
See the [test audit](test_suite_audit.md) for tier definitions and dated evidence.
Extended deterministic enterprise-fixture checks must also use current source:

```sh
Rscript --vanilla -e 'pkgload::load_all("."); testthat::test_file("doc/market_intercept_enterprise/fixture_checks.R")'
```

For statistical recovery, use the explicit prerequisites, output configuration,
and gates in [current CmdStanR validation](validation/README.md). The archived
124-task [enterprise runner](market_intercept_enterprise/README.md) needs a
reviewed backend port; its passing fixture checks do not validate current HMC.
For documentation-only work, resolve changed local links and run
`git diff --check`; do not rerun sampling or package installation.

## Consolidated agenda and topic evidence

Select follow-up work from the user's task. These statuses are evidence bounds,
not authorization to launch unfinished experiments. They replace the operational
checklists in the historical umbrella and recovery notes below.

| Topic | Recorded status and next decision | Detailed source |
| --- | --- | --- |
| Market/outside-share rotation | Implemented experimental coordinates with bounded sampler evidence; scaling and broader equivalence remain incomplete. | [rotation results](market_s0_rotation_results.md) |
| Latent-quality policies | Implemented; frozen observed/latent benchmark completed. Broader statistical coverage is unestablished. | [policy results](latent_quality_policy_results.md), [usage](statistical_toolkit_usage.md) |
| Hybrid margins | Frozen seven-fixture benchmark completed; outlier/predictive evidence remains pending. | [hybrid results](collapsed_gaussian_hybrid_results.md) |
| Model universe | Public descriptors and evidence aggregation implemented; restricted-family policy BMA is outside that interface. | [usage](statistical_toolkit_usage.md), [design contract](model_universe_specification.md) |
| Hierarchy priors and recovery | Approved defaults and bounded comparisons recorded; existing fits retain frozen priors. | [default history](default_model_history.md), [paired results](validation/hierarchical_prior_recovery_results.md) |
| Market-intercept collapse | Implemented opt-in; fixed-target validation differs from broad sampler certification. Enterprise port remains required. | [validation](market_intercept_collapse_validation.md), [enterprise status](market_intercept_enterprise/README.md) |
| Tariff posterior reuse | Implemented initial route subject to provenance/equivalence conditions; heterogeneous post tariffs remain a separate extension. | [results](tariff_posterior_reuse_results.md), [contract](tariff_posterior_reuse_plan.md) |
| Sampler presets and diagnostic verdicts | Proposed, unimplemented; current controls remain explicit. | [proposal](sampler_usability_followup_plan.md) |
| Prior-predictive diagnostics | Historical requested follow-up; no general public prior-predictive interface was found in this review. Scope separately before implementation. | [historical roadmap](statistical_toolkit_followup_plan.md), [inference contract](market_prediction_inference_design.md) |
| Recovery and model discrimination | Current bounded CmdStanR gates exist; broad repeated-sample coverage and general conduct identification remain open. | [validation entrypoint](validation/README.md) |
| Single-market Logit normalization | Corrected in both Logit kernels: one direct `Normal(0,1)` inside intercept, no single-market hierarchy SD; multi-market hierarchy, outside-share pooling and CES kernels unchanged. Deterministic 10%/20%/40% support and B/C/CF/CS posterior-policy identities pass. Full ordinary suite: 6,454 assertions; CmdStan: 590; package check: zero errors/warnings, one existing NOTE; hosted ordinary CI passed at `112c448`. Fresh validation: 25/25 tasks, 138 paired target states, 18/18 finite bridge repetitions; evidence differences agree within MC uncertainty, with all three precision caps inconclusive and actual mixing warnings retained. Local revised-paper B M1T1 follow-up: 56/56 fits, 168/168 finite bridge repetitions, 60.067%/59.527% Bertrand-family posterior from 25% prior; merger/tariff results and small explicit cost attrition retained. PNB M1T1: 7/7 fits, 21/21 finite bridge repetitions, complete 200-state policies per model; family probabilities 40.621% MC/37.946% Bertrand/21.433% Cournot from equal thirds. Matched Nutpie PNB rerun: 7/7 fits, 21/21 finite bridges, exact source/data/prior identity, maximum posterior model-probability difference .576 percentage points, 3 divergences versus 9 CmdStan; Nutpie C has weaker mixing (maximum R-hat 1.016, minimum bulk ESS 469). All merger states retained. Nutpie PNB acceptance .999 sensitivity: 7/7 fits, 21/21 finite bridges, zero divergences/depth hits, maximum R-hat 1.00724 and minimum bulk ESS 680; native time rises 52.8%, with stable probabilities/policies. Matched CmdStan .999 rerun: 7/7 fits, 21/21 finite bridges, zero divergences/depth hits, maximum R-hat 1.00663 and minimum bulk ESS 977; matched .999 backend model probabilities differ by at most .558 percentage points and BMA merger means by .0121 percentage points. These bounded runs do not establish repeated-sample recovery or universal sampler superiority. Existing fits require refitting under specification 12; downstream banking runner still needs its explicit specification-11 pins updated before production. | [contract](single_market_utility_normalization.md), [results and local follow-ups](validation/single_market_intercept_results.md), [fresh native checks](validation/single_market_intercept_validation_results.md) |
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
