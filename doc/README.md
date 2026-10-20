# Development and evidence index

Reviewed 2026-10-20 against baseline `fff419b` and the optional Nutpie backend
work, package version `0.2.0.9005`.
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
| Estimation, controls, priors | [estimation](../R/estimate.R), [priors](../R/priors.R) | `stan-control-policy`, `hierarchical-sd-priors`, `market-intercept-sd-prior`, or `residual-sd-prior` |
| Sampler backends and paired validation | [backend adapters](../R/sampler_backend.R), [Nutpie status/results](validation/nutpie_backend_results.md), [evidence contract](validation/nutpie_evidence_design.md) | `sampler-backend|nutpie|backend-policy-equivalence` |
| Quality/intercept reconstruction | [quality](../R/posterior_quality.R), [collapse contract](market_intercept_collapse_plan.md) | `posterior-quality|market-intercept-collapse` |
| Evidence and exact-k comparison | [comparison](../R/bma.R), [model universe](../R/model_universe.R), [conditional k](../R/conditional_k.R) | `bma|model-universe|restriction-comparison|conditional-k` |
| Policy states and counterfactuals | [draw plans](../R/posterior_draw_contract.R), [usage](statistical_toolkit_usage.md), [tariff contract](tariff_posterior_reuse_plan.md) | `posterior-policy-models|latent-quality-simulation|simulate-dispatch` |
| Conduct and cross-package models | [conduct](../R/conduct.R), [Stackelberg map](stackelberg_implementation_map.md) | `conduct|stackelberg-r-contract` |

`antitrust` owns common structural/counterfactual machinery; `coordination` owns
fringe/Stackelberg models. Read [DESCRIPTION](../DESCRIPTION) for dependency
minima and [CI](../.github/workflows/qa.yml) for the current development pin and
canonical R version. CI pins `antitrust`; it does not pin `coordination`.
The pin is version `0.99.39`, while `coordination` `0.0.0.9006` requires
`antitrust >= 0.99.41`; those revisions cannot be assumed compatible together.
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
| Optional Nutpie sampler | Completed 34 fits plus the unchanged-fit CES policy replay: all 1,385 scored checks within declared MCSE tolerances; diagnostics, attrition and tail uncertainty retained. All 15 evidence tasks passed numerical checks, but divergences leave certification inconclusive and evidence/BMA disabled. CmdStanR remains default; Nutpie experimental. Opt-in CmdStan regressions: 590 assertions passed; focused runner: 61 passed. | [executed results](validation/nutpie_backend_results.md), [evidence results](validation/nutpie_evidence_results.md), [paired runner](validation/nutpie_backend_validation.R) |

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
