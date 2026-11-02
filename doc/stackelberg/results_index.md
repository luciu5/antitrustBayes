> **Historical checkpoint: September 14, 2026 evidence index.**
> For current development and the consolidated agenda, use the [documentation index](../README.md).
> Empirical pilot progress and package-check counts below belong to the recorded run. They are not current process status or current-suite results.

# Stackelberg evidence index

This index distinguishes derivation, executed numerical validation, posterior
estimation, and unresolved research questions. Historical failed checks remain
in the artifact directory; the final gates identified here supersede them.

## Analytically established

[The mathematical audit](../stackelberg_math_audit.md) derives multiproduct
Logit and CES Bertrand/Cournot follower FOCs, implicit reaction derivatives,
and simultaneous noncooperative leader FOCs. Logit has common dollar markups
within firms; CES has common Lerner margins. Shares in the FOCs are
unconditional; firm sectors count each firm once. The naive fractional-role
interpolation is rejected. Exact ordered partitions are marginalized jointly
over all observations belonging to their equilibrium unit.

[The implementation map](../stackelberg_implementation_map.md),
[ordered-k addendum](ordered_k_addendum.md), and
[nonempty-core addendum](nonempty_core_addendum.md) specify ownership, timing,
Bayesian positive-k support, cutoff conditioning, and posterior simulation.
Native specification 7 retains prior constants and support normalizers.
The inherited optional HMT penalty is not a normalized prior, so bridge
evidence rejects fits using that active penalty; the recovery and PNB
protocols disable it.
No external marginal-likelihood correction or legacy evidence mode remains.

## Numerically verified

- Independent primitive audit: 42 layouts and 26 boundary cases, including
  input-market Logit and multiproduct leaders/followers. See
  `derivative_validation.csv`, `optimization_validation.csv`, and
  `boundary_validation.csv`. The grouped maxima are in
  [derivative_validation_summary.csv](derivative_validation_summary.csv).
- Deterministic final API: baseline calibration, off-equilibrium reaction
  derivatives, ordinary-conduct nesting, and ownership/cost/combined shocks.
  See `coordination_final_*_validation.csv` and the independent `redteam/`
  reports. Both final source-package checks report **Status: OK**.
- Native Stan: 180 prior/flag cases, 72 target comparisons, 96 whole-unit
  density oracles, 20 gradient checks, 24 normalized nesting cases, 12
  singleton cases, and 24 compiled posterior-reconstruction cases.
  See `native_prior_flag_normalization_validation.csv`,
  `compiled_likelihood_reference_unit_validation.csv`, and the compiled
  validation scripts/tables.
- Final Bayesian source tests: **901 passed, 0 failed, 0 skipped**, with
  12 expected existing test warnings. Package check has no errors, warnings,
  or notes; complete output is in `native_release_source_tests.Rout` and
  `native_release_source_check.log`.
- Final deterministic package suite: 383 passed, 24 expected warnings,
  no failures/errors/skips. Its independent package check has no errors,
  warnings, or notes.

## Empirical work

The independent final runner-control review passed against frozen commit
`7f7fc7a`. Repeated-seed recovery started at 06:26 UTC on 2026-09-14. The initial
short run was stopped at the sampling gate: all ten completed fits failed
diagnostics; six were interrupted and 128 had not started. The retained
[initial failure record](initial_sampling_failure/STOP_RECORD.md) separates
these outcomes. The separately labeled [computational rerun](computational_rerun_plan.md)
keeps the targets, priors, data and gates fixed. Its frozen 12-fit pilot
started at 06:53 UTC; two completed fits currently fail Rhat and the other
fits remain in progress or planned. PNB has not started. No posterior recovery or model discrimination is
claimed from specified-draw, smoke, or mathematical tests.
The fixed pre-fit protocol is in [recovery_spec.md](recovery_spec.md).

## Still unresolved

The audit establishes interior FOCs and local reaction regularity, not global
leader-equilibrium existence, uniqueness, or optimality for arbitrary data.
Boundary or inadmissible solver branches remain failures. Initial CES support
is output markets with positive outside revenue and positive recovered costs.
Posterior reconstruction anchors demand to observed equilibrium shares; it
does not propagate the complete noisy demand-quality state. Threshold and
adjacent-k sensitivities recover costs under alternative hard games, whereas
the default mixture uses conditional posterior game probabilities.
Empirical identification, HMC performance, model discrimination, and policy
coverage must be assessed from the executed repeated-seed runs.
