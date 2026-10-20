# Nutpie evidence contract

Written 2026-10-19 before changing evidence dispatch. Baseline: package
0.2.0.9005, Git `fff419be8f4d9096eeff7b7f39dcc5e4481ef8ab`.

## Decision

The optional Nutpie backend samples the existing Stan program through
BridgeStan. Sampling the same posterior does not by itself establish a valid
marginal-likelihood calculation. This implementation therefore withholds
Nutpie bridge evidence and rejects mixed-backend comparisons until the
normalized target and equivalent evidence have been executed and validated.
It does not substitute a custom likelihood, a prior approximation, or a reduced
model/k universe. Generated log-likelihood arrays can support predictive
scoring through the common draw contract where ordinary LOO/WAIC is valid,
including predictive-only comparisons across backends. Mixed-backend normalized
evidence and posterior model probabilities remain prohibited.

## Required density and coordinates

Let `theta = T(q)` be Stan's inverse constraint transformation. A bridge
calculation on the unconstrained coordinates must evaluate

```
log f(q) = log p(y | T(q)) + log p(T(q)) + log |det DT(q)|,
Z = integral f(q) dq.
```

Every likelihood/prior normalization constant and truncation normalizer must
be retained. The Jacobian occurs exactly once. BridgeStan's
`StanModel.log_density(q, propto=False, jacobian=True)` supplies the appropriate
evaluation for a normalized Stan program. Its defaults use `propto=True`, so
relying on the default is insufficient. Nutpie's sampler `lp` is also
insufficient evidence of an absolute normalization convention.

The `bridgesampling::bridge_sampler.matrix()` route is mathematically suitable
when its columns contain **only** the Stan unconstrained coordinates, in the
exact evaluator order, with named bounds `lb = -Inf`, `ub = Inf` and real
parameter types. Constrained draws, transformed parameters, generated
quantities, and `lp__` must never be passed as additional coordinates. Already
unconstrained coordinates must not receive a second bound transformation.
The matrix interface treats rows as one ordered chain for ESS. This matches
the existing CmdStan bridge dispatch, which calls the same matrix method after
`unconstrain_draws(format="matrix")`. The Nutpie adapter retains the iteration
and chain index alongside those rows. Executed validation must check chain
mixing independently; concatenation alone cannot establish adequate precision.

## Existing CmdStan route

`bridgesampling` 1.2.1 calls the saved fit's `unconstrain_draws()` and evaluates
`fit$log_prob(q, jacobian=TRUE)`. Installed CmdStanR 0.9.0's
`include/model_methods.cpp` calls `stan::model::log_prob_propto<true>` for
that evaluation. The present four package kernels explicitly accumulate
normalized `_lpdf`/`_lpmf` densities and support normalizers rather than relying
on sampling statements to retain constants. Consequently the CmdStanR
`propto` implementation alone does not prove a discrepancy; the actual
compiled targets must be compared. The existing CmdStan bridge route is
preserved by this backend addition.

## Validation needed before enabling Nutpie evidence

1. Preserve source hash, exact Stan data/hash, kernel/specification, active
   priors, model and full ordered-k support, parameter names/order, toolchain
   versions, inference method, chains and draw provenance.
2. Recover unconstrained coordinates from retained draws and verify
   constrain/unconstrain round trips at nontrivial interior points.
3. At identical coordinates compare BridgeStan `propto=False,
   jacobian=True` with the existing CmdStan bridge evaluator, and independently
   compare Jacobian-on/off and constants-on/off differences. Include positive
   scales, bounded parameters, and every kernel/configuration being certified.
4. Execute a normalized tractable reference and same-model bridge comparison
   with uncertainty/repetitions, plus a nontrivial candidate-set comparison.
   Verify save/reload and original support integrity. Passing mocked dispatch
   tests or an algebra identity does not satisfy this requirement.
5. Reuse the current normalized-specification and HMT exclusions. Missing or
   invalid positive-prior evidence withholds all model probabilities; never
   renormalize over successful candidates.

No executed BridgeStan/CmdStan target or bridge-equivalence validation is
claimed by this note. Nutpie evidence remains explicitly unsupported in this
implementation even if a caller inserts an apparent validation flag in saved
metadata. Conditional-k evidence supplied explicitly to the established
complete-support mixer retains its existing semantics.

## Primary API sources inspected

- [BridgeStan Python API](https://brianward.dev/bridgestan/latest/languages/python.html):
  `log_density(theta_unc, propto=True, jacobian=True)` and parameter transforms.
- [Bridge sampling dispatch source](https://github.com/quentingronau/bridgesampling/blob/master/R/bridge_sampler.R):
  CmdStanMCMC, matrix routes, bounds, and matrix chain assumptions.
- [Bridge sampling Stan evaluator source](https://github.com/quentingronau/bridgesampling/blob/master/R/bridge_sampler_tools.R):
  `.cmdstan_log_posterior` evaluates `log_prob(..., jacobian=TRUE)`.
- [CmdStanR log-probability API](https://mc-stan.org/cmdstanr/reference/fit-method-log_prob.html)
  and installed CmdStanR 0.9.0 `include/model_methods.cpp`: target invocation.

Source/API review is design evidence; executed validation is reported
separately and must identify actual versions and artifacts.

## Implemented capability and local checks

The optional `.bayes_nutpie_bridge_context()` now restores an evaluator from
the existing Nutpie shared-library path, verifies its binary/source/data and
serialized-JSON hashes, and uses `bridgestan::StanModel$new()` without invoking
compilation. It matches only the model's parameter-block names, reconstructs
unconstrained coordinates, checks constraint round trips, retains iteration
and chain indices, and supplies normalized evaluations with unbounded matrix
coordinates. The public bridge route is guarded by an internal false constant;
saved fit metadata cannot override it. Missing optional dependencies, stale
native artifacts or failed integrity checks stop restoration. The restored
native model's BridgeStan version must exactly match the R client version; a
mismatch fails rather than relying on a compatibility warning.

NutpieR's inspected current Rust interface has `bs_param_unconstrain()` and
`bs_evaluate()` but the latter fixes both `propto` and `jacobian` to true. The
separate BridgeStan R interface permits the required explicit flags. NutpieR's
`Cargo.lock` pins BridgeStan 2.8.0; a client using that native library should use
the matching `bridgestan` R package version, rather than assuming its current
main (2.9.0) has a matching ABI.

Focused no-HMC tests passed with the existing compatible paper library:
`nutpie-evidence`, `bridge-context-restoration`, `bridge-performance-migration`,
`bma`, `model-universe`, `conditional-k`, and `restriction-comparison`. The new tests exercise chain/unit identity, source
integrity, unvalidated/mixed-backend rejection, data tampering, and independent
positive/bounded parameter Jacobians and coordinate ordering with a small
evaluator double. They establish dispatch/integrity behavior and the specified
evaluator flags; they are not executed native target/evidence equivalence.

## Bounded executed-validation route

[`nutpie_evidence_validation.R`](nutpie_evidence_validation.R) defines an
explicitly opt-in protocol. Sourcing it performs no inference. Its default
command prepares/collects the complete 15-task design. Native targets reuse
the benchmark's B, CS, CES_B and CES_CS CmdStan/Nutpie-diagonal saved pairs,
covering all four kernels. It evaluates identical unconstrained points with
both Jacobian settings and both BridgeStan constant settings, including
posterior and deterministic interior points. Absolute/relative point tolerances
are fixed at `1e-7 + 1e-9 * max(1, abs(targets))`.

The MC/B/C benchmark fixtures have different observations and must never be
combined into one evidence universe. The evidence runner therefore reuses the
B pair and plans exactly four additional MC/C fits on **B's unchanged data and
priors**, using four chains, two cores, 1000 warmup plus 1000 retained draws,
seed 419771, target acceptance .95 and depth 13. These fits require a separate
opt-in flag and can start only after all four target checks pass. A failed
or timed-out attempted task is retained. Explicit retry archives the complete
previous attempt before beginning another, and completed tasks are reused.

After all shared-data candidates exist, the runner uses three normal bridge
repetitions per model/backend. It reports the larger of repetition SD and
available within-repetition MCSE as numerical uncertainty. Paired log evidence
must agree within four combined uncertainty estimates, and that tolerance
must itself be at most .25 log units. The explicit atomic B/MC/C universe has
equal validation priors, the same observation contract and the complete three
models on both backends. Probability differences must satisfy a conservative
softmax bound capped at .05. Evidence diagnostics require four complete chains,
no divergences/depth hits, BFMI at least .3, R-hat at most 1.01, and bulk/tail
ESS at least 100. These are protocol checks for this bounded exercise, not a
general guarantee of posterior accuracy. A normalized Exp(1) times Uniform
reference separately checks bridge integration and both transformation
Jacobians using exact iid draws without a new Stan kernel or HMC.

Each execution records a new durable invocation clock, capped at four hours
and adjustable with `AB_NUTPIE_EVIDENCE_SECONDS`. The historical backend clock
is preserved and does not block explicit resumption after a machine reset.
Scientific fixture and source manifests remain immutable. The runner retains
native fit hashes, source/data/support identities, R/Stan source hashes, exact coordinate
and evaluator settings, three repetition objects, diagnostics, and partial
failures under `.scratch/nutpie-validation/evidence`. Tracked CSVs contain
compact status, target differences, bridge differences and universe results.
Missing evidence never causes normalization over a subset. A protocol `passed`
result does not change the package gate automatically, and mixed-backend
comparisons remain explicitly unsupported.

With the backend saved pairs complete:

```sh
# Prepare/collect only: no compilation or inference.
Rscript --vanilla doc/validation/nutpie_evidence_validation.R

# Saved-fit target/bridge work only; no new HMC.
AB_NUTPIE_EVIDENCE_EXECUTE=true Rscript --vanilla doc/validation/nutpie_evidence_validation.R

# Opt in to the four additional shared-data fits after target checks pass:
AB_NUTPIE_EVIDENCE_EXECUTE=true AB_NUTPIE_EVIDENCE_SAMPLE=true \
  Rscript --vanilla doc/validation/nutpie_evidence_validation.R

# Explicitly retry interrupted/failed tasks; preserve their previous directories:
AB_NUTPIE_EVIDENCE_EXECUTE=true AB_NUTPIE_EVIDENCE_RETRY=true \
  Rscript --vanilla doc/validation/nutpie_evidence_validation.R
```

Select the existing compatible libraries (including nutpieR 1.9.1 and
bridgestan 2.8.0) and CmdStan location before running. Native compiled caches
are reused from the backend benchmark root. No executed result from this
protocol is claimed until its artifacts exist and have been reviewed.

## Executed analytic reference, 2026-10-19

The normalized Exp(1) times Uniform(0,1) reference completed using
`bridgesampling` 1.2.1, 4000 exact iid draws, and three normal bridge
repetitions. The median log evidence was `0.005866432` (exact value zero),
repetition SD `0.002834759`, and maximum within-repetition log-evidence MCSE
`0.006304421`; the prespecified four-MCSE tolerance was `0.02521768`.
The three estimates were `0.005866432`, `0.006985260`, and `0.001612453`.
Artifacts are retained under `.scratch/nutpie-evidence-audit/reference`.
This checks the integration API and the two transformation Jacobians. It does
not establish native Stan target or economic-model evidence equivalence.
