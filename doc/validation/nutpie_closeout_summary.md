# Nutpie closeout: public evidence and mixed backends

Executed 2026-10-20 EDT / 2026-10-21 UTC, package `0.2.0.9006`.
All **42/42 tasks completed**, with zero execution failures. This validates the
scoped feature; the numerical precision limitations below remain inconclusive.
Starting `main`: `fff419be8f4d9096eeff7b7f39dcc5e4481ef8ab`.
Starting `nutpie`: `8d6daa115bf2fb7888a26d1c3ca4ba6666355c64`.
The user fetched both refs during final validation; both remained unchanged.
Production implementation: `c72b1fc`; frozen validation runner: `e8e16f3`;
independently reviewed saved-policy uncertainty audit: `413e6f0`.
Subsequent report commits preserve that production source. The final branch
identities and push result are recorded in the closeout handoff/Git history.

## Implementation

The existing thin backend seam, `BayesLogit`, posterior arrays, generated
quantities and economic reconstruction were retained. This closeout removes
the global `.bayes_nutpie_bridge_validated` flag and its permission checks;
there is no replacement certification flag or diagnostic threshold.

Normalized bridge evidence uses the same Stan target in unconstrained
coordinates, with `propto = FALSE, jacobian = TRUE`. Exact source, saved data,
proper priors/model specification, completed HMC draws, parameter coordinates,
round trips and finite evaluations remain mandatory. Divergences, R-hat and
ESS describe reliability and remain available; they do not suppress evidence.

`bayes_compare()` accepts homogeneous and mixed CmdStan/Nutpie candidates and
records backend, adaptation, controls, diagnostics, source/data hashes,
specification and bridge uncertainty by model. Backend is computational
provenance, not scientific model identity. Existing outcome, likelihood,
market-side and shared-prior contracts remain enforced. Duplicate scientific
hypotheses cannot acquire extra prior mass merely by changing their sampler.

Every requested bridge repetition and its status is retained. Finite successes
are summarized even when another repetition fails; repetition SD requires at
least two successes. Every positive-prior candidate still requires finite
evidence before model probabilities exist. Bridge weights never silently fall
back to stacking or renormalize a surviving subset. Conditional-k requires all
original states, restores their original prior once, and accepts mixed backends
without a new `backend_by_k` interface.

Saved native libraries are hash-verified and loaded through content-named
copies in the external cache, avoiding R DLL basename collisions. Missing,
corrupt or incompatible libraries can be rebuilt using the exact matching
installed source/data and a compatible toolchain. Every restoration checks all
saved sampler `lp__` values. A density disagreement fails without recompilation
rescue; source/data mismatches fail before native work. Restoration provenance
is separate from the original fit, which is not rewritten.

Public API remains `bayes_est(..., backend = "cmdstanr")` by default, or
`backend = "nutpie", nutpie_options = list(adaptation = "diag")`.
`adaptation = "low_rank"` remains an expert option. Common controls remain in
`stan_options`; HMC and CmdStan Pathfinder behavior is preserved. Optional
dependencies remain `nutpieR >= 1.9.1` and BridgeStan; this closeout adds none.
Python/PyMC/flow adaptation and custom density kernels were not introduced.

## Evidence and full conditional-k validation

The [frozen protocol](nutpie_closeout_protocol.md) executes public comparisons,
model universes, fixed CF/CS evidence, eight new singleton fits, complete CF
support `0:4`, all component policies and four empty-cache rebuilds. Original
k=2 fits were reused only after exact source/data/prior identity checks.
Each fit has four chains and 1,000 retained draws per chain. New fits use 1,000
warmup draws, target acceptance `.95`, depth13, two cores and seed `420701+k`.
All eight new fits had zero divergences/depth hits, maximum R-hat `1.006203`,
minimum bulk ESS `983.94` and tail ESS `815.02`. These are reported outcomes,
not availability rules. Reused CF CmdStan/Nutpie fits retain 22/12 divergences;
CS retains 19/7. No settings or fits were replaced to remove warnings.

All 90 requested public bridge repetitions returned finite evidence. Every one
of 21 evidence agreement comparisons is within four combined numerical errors.
Twenty meet the precision cap; the mixed-universe MC comparison has difference
`-0.059928` and conservative four-error bound `.271119`, above the `.25` cap,
so its precision result is inconclusive. The ordinary public mixed comparison
passes all three precision checks. Model probabilities remain complete; their
largest B/MC/C backend difference across routes is `.010634`, below `.05`.

| Case | CmdStan log evidence | Nutpie log evidence | Nutpie minus CmdStan | Four combined errors |
| --- | ---: | ---: | ---: | ---: |
| B | -3.396304 | -3.401170 | -0.004866 | 0.180565 |
| MC | -5.569507 | -5.623483 | -0.053977 | 0.180809 |
| C | -4.762691 | -4.789774 | -0.027083 | 0.120376 |
| CF, k=2 | -5.855591 | -5.867276 | -0.011685 | 0.186205 |
| CS, k=2 | -3.136000 | -3.105209 | +0.030791 | 0.148614 |
| CF, full 0:4 | -7.173344 | -7.185460 | -0.012117 | 0.173971 |

The original CF prior is `(.10,.15,.25,.30,.20)`. Independent log-sum-exp
arithmetic verifies `Z = sum_k prior_k * Z_k` and every posterior-k weight.
The largest CmdStan/Nutpie posterior-k difference is `.000496227`, or
`.000506797` when the mixed route is included. Mixed family log evidence is
`-7.184984`, difference `-.011640`, four-error bound `.245209`: agreement and
precision both pass. No component or prior mass was dropped.

Exact repetition values/SDs, within-repeat uncertainty, diagnostics and weights
are in [evidence](nutpie_closeout_evidence.csv),
[comparisons](nutpie_closeout_comparisons.csv),
[diagnostics](nutpie_closeout_diagnostics.csv) and
[conditional-k](nutpie_closeout_conditional_k.csv). Earlier normalized-target
agreement across all four kernels was approximately `1.1e-13`; its design and
analytic integration check remain in the [original evidence report](nutpie_evidence_results.md).

## Merger and tariff distributions

The run retains 40,000 baseline states and 80,000 policy records; 77,418 records
are valid. All 1,291 baseline failures are classified finite nonpositive costs,
with unchanged support rules. The maximum physical-cost identity error among
valid policy records is exactly zero. Draw plans, quality identity, original
chain positions and all failures remain in the saved results.

The table shows price percent changes, conditional on valid states within each
k, retaining the original evidence-based k weights. It is not an unconditional
policy distribution. Weighted failure mass is `.004153` CmdStan, `.003744`
Nutpie and `.003648` mixed. The small-weight k=1 component still has 484/4,000
CmdStan and 518/4,000 Nutpie failures; this is disclosed rather than hidden by
its small posterior weight. See [failure categories](nutpie_closeout_realization_failures.csv).

| Policy | Backend | Mean | SD | Median | 2.5% / 97.5% | P(effect >5%) |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| Merger | CmdStan | 7.16007 | 1.24166 | 7.09754 | 5.23106 / 9.40755 | .981363 |
| Merger | Nutpie | 7.11348 | 1.25584 | 7.08058 | 5.08703 / 9.37196 | .977534 |
| Merger | Mixed | 7.11534 | 1.26631 | 7.08058 | 5.08703 / 9.37455 | .977487 |
| Tariff | CmdStan | 4.95630 | .51875 | 5.07827 | 3.02027 / 5.37093 | .745202 |
| Tariff | Nutpie | 4.97295 | .52854 | 5.08293 | 3.01168 / 5.44641 | .752590 |
| Tariff | Mixed | 4.97323 | .52708 | 5.08292 | 3.02039 / 5.44701 | .752277 |

The original helper shares its latent-quality seed across components/backends.
A separate saved-row audit therefore sums component-at-a-time block-bootstrap
SEs, adds evidence-weight uncertainty, then sums compared-route bounds without
assuming independent policy errors. This is a first-order uncertainty bound,
not an exact finite-sample quantile confidence interval. No fit, policy value,
prior weight or original estimate changed. The [audit report](nutpie_closeout_uncertainty_results.md)
describes its limitations and links every contribution.

All 28 revised comparisons are compatible with their approximate four-error
bounds; 13 remain precision-inconclusive and two are exact constant identities.
Both means meet the `.25` price-point cap: merger difference `-.046593` with
bound `.184220`; tariff difference `+.016645` with bound `.067686`.
SDs/tails and some threshold probabilities do not have strong precision
validation at this budget. For example, merger lower/upper quantile differences
are `-.144026`/`-.035595`, but bounds `.891789`/`1.106703` exceed `.25`.
These inconclusive results are not relabeled passes or a blanket policy-tail
equivalence claim. Original approximate outputs remain preserved separately.

## Portability, regressions and recommendation

Native restoration results are in [portability](nutpie_closeout_portability.csv),
with fresh compile/cache preparation reported separately from sampling. Every
rebuild checks 23 normalized-target points, all 4,000 saved log densities and
all constrained/unconstrained round trips at tolerance `1e-8`; original fit
bytes and scientific source/data identity must remain unchanged.

| Kernel | Fresh rebuild seconds | Normalized points / maximum difference | Saved log densities / maximum difference |
| --- | ---: | --- | --- |
| logit_supply | 99.013 | 23 / 0 | 4,000 / 1.78e-15 |
| logit_stackelberg | 97.086 | 23 / 0 | 4,000 / 1.78e-15 |
| ces_supply | 84.009 | 23 / 0 | 4,000 / 1.78e-15 |
| ces_stackelberg | 97.629 | 23 / 0 | 4,000 / 1.78e-15 |

All four original fit byte hashes remain unchanged. These are target/library
restoration checks, not new sampling-performance benchmarks.

[Regression checks](nutpie_closeout_regressions.md): 5,166 ordinary assertions,
5,159 without optional native dependencies, and 590 existing opt-in CmdStan
assertions pass. The source package builds/checks with zero errors/warnings and
one pre-existing `tail` namespace NOTE; installed tests pass 4,938 assertions.
All four installed Stan byte hashes match the frozen fits; installed Nutpie
evidence and help/default-route smoke checks pass. Independent source and result
reviews found no substantive scientific/code blocker; final native results
passed and the reviewer recommends the scoped feature merge. Hosted CI and a
clean native CI installation remain unrun.

The [original 34-fit benchmark](nutpie_backend_results.md) and
[stricter nine-fit pilot](nutpie_certification_pilot_results.md) remain intact.
The pooled C/CF/CS benchmark recorded Nutpie diagonal divergences 1/1/2 versus
low-rank 3/13/9; this does not compare against later tuned production CmdStan
profiles. No universal sampler-quality or superiority claim follows. Nutpie's
native `max_energy_error = -Inf` sentinel is not itself a divergence; native
divergence flags remain authoritative. No upstream Nutpie target-correctness
defect was found.

Keep CmdStanR as the default and mature fallback. Nutpie diagonal is a validated
same-target alternative, particularly worth trying for difficult hierarchical
C/CF/CS fits; low-rank remains an expert/experimental choice. Finer policy-tail
precision and broader sampler qualification remain outside this closeout.

No Stan kernel, prior/default, likelihood, conduct equation, policy-role rule,
counterfactual definition or cost-support rule was changed. No `mergerBayes`
production analysis, DGP, sampler profile or scientific specification was changed.
The initial unrelated `.scratch/` work remains untouched and uncommitted.
