# Nutpie public-API closeout protocol

Prepared implementation, 2026-10-20. No closeout inference or native
recompilation has run. Execution awaits the reviewed production/API freeze.

The bounded protocol has 42 checkpointed tasks:

- Six public B/MC/C comparisons: homogeneous CmdStan, homogeneous Nutpie and
  mixed backends, through ordinary `bayes_compare()` and the explicit universe
  route. Every result also exercises public `bayes_weights()`.
- Four saved fixed-k comparisons: CF and CS, k=2, both backends.
- Eight new CF singleton fits covering k=0,1,3,4 on both backends; original
  k=2 fits are reused. Four chains, 1,000 warmup plus 1,000 retained draws per
  chain, target acceptance .95, depth13, seed `420701+k`, two cores. These are
  new explicitly labeled component runs, with the original target otherwise
  unchanged.
- Eight new public singleton evidence calculations, plus complete conditional
  evidence aggregation. Every public bridge calculation requests three normal
  repetitions; this is 90 planned model/repetition calculations across the
  public routes and singleton checks.
- Ten full 4,000-state component policy calculations and one conditional
  policy mixture calculation, using the existing fixed-plan merger/tariff
  helper and preserving every failed realization.
- Four portability checks: one saved fit per maintained kernel, missing native
  library path, a new empty cache and a genuine recompilation, followed by
  constrained round-trip and normalized-target comparisons at 23 states.

The CF family's complete feasible J=4 support is **0:4**, including both
endpoints. Its declared original prior is (.10,.15,.25,.30,.20). Component
fitting uses singleton prior1; `bayes_combine_conditional_k()` restores the
original prior once. A separate arithmetic identity verifies the normalization.
An explicit mixed-backend combination uses Nutpie at even k and CmdStan at odd
k. No unavailable positive-prior component can be discarded.

Log-evidence agreement uses four combined numerical uncertainties. Each
uncertainty is the larger of bridge repetition SD and available within-repeat
log-evidence MCSE; the coefficient-of-variation delta estimate from
`bridgesampling::error_measures()` is also retained conservatively. Precision
claims require all three repetitions and a tolerance at most .25 log units.
Finite evidence and sampling-quality warnings are separate results: divergence
counts are reported and do not suppress a calculation.

Every public comparison task uses its own predeclared R seed:
`420811 + 1000 * match(task_id, frozen_task_ids)`. All 18 public calls therefore
use distinct streams across components, backends and comparison routes. Models
and bridge repetitions within one public call advance its RNG sequentially.
No shared-seed common-random-number comparison is used; the family uncertainty
propagation assumes independent bridge numerical errors across its separate
component calculations. Exact task seeds remain in the evidence CSV and raw
public return wrappers.
Mixed-versus-CmdStan comparisons share some component fits and bridge objects.
Their reported difference-error bound adds the two source uncertainties,
retaining unknown covariance conservatively. Homogeneous CmdStan-versus-Nutpie
comparisons use root-sum-square errors from their distinct fitting/bridge
streams. The declared precision caps apply to both routes.

Policy means, SDs, medians, central intervals and fixed threshold probabilities
use weighted empirical component distributions with evidence-based k weights.
The original k weights remain unchanged after within-component policy failures;
the resulting policy distribution is explicitly conditional within each k.
Separate failure mass is `sum_k posterior_k * failed_k / planned_k`. A positive
weight component with no valid states makes the mixture unavailable.

Approximate mixture Monte Carlo errors use 200 fixed-seed circular moving-block
bootstrap repetitions, block length32 separately within every retained chain.
NA failure states are resampled with their chain positions. Independent normal
perturbations of component log evidence propagate the retained bridge numerical
uncertainty into k weights separately. The reported total MCSE bound is state
SE plus evidence SE, retaining their unknown dependence conservatively rather
than assuming a zero covariance. Both contributions remain in the results.
This is an uncertainty approximation, not an exact
posterior uncertainty distribution or evidence-certification claim. Exact
constant indicators are identified separately and do not acquire fabricated
MCSEs; matching constants are exempt exact identities, while different constants
fail agreement. Four combined Monte Carlo errors must fit predeclared caps before a
precision claim: .25 price percentage points for means/quantiles, 10% of
reference policy SD for SDs, and .03 probability for thresholds.

Run after freeze:

```bash
bash doc/validation/run_nutpie_closeout.sh --prepare
bash doc/validation/run_nutpie_closeout.sh --run
```

The default invocation budget is four hours and every attempt is retained in
`.scratch/nutpie-closeout`. Repeating `--run` reuses completed tasks and saved
fits. Failed tasks require explicit `AB_NUTPIE_CLOSEOUT_RETRY_FAILED=true`;
previous attempt directories remain intact. No replacement seed or automatic
tuning is provided. Original fit/manifest byte hashes, current production R
sources, helper/runner scripts and dependency provenance are sealed. Preparation
and collection do not sample or compile.
Ordinary tasks reuse the original validation's workspace compilation caches;
the existing backend helpers verify source/compiler/cache keys before reuse.
Resolved cache paths and the CmdStan version/path are sealed in the runtime
contract. Each portability task overrides the Nutpie cache with its own new
empty attempt directory and must genuinely rebuild its maintained kernel.
The new manifest records the starting Git commit and branch separately from
runtime identity, allowing later report-only commits without invalidating fits.
It also records all four Stan source hashes directly and rechecks them before
every worker; the original manifest and fit seals remain unchanged.
