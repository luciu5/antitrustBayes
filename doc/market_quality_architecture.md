# Market-quality implementation map

Status: native specification and R integration independently reviewed; the
package suite passes. Sampler and application gates remain incomplete, with
initial sampler failures retained. Defaults are not yet changed.

## Identity and persistence

Economic ownership/ranking IDs (`group`, `N_group`) retain their meaning. A
separate quality incidence map identifies either historical global groups or
market-group units. Units are contiguous by market under market scope, with
groups ordered by their original integer group ID within each market. This
preserves the historical Helmert basis exactly for one market. An explicit
lookup retains original market and group labels alongside a unique quality ID.
Periods do not enter the quality ID.

The affected R entry points are `bayes_est()` and the Logit/CES data preparers
in `R/estimate.R` and `R/helpers.R`; quality recovery lives in
`R/posterior_quality.R`. `R/latent_quality_simulation.R` and the selective draw
contract consume the new incidence map. Model-universe fingerprints,
structural reconstruction metadata, printed fits, and comparison tables carry
scope separately from economic ownership. The four Stan sources are
`logit_supply`, `logit_stackelberg`, `ces_supply`, and `ces_stackelberg`;
their generated C++ headers are regenerated from those sources.

## Stan data contract

The four existing Logit/CES ordinary/exact-k programs receive:

* `quality_scope`: 1 global, 2 market.
* `N_quality`, `obs_to_quality[N]`: separate quality units/incidence.
* `N_quality_block`, `quality_block_start`, `quality_block_n`: contiguous
  centering blocks. Global scope has one block; market scope has one per market.
* `N_quality_coord = N_quality - N_quality_block`, `quality_coord_start`:
  contrast-coordinate offsets. Singleton blocks have no coordinates.
* Complete-market observation and matched-margin CSR indices. Ordinary and
  exact programs use the same maps; data preparation must not duplicate keys.

The existing `b_group_centered` sampling-coordinate name and `b_quality_all`
economic output name remain readable, with dimensions determined by the new
quality maps. Each block uses the already validated linear-time Helmert
transform. Global scope uses the same one-block operation as before.

## Collapsed target and generated quantities

Share statistics are accumulated in one pass using quality incidence and
market CSR. Market scope evaluates independently centered blocks. Exact-k
components add only the local matched margin innovations before normalized
Gaussian integration; their log scores contain joint-minus-share block density.
Each market's share block is added once. No component traverses all markets or
allocates global-size quality statistics.

The number of exact-state evaluations is proportional to the sum of the
market support sizes, not their Cartesian product. Each evaluation still
requires local observation, margin, group and firm-period work. The retained
generated-quantity median calculation also performs quadratic work in local
support size and uses local observation-by-support scratch arrays. These are
not dense Gaussian covariance matrices, but prevent claiming linear cost
independent of a growing support size.

Global scope retains its historical integration and restrictions. In
particular, global collapsed quality plus active rho and multi-market exact k
remains rejected because it does not factorize market by market. Market scope
removes that restriction only when the newly derived factorization applies.

Generated component probabilities use exactly the model's normalized component
scores. Market-level joint log likelihood includes its collapsed share block
and marginalized conditional margin contribution. Observation-wise LOO/WAIC
remains unavailable for collapsed dependent outcomes. Existing demand/conduct
functions and prediction formulas are untouched.

## R reconstruction and policy draws

`posterior_quality()` reads quality incidence and block labels. For active
rho plus exact k it realizes k separately by market, then conditions quality
on that component's margin mean. A shared pure margin-prediction helper must
match the existing Stan equations exactly; it introduces no equilibrium solver.
Outputs identify market, original group, quality unit, structural draw, and
realized k. Conditional means/variances are conditional on the realized state;
mixture summaries must average appropriately over state uncertainty.

`posterior_policy_draws()` realizes roles before quality and passes those exact
role realizations into quality reconstruction. The reusable plan carries scope,
unit lookup, and quality vector with each structural draw. Latent demand
anchoring uses observation-to-quality mapping, and keeps the realized vector
fixed pre/post policy. Selective extraction remains in place.

## Version and defaults

Model specification version advances from 9 to 10. Scope is stored in fit
metadata, Stan attributes, compatibility fingerprints, and printed summaries.
An absent scope in a historical fit means historical global scope; it never
means the new market scope. Unsupported old native targets require an explicit
refit rather than silent migration. New runtime data contracts reject missing
scope/index maps; tiny historical test fixtures may resolve global scope
explicitly through dedicated compatibility helpers.

Before all default-change gates pass, callers select market scope explicitly
and the existing function defaults remain intact. After gates pass, signatures,
help pages, tests, and history change together to Gaussian shares/margins,
collapsed market quality, and active rho. No unrelated pooling, prior, sampler,
or structural default changes are part of this task.

## Supported statistical combinations

These entries describe current implementation support, not a passed sampler
gate. Exact families use one ordered-k state per complete market across
periods; Stackelberg/fringe conduct support is unchanged.

| share / margin likelihood | quality | scope | rho | support |
|---|---|---|---|---|
| Gaussian / Gaussian | sampled | market or global | off or active | ordinary and exact-k |
| Student-t / Student-t | sampled | market or global | off or active | ordinary and exact-k |
| Gaussian / Gaussian | collapsed | market | off or active | ordinary and multi-market exact-k |
| Gaussian / Gaussian | collapsed | global | off | ordinary and exact-k |
| Gaussian / Gaussian | collapsed | global | active | ordinary; exact-k restricted to one market |
| Gaussian / fixed-df Student-t | collapsed | market or global | off | retained hybrid |
| Student-t shares | collapsed | either | either | rejected: Gaussian collapse is not applicable |

Collapsed quality still rejects market-level residual-scale hierarchies.
The hybrid rejects active rho and estimated degrees of freedom. Fully Gaussian
models have no df parameter. `posterior_quality()` reconstructs collapsed
fits; sampled quality remains directly available in Stan's economic outputs.
Auction/bargaining latent-quality counterfactual anchoring remains unavailable;
complete-universe comparisons use common observed anchoring. Observation-wise
LOO/WAIC remains disabled for collapsed dependent outcomes.

Multi-market global collapsed quality with active rho and exact-k is a genuine
factorization restriction. It is not waived for speed; callers select the
distinct market-quality model explicitly if they want market factorization.

## Worker boundaries

Independent reconnaissance and Gaussian oracle precede implementation. Stan
implementation and R API/reconstruction use separate private clones, owning
disjoint files. Root reviews equations and integrates sequentially. A fresh
reviewer checks normalized targets, single-market equivalence, exact mixture
semantics, scope/version handling, and policies before HMC benchmarks.
