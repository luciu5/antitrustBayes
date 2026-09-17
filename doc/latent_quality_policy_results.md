# Frozen latent-quality policy benchmark results

This is the review report for archive `latent-quality-policy-attempt002`. The
run used the corrected policy runner at
`doc/statistical_toolkit_runs/20260916/latent-quality-policy/run_latent_quality_policy_benchmark.R`
with `antitrustBayes` 0.2.0, `antitrust` 0.99.39, and `coordination`
0.0.0.9003 from `/tmp/toolkit-reviewed-Rlib:/tmp/cga-Rlib`. The runner SHA256
was `9dd4a9f162155fadf733fd8ab524e283f64fa61b806be1b37fe0c369fdca202a`.

The benchmark completed 84 scenario groups: 14 frozen B/D cases, observed
and latent-quality anchors, three policies (`merger`, `cost10`, and
`merger_cost10`), and 100 selected posterior draws per group. It evaluated
8,400 draw scenarios and wrote 478,800 metric rows, 4,788 summary rows, and
239,400 observed/latent paired rows. Every scenario group was complete. Every
draw had finite pre and post prices, with zero solver failures, zero partial
failures, and zero top-level simulation errors. Every summary denominator is
100.

## Frozen draw identity and pairing

The regenerated attempt002 plans matched the preserved attempt001 plans for
all 14 cases: selected draw indices, quality matrices, k caches, group labels,
and role fields were each exactly identical. The audit is in
`latent_quality_policy_plan_identity.csv`; source and package hashes are in
`plan_identity_source_hashes.csv` and `policy_results_provenance.csv`.

B and D remain separate posterior cases throughout the output keys. There are
14 unique `case_id`/`mode` combinations and 14 unique design/mode combinations.
Although the deterministic draw schedule uses the same positions for B and D,
their quality vectors are different; across the seven designs the maximum
absolute B/D quality difference ranges from 0.3138303 to 0.4368509. No B/D
rows were pooled as if they were the same theta draw.

The selected outcome metrics have 21,000/21,000 joint observed/latent paired
draws. The all-metric paired file has 8,400 expected joint failures, all from
the two latent-pre diagnostic metrics in the observed anchor, where those
diagnostics are unavailable by construction. The corresponding 84 summary
rows are explicitly marked `metric_unavailable` with requested denominator
100 and valid denominator 0; they are not solver failures.

## Paired latent-quality shifts

The comparison retains each design's units. The reviewable files are
`latent_quality_policy_comparison.csv`,
`latent_quality_policy_paired_shifts.csv`, and
`latent_quality_merger_price_comparison.csv`. The last reports the
fixed-reference-weight price-index change, a fraction rather than percentage
points, for each case and anchor.

For Gaussian-margin fits (B), adding latent relative-quality draws changes the
mean merger effect by +3.8% in the first noisy Logit seed, -11.6% in the second,
and -20.1% in CES Bertrand. The corresponding posterior standard deviations
increase by factors of 4.5, 4.2, and 6.3. The two Logit seeds remain separate
replications. Fringe and Stackelberg fixtures also show wider distributions:
SD factors 1.5 and 2.6 respectively. The exact-share stress fixture increases
SD by a factor of 4.3.

These percentages describe changes relative to small mean merger effects;
they are not percentage-point increases in prices. For example, the first
noisy Logit mean price-index effect changes from 0.0533% to 0.0553%, while the
CES mean changes from 0.1887% to 0.1508%. Latent quality matters particularly
for uncertainty, even where the mean effect remains close.

Hybrid-margin fits (D) show the same broad pattern in these fixtures. B and D
are separately fitted posteriors, so their draw positions do not constitute
paired theta draws. The synthetic margins have no prespecified outlier
contamination; similarity here does not establish robust performance under
accounting outliers. No default is changed by this evidence.

## Latent pre-equilibrium checks

The latent anchor solved a new pre-policy equilibrium while retaining the
observed-anchor cost primitive. Across case/policy summaries, the mean of
`latent_pre_foc_max` ranges from `2.447764e-15` to `5.671635e-9`; its 95th
percentile ranges from `3.268219e-15` to `1.534198e-8`. The corresponding mean
maximum absolute pre-price discrepancy ranges from `2.742386e-3` to
`4.573858e-2`, with 95th percentiles from `4.016192e-3` to `5.551295e-2`.
The full per-case summaries are in `latent_quality_policy_comparison.csv`.

The nonzero price discrepancy is expected under the retained-cost
specification: demand is replaced by the latent-quality state and the
pre/post equilibria are solved at the original cost primitive. Costs are not
re-inverted after the demand swap.

An independent read-only audit then loaded all 84 stored scenario objects and
their 100 draws (8,400 post states). It evaluated the generic output residual
`pricePost - mcPost - calcMargins(level=TRUE, preMerger=FALSE)` using each
object's post-dispatched ownership and native post-role fields. All 8,400
states passed the `1e-6` threshold, with zero method errors or nonfinite
residuals. The largest absolute residual was `3.879225e-7` (LogitCournot,
latent-quality, `cost10`); the largest class-specific residuals were
`1.488649e-8` for Logit, `1.488084e-8` for CES, `4.718448e-15` for native
core-fringe Logit, and `8.679280e-12` for native Stackelberg Logit. The full
per-draw and compact summary audits are
`post_foc_diagnostics_per_draw.csv` in the ignored attempt002 archive and
`latent_quality_policy_post_foc.csv` in the reviewable doc directory.

## Normalization and remaining limits

The latent construction uses the documented outside-share-preserving softmax
projection. It reanchors the common utility level and propagates relative
quality uncertainty; it does not claim a fully generative compositional share
likelihood or full absolute-quality-level propagation. CES quantity shares and
revenue conversions remain distinct, and the fixed market size/budget
convention is retained.

The benchmark covers output-demand Logit/CES and the frozen B/D fixtures.
It is not input-procurement HMC evidence or evidence for a new default.
A separate narrow input-Logit extension preserves the native value primitive
and positive markdown, with B/C and native-role unit/reference checks. Its
independent final review is pending; input CES remains unsupported.

## Audit artifacts

`latent_quality_policy_audit.csv` records the row counts, denominators, statuses, and
pairing checks. `policy_results_provenance.csv` records the runner hash,
package versions, library paths, and archive. Raw result RDS files are in the
temporary run directory referenced by `policy_progress.csv`; the compact CSVs
are retained in the ignored attempt002 archive and in the review scratch
directory.

The reviewable runner `doc/latent_quality_policy_benchmark.R` is a byte-identical
copy of the executed runner. Original failed preparation attempts and raw
objects remain locally archived and excluded from Git.
