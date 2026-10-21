# Nutpie public-API closeout

See the [final closeout summary](nutpie_closeout_summary.md) for interpretation.
The frozen policy estimates and original approximate uncertainties below remain
preserved. Use the [separate conservative uncertainty audit](nutpie_closeout_uncertainty_results.md)
for the final policy precision assessment: 13/28 comparisons remain inconclusive.

Collected 2026-10-21 02:14:18 UTC.

42 of 42 tasks completed; 0 failed.

Execution completion, finite evidence, numerical precision, sampler diagnostics and policy support are separate outcomes.
No divergence threshold gates calculations. Original targets, priors, controls and failed states remain recorded.

Protocol: J=4 CF full support 0:4, prior (.10,.15,.25,.30,.20), eight new singleton fits; original k=2 fits reused.
Each public comparison requests three normal bridge repetitions. Agreement uses four combined numerical uncertainties, capped at .25 log units for precision claims.
Policy calculations retain all 4,000 draws/component with fixed plans; conditional mixtures retain original evidence-based k weights after within-component omissions.
Mixture MCSE bounds add approximate state block-bootstrap and bridge-weight SEs without assuming their dependence is zero; failures are retained and resampled, not replaced.

See nutpie_closeout_status.csv, evidence/diagnostics/comparisons/probabilities, conditional_k, component_policy/policy_failures, policy/policy_comparisons and portability CSVs beside this report.
Raw public return objects, repetition objects, plans, draw failures, source/fit seals and every attempt remain in .scratch/nutpie-closeout.
