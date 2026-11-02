# Latent-quality policy comparison: frozen protocol

This is a substantive counterfactual mapping experiment, separate from the
exact market/outside-share coordinate experiment. Neither estimation nor a
simulation default changes. Freeze this protocol before executing simulations.

Use the seven original prespecified designs and their existing usable B fits
(collapsed Gaussian shares, Gaussian margins), and the corresponding seven
usable D fits (collapsed Gaussian shares, fixed-t5 margins). Do not refit any
model. A and C fits failed the original HMC gate and are not evidence for a
posterior-policy comparison. Keep their failures in the earlier results.

For each fit, select 100 posterior parameter indices using the original
benchmark seed helper with the design identity and the purpose
`latent_quality_policy_plan`. Construct one `posterior_policy_draws()` plan:
one global conditional quality draw per selected parameter row and one exact
ordered-k realization per market/row, reused across years and interventions.
Save the plan and frozen inputs before the first equilibrium calculation.
B and D are different fitted distributions; shared design seeds do not make
their posterior parameter rows paired draws of the same distribution.

Run both `demand_anchoring="observed"` and `"latent_quality"` with that plan
for each of these three interventions, using one worker and insideSize=1:

1. Merge the two largest observed firms, without a cost change.
2. Increase all marginal costs by 10%, without changing ownership.
3. Combine that merger and the 10% marginal-cost increase.

Firm ranking and merger identities use predetermined observed firm shares.
The latent mode uses the outside-preserving reanchoring defined in
`latent_quality_policy_math.md`. It retains the observed-anchor FOC-inverted
cost primitive and physical Logit market size/CES expenditure budget, then
solves a new pre-policy equilibrium. The same quality and cost primitive enter
post-policy demand; only the declared intervention changes ownership/costs.
This compares a documented normalization of latent relative demand, rather
than claiming the fitted absolute-log-share model is a simplex likelihood.

Record every requested draw and retain all solver errors. Report requested
and valid counts and failure rates for each fit, anchor, and intervention.
Save full simulation objects outside Git. Summaries include product price and
share changes, the fixed observed-quantity-weight price index, inside share,
consumer change (-CV), producer surplus change, and their sum. Report missing
welfare methods as unavailable. For each metric return mean, SD, and 5th,
50th, and 95th percentiles. For observed versus latent comparisons within a
fit, also report paired differences on explicitly labeled jointly successful
draws; do not hide failures by changing denominators. Report the latent
baseline discrepancy and maximum finite FOC residual.

No ESS-per-second inference is made from a randomly selected 100-draw policy
subset. This is a bounded distributional and solver validation, not a recovery
study or an outlier-robustness experiment. No diagnostic gate, prior, HMC
setting, or policy definition may be changed to make a failed comparison pass.

Store heavy artifacts under the ignored
`doc/statistical_toolkit_runs/20260916/latent-quality-policy/` directory. Track
only the runner, compact CSVs, and assessment. Do not change the package default
from these comparisons alone.
