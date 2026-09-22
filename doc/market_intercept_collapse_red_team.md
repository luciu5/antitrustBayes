# Independent red-team review: market-intercept collapse

The reviewer did not author the production implementation. The review covered
the four Stan kernels, R data/API plumbing, posterior reconstruction, exact-k
generated quantities, model-universe identity, bridge contracts, and latent
policy plans.

## Findings resolved during review

1. The first integration exposed CES as supported before its Stan path was
   complete. The CES simultaneous and exact-k kernels were then implemented,
   parsed, and checked against the same dense oracle as Logit.
2. The first model-universe fingerprint treated sampled and integrated
   intercept coordinates as different scientific hypotheses. The coordinate
   flag was removed from the statistical-hypothesis fingerprint. It remains in
   fit/cache provenance and policy-plan identity, where it belongs.
3. CES initially exposed a zero-valued `b_market` vector in integrated mode.
   It now uses a zero-length sampled-intercept vector and a separately named
   integrated intercept location, matching the Logit output contract.

## Final assessment

No economic-target blocker was found.

- Each market integrates its common intercept and centered quality once.
- Exact-k component scores add the market's integrated joint/share likelihood
  ratio before the market-level `log_sum_exp`.
- Missing-margin observations retain the unconditional share variance.
- A one-group market fixes centered quality to zero.
- Market-specific residual scales use the corresponding market scale.
- Total/fraction quality-scale coordinates remain restricted to a common
  residual scale.
- Posterior reconstruction draws the intercept first and quality conditional
  on the same intercept and realized k.
- The policy score restores the intercept offset once. Its later cancellation
  under normalized inside shares and fixed fitted outside share is the existing
  latent-quality anchoring convention.
- Sampled and integrated coordinates are one model-universe hypothesis, while
  policy plans fingerprint their reconstructed latent state.

The reviewer did not compile Stan, run HMC, or execute the full package suite.
Those gates were assigned to the lead and are reported separately.
