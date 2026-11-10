# antitrustBayes 0.3.0.9002

- Specification 15 adds `share_residual_model = "demand_shock"`, the default for new Gaussian relative-share Logit/CES fits. Persistent quality and period demand are distinct from independent measurement error; the measurement SD has a configurable Lognormal prior with median .10 and log SD .5. The variance allocation can remain prior-dependent.
- `share_residual_model = "combined"` retains the original statistical model. Absolute-share fits continue using that model. Specification-14 saved fits retain their original meaning and verified frozen-kernel provenance; no saved posterior or evidence is relabeled.
- Version-4 policy plans retain jointly recovered quality, period demand, measurement precision, outside shares and role realizations. Counterfactuals carry quality plus period demand through both equilibria; replicated observations include measurement noise. Explicit observed anchoring remains available.
- New-model inference/evidence must be computed before adopting updated model probabilities. Historical results and raw experimental records are preserved.

# antitrustBayes 0.3.0.9001

- Accept proven LF/CRLF-only Stan source differences for saved-fit evidence while preserving original raw hashes and checking embedded source. New fits also record a versioned canonical source identity; substantive source changes remain incompatible.
- Read package revision provenance from the loaded package metadata rather than the caller's Git repository. Document a complete supplied-state example without a Stan toolchain and separate estimation setup; avoid irrelevant baseline-construction warnings.
- Estimated-fit counterfactuals now default to the fitted posterior demand state (`demand_anchoring = "latent_quality"`) in `simulate()`, `bayes_sim()`, and `posterior_policy_models()`. Reuse a `posterior_policy_draws()` plan to hold parameters, recovered conditional quality, and latent roles fixed across policies. Request `demand_anchoring = "observed"` explicitly to reproduce observed-share calibration.
- Preserve the economic quality-scale alias through selective posterior extraction so nested collapsed-quality recovery works with the total/fraction parameterization.
- Specified fits retain their supplied-parameter/calibrated-market default; they do not have an estimated quality posterior. The estimation model specification remains 14, and this simulation-default change does not invalidate compatible inference.
- Support native simultaneous auction/bargaining posterior states, including input Logit auctions through their cost-based value index. Smooth-cutoff Logit auctions still lack a valid latent value-index reconstruction and now fail explicitly; explicit observed calibration remains available for that historical route. Require the paired `antitrust >= 0.99.43.9000` sequential-auction state fix.

# antitrustBayes 0.3.0.9000

- Add a `calcPriceDelta` method for `BayesSim` objects that computes a per-draw Paasche (default) or Laspeyres price index by delegating each draw to `antitrust::calcPriceDelta` on that draw's terminal structural model, after requiring finite prices, nonnegative quantity weights, and a positive denominator. Invalid draws are recorded as `NA` with a reason rather than dropped. Supports `market`, `party`, `levels`, and a `posterior = TRUE` per-draw mode.

# antitrustBayes 0.3.0

- Make relative shares the specification-14 default for Logit/CES estimation and require refitting older model specifications.
- Validate two bounded relative-share known-truth fixtures and independent baseline/merger FOCs; these checks do not establish general repeated-sample coverage.
- Align posterior policy realization and installation requirements with the coordinated structural packages.

# antitrustBayes 0.2.0.9009

- Specification 14 makes relative shares the default for Logit and CES. Existing
  fits from earlier specifications must be refitted with `bayes_est()` before
  comparison or simulation; changing saved metadata cannot migrate the target.
  Use `share_basis = "absolute"` explicitly when refitting an absolute-share
  analysis, and record the chosen likelihood and priors with the new fit.
- Require the compatible `antitrust >= 0.99.42` and
  `coordination >= 0.0.0.9006` API stack; document CI's frozen source revisions.
- Add a separate, bounded specification-14 known-truth CmdStanR release gate for
  the current relative-share Logit panel and finite-noise CES panel. Retain all
  fits and failures, check structural recovery and active-parameter diagnostics,
  and independently reconstruct baseline and merger first-order conditions.
- Qualify `utils::tail()` in CmdStan diagnostics to remove the package-check NOTE.
