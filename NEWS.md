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
