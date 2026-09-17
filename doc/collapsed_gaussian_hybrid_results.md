# Collapsed Gaussian shares with fixed-t5 margins

Executed 2026-09-16 against baseline `12b26c0`, legacy market/outside-share
coordinates. This experiment changes only the margin likelihood relative to
the original collapsed Gaussian benchmark B. Priors, seven realized fixtures,
original B-derived sampler seeds, and HMC budgets are unchanged. Rho is off.
The manifest was frozen before execution; complete evidence is in the ignored
`doc/statistical_toolkit_runs/20260916/hybrid/` archive.

## Sampler evidence

All seven D fits pass the original usability gate: zero divergences and depth
hits, BFMI >= 0.3, maximum rank R-hat <= 1.01, and minimum monitored bulk/tail
ESS >= 100. The minimum BFMI across fits is 0.671; maximum R-hat is 1.0075.
One Cournot tail-ESS warning is retained in the archive. It does not invalidate
the separately computed gate result.

Root independently recomputed paired efficiency ratios from the delivered
tables. Median native Stan runtime D/B is 1.288; median structural bulk
ESS/native-Stan-second D/B is 0.814. Each design has one D fit, with two
prespecified noisy Bertrand data seeds among the seven designs. These are
controlled small benchmarks, not a broad Monte Carlo recovery study.

## Statistical and economic evidence

The seven designs cover Logit simultaneous Bertrand/Cournot, exact BF and BS,
exact-share stress, noisy shares, and CES Bertrand. Gaussian-share collapse
remains the same mathematically verified likelihood. Fixed-t5 margins change
the statistical model; B/D posterior equality is neither expected nor claimed.

Three bridge repetitions completed for the noisy seed-271828 simultaneous
Bertrand fit: median log evidence 25.3884, range 0.0226. This is evidence for
that D specification, not proof of equivalence to B or a full conduct-model
comparison. The fully normalized density remains available.

## Pending evidence and recommendation

The frozen policy-counterfactual comparison is now complete; see
`doc/latent_quality_policy_results.md` and its compact CSVs. All 84 observed/latent
B/D scenario groups completed and all 8,400 post-policy states passed the FOC
audit. B/D results show similar broad patterns on these fixtures, including
increased uncertainty under latent-quality simulation. A prespecified
margin-outlier experiment remains pending. The original noisy-margin DGP is Gaussian, so this benchmark
does not establish the hybrid's robustness advantage with real outliers.

The hybrid is computationally viable on these fixtures and merits a prominent
robustness option. No default is changed. A routine-default decision still needs predictive and outlier evidence and
assessment of the latent-quality normalization. Saved-fit bridge contexts also
require the separately validated reconstruction adapter after a session restart.
