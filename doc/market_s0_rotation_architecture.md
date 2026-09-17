# Market/outside-share rotation: implementation map

Baseline: `12b26c0`, statistical specification 9. This is an experimental
coordinate change. The legacy path remains the default. Independent
reconnaissance and native target/gradient validation are complete. The
experimental API and all four coupled Stan kernels are integrated locally;
paired HMC benchmarking is running without a default change.

## Scope and API

Add `market_s0_parameterization = c("legacy", "rotated")` to `bayes_est()`,
`prepare_stan_data()`, and `prepare_ces_stan_data()`. Store the requested mode,
the effective rotation indicator, and the coordinate convention in fit
metadata. Encode legacy/rotated as 1/2 in Stan data. The existing selected
`s0_pooling`, likelihoods, priors, and conduct remain unchanged.

All single-market specifications normalize `b_market` to zero. The requested
rotation is then a no-op, with the original sampling dimensions and names.
For multiple markets, directly sampled outside levels use `c_s0 = v + b`.
Noncentered outside hierarchies use the standardized rotated coordinate
`c_s0_raw = z + b / sigma_s0`. The latter preserves a unit native-coordinate
Jacobian; it is not an independent standard-normal random effect.

Expose reconstructed base locations, unbounded market-year sums, and effective
bounded-link sums as diagnostic transformed parameters. Retain familiar
`b_market`, `s0_prob`, `alpha_structural_market`, and CES revenue-share/elasticity
outputs. Keep a standard `stanfit` object and its ordinary sampler diagnostics.

## File responsibilities

| Files | Authorized change |
|---|---|
| `R/estimate.R` | Argument validation, data forwarding, coordinate metadata |
| `R/helpers.R` | Both demand data-preparation paths and mode resolver |
| `inst/stan/logit_supply.stan` | Direct/pooled coordinate reconstruction and prior evaluation |
| `inst/stan/logit_stackelberg.stan` | Same reconstruction; preserve exact-k mixture |
| `inst/stan/ces_supply.stan` | Same reconstruction; preserve CES normalization |
| `inst/stan/ces_stackelberg.stan` | Same reconstruction and unchanged exact-k mixture |
| `src/stanExports_*` | Regenerate using the normal package build workflow |
| `man/` | Regenerate documentation for the new argument |
| `tests/testthat/test-market-s0-rotation-*` | API, source, and paired-state controls |
| `doc/market_s0_rotation_*` | Independent references, executed validation, benchmarks, assessment |

No simulation, coordination, FOC, cost-recovery, model-space prior, or quality
covariance changes belong in this patch. Latent-quality policy propagation is
a separate substantive change with a separate specification and review.

## Reconstruction order and existing protections

Compute the unchanged noncentered market intercept before outside-share
reconstruction. For rotated sampling, reconstruct the old outside level/raw
effect, add the unchanged `X_s0 * beta_s0`, and call the existing
`bounded_inv_logit`. Retain the original share-mean expression using
`logit(s0_prob) + b_market`. The bounded-link correction must not disappear.

Evaluate existing priors on reconstructed old levels or standardized effects.
Keep prior constants and support normalizers. There is no prior on a new sum
coordinate independent of the old outside-share prior.

Keep the existing Logit Cournot kappa coordinates and implied-alpha priors.
They consume reconstructed `s0_prob` in exactly the same way as legacy.
CES retains its physical outside revenue-share calculation, outside price,
standardized coefficient, and elasticity transformation.

`gamma_market` in Logit and `b_market_rot` in CES are presently derived
diagnostics, not sampling coordinates. Preserve their historical convention:
multi-period markets use `b_market + mean(logit(s0_prob))`; one-period markets
use `b_market`; single-market normalization gives zero. New diagnostic sums
have distinct names.

## Validation contract

At matched economic states, native unconstrained targets must agree directly.
The standardized-shift map has determinant one, including when its shift
depends on hierarchy scales. Automatic gradients must obey the chain rule
through the full coordinate map. Dense covariances are permitted only in the
independent small reference, never in production.

Paired validation covers Logit/CES, Bertrand/Cournot, singleton and panel
dimensions, all pooling modes, centered/noncentered hierarchies, nonzero
outside covariates, sampled/collapsed quality, missing margins, supported rho
modes, and BF/CF/BS/CS ordered-k states. Compare economic outputs, the share
score, margin predictions, k probabilities, and normalized bridge targets.
Existing unsupported rho/likelihood combinations remain explicitly rejected.

Bridge sampling already evaluates Stan's normalized target. No separate
approximate bridge density is introduced. Statistical specification 9 denotes
the unchanged economic/statistical contract; coordinate provenance is recorded
separately. Compatibility with serialized legacy fits must be tested rather
than inferred from that version number.

## Isolation and evidence

Implementation uses `/tmp/market-s0-rotation-implementation`, an isolated
worktree whose Git metadata is writable. A fresh worker reviews the equations
and runs independent paired-target checks before HMC. Review and integrate
files into the active workspace only after those checks.

Freeze benchmark settings and starts before execution. Preserve unsuccessful
fits and warnings. Store heavy artifacts under an ignored Windows-drive
archive; only compact scripts, tables, and notes are versioned. Do not change
the default or add second-stage panel decomposition on the strength of code
inspection or a single attractive run.
