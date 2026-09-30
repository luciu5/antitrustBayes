# Sampler usability follow-up: presets and fit diagnostics

Status: proposed, not implemented. Recorded 2026-09-30 from the mergerBayes
banking ladder and calibrated information ladder, where getting fits past
strict HMC gates took repeated, hand-chosen sampler settings. No default,
prior, likelihood or kernel is changed by this note.

## Why

`bayes_est()` defaults to `adapt_delta = 0.95` and `max_treedepth = 13`
(`R/estimate.R`), and callers may override `adapt_delta`, `max_treedepth`
and `metric`. That is standard, but the settings the downstream work
actually needed were found fit by fit:

- Supreme Court Rung 3: per-model `adapt_delta` .999 or .9999, `dense_e`
  for Cournot, and more retained draws for stable bridge evidence.
- Supreme Court Rung 4: 16,000 retained draws and `warp3` bridging to bring
  the bridge log-evidence SD down, and `max_treedepth` 15.
- Information ladder: `adapt_delta` .999 and `max_treedepth` 14 or more for
  multi-market panels.

Most failures were a few events in thousands of draws, flagged by
deliberately strict gates. The fixes that mattered most were modeling
choices (a non-centered pooled-alpha kernel, observing all margins, a
tighter preset prior on the across-market alpha spread), not sampler tuning.
Sampler presets therefore cannot replace good parameterization or priors;
they only remove repeated trial and error once the model is sound.

## 1. Named sampler presets

Add an optional `sampler_preset` to `bayes_est()` that resolves to a
documented bundle of `adapt_delta`, `max_treedepth`, `metric`,
`iter`/`warmup`, and (through `bayes_compare()`) `bridge_repetitions` and
`bridge_method`.

- Presets are shorthand, not magic: explicit arguments always override them,
  and the resolved values are recorded in fit metadata so a fit states the
  settings it actually used.
- Candidate presets, to be validated before release: `"default"` (current
  behavior), `"careful"` (higher `adapt_delta` and tree depth for
  hierarchical or multi-market fits), and `"evidence"` (more retained draws
  and `warp3` when bridge evidence feeds model weights).
- Each preset needs a short validation record on the package's own fixtures
  (divergences, tree-depth hits, R-hat, ESS, bridge SD, run time). Settings
  chosen from one application's data are not a package default.

## 2. One clear diagnostic summary per fit

Today `show()` prints key parameter estimates with R-hat and bulk ESS, and
`summary()` prints the full CmdStan table. Neither gives a single verdict.
Add an exported `bayes_diagnostics(fit)` (and a short line in `show()`)
that reports, in one table:

- divergent transitions as a count and a fraction of post-warmup draws;
- maximum tree-depth hits against the configured limit;
- minimum per-chain E-BFMI;
- maximum R-hat and minimum bulk/tail ESS across sampled parameters;
- the sampler settings used;
- an overall status (pass / review) against documented default thresholds
  that callers can override.

When a check fails, print one line on what usually helps, for example:
divergences, raise `adapt_delta` and look for weakly identified scale
parameters; tree-depth hits, raise `max_treedepth`, which is usually cheap
when hits are rare; low ESS or high R-hat, run more iterations; unstable
bridge evidence, use more retained draws or `warp3`. Advice must never
suggest ignoring a failed check.

## Scope

Both items change user-facing API only. They must not alter the default
sampler settings of existing callers, the Stan kernels, priors or the
evidence calculation, and they need tests showing that explicit arguments
still win over a preset.
