# Coordinated release gate, 2026-11-01

**Passed: 20/20 frozen inference gates and 32/32 independent policy FOCs.**
Exactly two fresh CmdStanR fits completed; no failed fit was removed, refitted,
or replaced. Package `0.2.0.9009` and model specification 14 are unchanged.
The source baseline was `fb1f026c0e301f17343ebd611afe657cbaafc5f2`.

The [sealed runner](coordinated_release_gate.R) uses the current
`relative_validation_fixtures()` Logit C M4T4 (112 observations) and
`relative_noisy_ces_panel_arguments()` applied to CES C M4T4 (64 observations).
Both explicitly use relative shares: four chains, 1,000 total iterations,
500 warmup, `adapt_delta=.99`, depth 13, diagonal metric, 18-digit CSV output;
seeds are 726041 and 726042 respectively. CmdStanR 0.9.0 and CmdStan 2.39.0
completed every chain, retaining 2,000 draws per fit. Native metadata and draw
shape independently confirm the frozen settings. Including compilation,
Logit took 736.058 seconds and CES 192.047 seconds.

| Frozen criterion | Logit | Finite-noise CES | Requirement |
| --- | ---: | ---: | ---: |
| Market structural truths in 95% intervals | 4/4 | 4/4 | 4/4 |
| Maximum structural median relative error | 0.0829% | 1.2640% | <=8% |
| t1 outside-share truths in 95% intervals | 4/4 | 4/4 | >=3/4 |
| t1 outside-share median MAE | .000287 | .006659 | <=.06 |
| Maximum active-parameter R-hat | 1.010110 | 1.010796 | <=1.05 |
| Minimum active-parameter bulk ESS | 603.08 | 584.16 | >=100 |
| Divergences | 0 | 0 | 0 |
| Share-scale median / actual error RMS | 1.160 | 1.415 | [.5,2] |
| Margin-scale median / actual error RMS | 1.203 | 1.130 | [.5,2] |

Diagnostics cover all 26 active sampled Logit parameters and all 16 CES
parameters, excluding transformed constants. The four-market outside-share
gate fixes t1 before sampling because Logit outside shares are unpooled across
years; [all 16 cells per family are reported](coordinated_release_gate/recovery.csv).
CES margin-error RMS is .0141818 and its posterior scale median .0160192;
no historical Logit noise band is applied to CES. Share RMS removes the common
log-share direction within each market-period cell and uses the contrast
count as its denominator. [Noise evidence](coordinated_release_gate/noise.csv)
and [every gate](coordinated_release_gate/gates.csv) retain the frozen decisions.

On the final reviewed core in `/tmp/coordinated-release-fixed-lib`, the true
state and retained draw 83 (public plan seed 726051, fixed before fitting)
pass baseline and first-two-firm merger FOCs in all four markets, for both
families. The [32 residuals](coordinated_release_gate/policy_focs.csv) have maximum
absolute error `9.11085e-9`, below `1e-6`; no policy realization failed. The
oracle reconstructs demand and ownership equations directly, without
`calcMC()`, `calcMargins()`, or `calcShares()` as an oracle. CES preserves the
legacy Cournot convention of fixed calibrated outside-revenue alpha across
policies. Independent reviewer `review_economics` audited the unchanged
runner and found no substantive issue or circular margins oracle; the reviewed
runner SHA-256 is
`05184a3561a0f54b38f478d1b99bb1405c9a10d9aa0d8d74e23ead2eab869870`.

Ordinary source verification on a stable corrected-stack snapshot: **7,428
passed, zero failures/errors, 16 skips, 22 retained warnings**. Focused sampler
checks: 124 passed, no failures/errors/skips/warnings. Full-vignette package
build and `R CMD check --no-manual`: **zero errors, warnings, and NOTEs**;
`utils::tail()` qualification removes the former code-check NOTE. After the
final core's additional ownership/numerical fixes, targeted policy regressions
passed 803 expectations, with zero failures/errors, one opt-in skip and 12
retained warnings, alongside the 32 actual FOCs above.

[Exact source hashes](coordinated_release_gate/source_hashes.csv),
[fit/data/kernel/settings identities](coordinated_release_gate/fit_identity.csv),
[all active-parameter diagnostics](coordinated_release_gate/active_parameter_diagnostics.csv),
[runtime artifact hashes](coordinated_release_gate/runtime_hashes.csv), and
[durable artifact seals](coordinated_release_gate/artifact_hashes.csv) bind the
execution. Inference used an unchanged isolated library snapshot; final policy
checks used the final reviewed core, with separate recorded runtime hashes.

Raw fits, every chain CSV, native caches, manifest, plans, warnings, all execution
logs and the compact exporter remain in the ignored durable directory
`/mnt/d/Projects/antitrustBayes/.scratch/coordinated-release-gate`.
Preparation and packaging failures are retained there: a shared dependency
reinstall race; an initial no-vignette packaging attempt whose tutorial stopped
before HMC at unavailable default-cache creation; and a shared-library check
that fell back to old `trade`. The first compact exporter also misinterpreted
CmdStan's per-process `num_chains`; its reporting assertion was corrected to
use the fit's four-chain count. None changed a fit, source seal, scientific
setting, or threshold. The final package check and export succeeded.

This is bounded known-truth recovery for two current relative-share fixtures.
It establishes no repeated-sample coverage, general conduct identification,
bridge evidence result, or completion of the archived 124-task enterprise grid.
