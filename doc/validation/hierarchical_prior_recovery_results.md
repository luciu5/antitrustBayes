# Paired hierarchy SD prior recovery

Executed locally on 2026-10-11 with antitrustBayes 0.2.0.9003, commit
`74e4139`, and CmdStan 2.39.0. Four planned fits completed in 20 minutes
including reporting and an interruption. No production jobs were submitted.

## Assessment

Both prior families recover market alpha levels and give similar merger and
tariff predictions in these fixtures. Gamma does not uniformly improve
sampling: divergences decrease in the small case and increase slightly in
the moderate case; ESS, BFMI and runtime also have mixed results. Outside-share
heterogeneity remains weakly estimated, particularly in the small case.

These results support continuing the approved production experiment with its
specified priors and normal diagnostic reporting. They do not establish a
universally superior prior or motivate another pilot or prior adjustment.

## Fixed design

Each independently constructed Bertrand fixture has eight markets, four
periods and five single-product firms. Both fixtures use identical primitive
and measurement innovations. Within each fixture, the two prior fits use
identical observations, sampling seeds and other priors. An independent audit
verified that their Stan data differ only in the two hierarchy prior family
codes.

- Gamma: shape 2, rate 5 for log-alpha SD and logit-outside-share SD.
- Historical Half-Normal: scale .5 for both SDs.
- In these fits, the market-intercept SD prior is explicitly fixed at Half-Normal(.5).
- Small generating SDs: log-alpha .05, logit-outside-share .08.
- Moderate generating SDs: log-alpha .30, logit-outside-share .40.
- Generating market-intercept SD: .25 in both fixtures.
- Compositional share error: .025; log-margin error: .035; prices observed
  accurately. Measurement error follows equilibrium construction.
- Four chains, 4,000 iterations, 2,000 warmup, four cores, dense mass matrix,
  adapt_delta .995, depth 14. Sampling seed 518921 for every fit.

This tests the **package** Gamma defaults. The paper explicitly uses alpha-SD
Gamma rate 10; this exercise does not change or directly compare that setting.
The structural hierarchy is matched, while correlated compositional share
measurement remains an approximation to the independent Gaussian log-share
working likelihood. Two fixed fixtures do not establish repeated-sample
coverage, simulation-based calibration, or eight-model conduct recovery.

All 64 market-period equilibria passed fresh solves from perturbed prices.
Maximum absolute reproduction errors were 1.13e-8 for price, 1.11e-16 for
shares and 4.44e-16 for physical costs. All physical costs were positive.

## HMC diagnostics

Each fit retained 8,000 transitions. Summary extrema cover all 1,987 reported
quantities, including all active parameters. There were no depth-limit hits.
The nonzero divergence counts remain warnings.

| Heterogeneity | Prior | Sampling seconds | Divergences | Max R-hat | Min bulk ESS | Min tail ESS | Min BFMI |
|---|---|---:|---:|---:|---:|---:|---:|
| Small | Gamma | 87.4 | 4 | 1.0021 | 964 | 758 | .521 |
| Small | Half-Normal | 85.1 | 7 | 1.0021 | 1,310 | 1,163 | .615 |
| Moderate | Gamma | 200.2 | 2 | 1.0078 | 826 | 545 | .472 |
| Moderate | Half-Normal | 163.6 | 1 | 1.0046 | 819 | 618 | .494 |

## Parameter recovery

Errors below use posterior medians. Every fit covered all eight true market
alphas and all eight true outside shares in its 95% marginal intervals. This
is a descriptive count for these fixed fixtures, not a coverage estimate.

| Heterogeneity | Prior | Mean absolute relative alpha error, % | Mean absolute outside-share error, percentage points |
|---|---|---:|---:|
| Small | Gamma | .460 | 1.982 |
| Small | Half-Normal | .499 | 2.011 |
| Moderate | Gamma | .487 | 1.987 |
| Moderate | Half-Normal | .477 | 2.080 |

Hierarchy SD medians and 95% intervals:

| Heterogeneity | Parameter | Generating SD | Realized market sample SD | Gamma | Half-Normal |
|---|---|---:|---:|---|---|
| Small | Log-alpha SD | .050 | .040 | .049 [.028, .109] | .045 [.026, .093] |
| Small | Logit-outside-share SD | .080 | .053 | .184 [.038, .439] | .154 [.009, .438] |
| Moderate | Log-alpha SD | .300 | .241 | .264 [.165, .480] | .264 [.163, .484] |
| Moderate | Logit-outside-share SD | .400 | .266 | .374 [.168, .717] | .399 [.169, .752] |

The generating SDs and realized sample SDs are distinct targets; the finite
set of eight markets is less dispersed than the generating population here.
All hierarchy SD intervals include their generating truths. The unchanged
market-intercept hierarchy also has broad intervals covering .25 in all fits.

The [interval figure](hierarchical_prior_recovery_intervals.pdf) shows Gamma
in blue and Half-Normal in gray. Red dashed lines mark generating SDs; green
dotted lines mark realized sample SDs.

## Policy recovery

Targets are fixed market-1/period-1 HH merger with zero efficiencies and a
statutory 25% tariff on the other firms. Reported effects are weighted price
changes in percent. Each fit used 200 paired parameter/quality states for
both policies. An independent audit verified matching draw indices, plan
fingerprints and quality hashes, and recomputed the results from saved draws.
There were **zero structural-realization failures and zero policy failures**;
all 200 pairs were used in every fit. Prior pairs also used identical draw
indices. The posterior states themselves differ with the prior.

| Heterogeneity | Prior | Merger mean [95% interval] | Merger oracle | Tariff mean [95% interval] | Tariff oracle |
|---|---|---|---:|---|---:|
| Small | Gamma | 5.400 [5.001, 5.984] | 5.689 | 6.085 [5.951, 6.260] | 6.170 |
| Small | Half-Normal | 5.416 [4.912, 5.854] | 5.689 | 6.090 [5.929, 6.235] | 6.170 |
| Moderate | Gamma | 4.628 [4.039, 5.197] | 4.767 | 6.085 [5.891, 6.260] | 6.119 |
| Moderate | Half-Normal | 4.587 [4.176, 5.175] | 4.767 | 6.069 [5.928, 6.262] | 6.119 |

All eight intervals include their oracle effects. Maximum absolute mean
errors were about .288 percentage points for merger and .085 for tariff.
The largest prior-induced shifts were about .042 points for merger and .016
points for tariff. These are predictions conditional on the true Bertrand
model, not eight-model BMA predictions.

## Evidence

Three normal bridge repetitions used the same fitted Bertrand model under
each prior. Changes in prior density are expected to change its evidence;
this is not a conduct-selection test.

| Heterogeneity | Prior | Log marginal likelihood | Bridge repetition SD |
|---|---|---:|---:|
| Small | Gamma | 475.336 | .009 |
| Small | Half-Normal | 475.795 | .005 |
| Moderate | Gamma | 460.721 | .071 |
| Moderate | Half-Normal | 460.431 | .084 |

Gamma minus Half-Normal log evidence is -.459 in the small case and +.291
in the moderate case. Both differences are modest. These bridge precision
summaries do not remove the divergence caveat.

## Execution record

The initial reporting step encountered a `draws_matrix` indexing error after
the first fit had been saved. Recovery extraction now converts to an ordinary
matrix. Explicit resume reused that posterior, checked identical fit arguments
and retained the original two-hour deadline; no HMC fit was repeated.

After all four tasks and the completed manifest were written, the resumed
process exited nonzero on a trailing parse error: the script had been edited
on disk while R was still reading it. Fresh parsing and independent artifact
checks passed afterward. The scientific outputs were retained and verified;
the process exit is not recorded as a clean end-to-end script run.

The process also reported more than 50 accumulated R warnings. An independent
replay of the small-Gamma policy calculation with its exact saved plan
reproduced every policy draw and summary. It recorded 400 benign API notices:
200 for unchanged tariff ownership and 200 for Logit ignoring BLP share-
inversion inputs. No solver or cost warning occurred in that replay. These
notices are separate from the retained HMC divergences and rejected early
warmup proposals. Other fits' notice counts were not reconstructed.

Raw posteriors, CSVs and logs remain outside Git at
`/tmp/hierarchical-prior-recovery-pJYUpF`. The completed manifest records start
2026-10-11 22:34:31 UTC and finish 22:54:32 UTC, before its original deadline.
Fixture SHA-256:
`48703dc430a4a38103b130c9c389dd3d47071582bfe37983224a60d2b67060ce`.

Reproduce with [hierarchical_prior_recovery.R](hierarchical_prior_recovery.R)
and the commands in [README.md](README.md#paired-hierarchy-sd-prior-validation).
Package defaults, the paper production scope and production configuration
were unchanged by this validation.
