# Executed package and numerical gates

The deterministic implementation and native-specification-8 Bayesian code
pass their package gates. No multi-seed statistical recovery or PNB result
is established by these checks.

| package | source check | test passes | expected warnings | failures | skips |
|---|---|---:|---:|---:|---:|
| coordination 0.0.0.9003 | Status: OK | 383 | 24 | 0 | 0 |
| antitrustBayes 0.2.0, native spec 8 | Status: OK | 971 | 12 | 0 | 0 |

The Bayes check was a normal fresh installation and full source check with
ANTITRUSTBAYES_RUN_STAN_TESTS=true, NOT_CRAN=true and MAKEFLAGS=-j1. Examples,
compiled code, namespace, documentation and both vignette rebuilds pass.
The original failed test-path check and unsuccessful installation-skipping
shortcuts remain local evidence; none constitutes the passing release check.

The reviewed coordinate source is worker commit
66cf923cf288e5d6adb7cea7304afa4d994c582f. The checked source archive SHA256 is
803d4e1decfb3e08b3d22edd29fefa54877974074de4605638b9b176b2635526.
The lead independently verifies all19 changed source paths byte-for-byte
against that archive. The native runtime is strict specification8.

[quality_coordinates.md](quality_coordinates.md) establishes the unchanged
physical-quality prior and normalized evidence. Compiled mapped-target,
prior, whole-unit likelihood, finite-difference gradient, deterministic
reconstruction, combined-shock and nesting gates pass; unchanged tolerances
and actual maxima are in [native8_gate_summary.csv](native8_gate_summary.csv).
The standalone Stackelberg FOC audit and derivative tables remain authoritative
for economics. No existing coordinated PriceLeadership class was repurposed.

The specification8 pilot completed at 2026-09-14 13:36:52 UTC. All twelve
first-seed, single-product fits are complete: one passes the sampling gate and
eleven fail. The four economics cases each have three candidate structures.
All 48 chains have 1,000 finite saved draws per chain. The
Logit Bertrand Stackelberg case has worst classic Rhat 1.2435 and minimum
effective sample size 19.26; CES Bertrand Stackelberg has Rhat 1.0466 and ESS
16.39. These failed fits do not establish recovery or valid policy inference.
Four matching native7/8 cases have identical data, priors, sampler settings and
seeds; coordinate equivalence has not established successful mixing.
The remaining 132 cases of the specification8 144-case recovery protocol
have not started. The separate older specification7 longer pilot still had
ten completed failures and two running cases at 2026-09-14 13:50 UTC; its
counts are not combined with the specification8 protocol.

The larger-data helper preserves the original exact-share DGP but has not
been fitted. Increasing its periods does not automatically resolve the local
location/scale geometry. The larger-data screen remains held; see
[exact_share_scale.md](exact_share_scale.md). Fits, binaries, full logs, raw
research history and large delivery archives remain local and are excluded
from the source branch intended for GitHub.
