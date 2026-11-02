# CF/CS rate-domain investigation

The explicit real-rate implementation has passed the production replay,
package regression tests, and independent integration review. Changes are in
both working trees. The agent's attempts to apply the coordination patch were
blocked by its read-only filesystem; the user subsequently applied it from
their terminal. All seven sibling files match the tested source snapshot,
and reverse-patch and whitespace checks pass. Neither repository has been
committed or pushed in this development step.

## Analytically established

The input Logit Cournot CF/CS interior FOCs admit the firm-level scalar-W
solution derived in `logit_rate_domain_math.md`. Its implied rates are finite
real values under finite, representable interior primitives; positivity is an
additional action-domain restriction. The extension retains the same FOCs,
cost/value primitives, roles, and within-firm product composition. It does not
derive a rate-floor equilibrium.

## Independently numerically verified

`logit_real_rate_reference.py` executed 96 quantity-profit checks for negative
rates, including multiproduct firms, one/several/all strategic leaders/core
firms, and a 10% value shock. Maximum errors were:

| Check | Maximum error |
|---|---:|
| Demand shares versus equilibrium quantities | 8.33e-17 |
| Finite-difference own/reduced profit derivative | 1.11e-10 |
| Direct one-coordinate best response versus equilibrium quantity | 3.53e-9 |

The existing observed-anchor audit retains 400 draws: 36 CF and five CS
positive-domain failures all coincide with a nonpositive rate in the unique
real interior solution. Previously successful solutions agree closely with
that independent oracle. Original tables remain unchanged.

The new latent-anchor audit likewise retains 400 draws, using the saved
permuted posterior indices, conditional quality draws, and exact-k states:

| Model | Requested | Saved failures | Nonpositive pre | Nonpositive post | Max original-FOC residual |
|---|---:|---:|---:|---:|---:|
| CF | 200 | 26 | 0 | 26 | 3.29e-12 |
| CS | 200 | 1 | 0 | 1 | 1.70e-12 |

All 27 latent failures occur after the merger, and all have k=2. The minimum
post rate belongs to the acquirer. The latent oracle agrees with saved
successful pre/post prices within 7.45e-13 and shares within 9.57e-14. Quality
and observed-FOC procurement values are held fixed; no model was refitted.

## Production validation

The production replay reuses the archived 200-draw plan for each CF/CS fit
under both observed and latent-quality anchoring. No fit, quality vector, k
state, or procurement value was replaced. The original inside quantity budget
is retained exactly. The no-cost-change merger is the archived acquirer
298816 / target 889018 scenario.

| Model | Anchor | Requested | Original positive-domain failures | Real-domain failures | Max original-FOC residual | Simulation seconds |
|---|---|---:|---:|---:|---:|---:|
| CF | Observed | 200 | 36 | 0 | 2.81e-13 | 2.993 |
| CF | Latent quality | 200 | 26 | 0 | 2.85e-13 | 3.773 |
| CS | Observed | 200 | 5 | 0 | 5.16e-14 | 4.490 |
| CS | Latent quality | 200 | 1 | 0 | 3.91e-14 | 5.405 |

All 68 original failures have a negative post-merger rate in the recovered
solution. There are no nonpositive pre-merger rates in these 800 draws. Across
all requested draws, maximum differences against the independent oracle are
3.27e-13 for prices, 6.26e-14 for shares, and 9.60e-14 for procurement values.
Across previously successful draws, maximum differences are 1.33e-12 for
prices, 4.19e-13 for shares, and 1.23e-12 for inherited welfare outputs.
This is a deterministic replay, not additional sampler evidence. The reported
times include simulation startup and assembly, but exclude loading archived
fits and writing replay artifacts. Old wall times were not remeasured, so
these values do not establish a speed ratio.

The real CF/CS routes reevaluate demand and markups at the returned prices,
requiring the absolute price FOC to be at most 2e-10. They also retain a raw
profit-FOC diagnostic gate of 2e-7. The new gate is additional coverage, not a
replacement by a weaker normalized criterion. Tests cover multiproduct
aggregation, ownership/cost and combined shocks, all leaders, subsets removing
all leaders, analytic/implicit CS agreement, exact zero rates, dollar markdowns,
welfare, typed positive-domain failures, and positive-route equivalence.

The antitrustBayes focused suite passed 57 assertions; the final full installed
suite passed 1,966 assertions with no failures/errors, 13 opt-in skips, and
13 existing expected warnings. The final installed coordination suite passed
572 assertions in 71 test blocks, with no failures/errors/skips and 24 existing
expected warnings. The exact-zero Stackelberg fixture has utility u(p)=p,
beta=1, and value=2: the production root returns price=0, share=1/2, and
dollar markdown=2. Undefined normalized residuals are NA with a product flag;
the raw FOC remains finite and both acceptance gates are unchanged.

Final test objects, counts, raw logs, and the last all-800 replay are retained
under the local run's `final-validation-01/`. The fresh independent reviewer
approved the final source and validation. Its earlier installed coordination
run independently passed 563 assertions; the final run adds the exact-zero
Stackelberg checks. Source expectation counts are not executed assertion counts;
reported totals above are computed from saved test result objects.

The Gaussian/Student-t, quality scope, priors, conduct, and Stan files are
unchanged. The installed native Stan library retains SHA-256
`50b8a7ee6b3d32a82df3c5368962f219f614247d5606e4c6900bf1011adbad7a`.
The default domain remains `positive`; the new domain requires an explicit
`price_domain = "real"` request and coordination >= 0.0.0.9004.

## Retained failed validation run

The first production replay stopped at its unchanged absolute welfare gate
(1e-7), despite matching CF prices/shares closely. The replay script had
incorrectly reconstructed `insideSize` from the realized latent inside share.
Latent Logit replacement changes shares while preserving the original quantity
budget; the erroneous replay input therefore rescaled welfare by up to 0.0254.
The script now reuses the saved `insideSize` and checks exact budget parity.
No production economics, residual tolerance, or welfare tolerance was changed.
The failed `production-replay-01` artifacts and its log remain preserved beside
the successful `production-replay-02` run.

## Reproducible artifacts

`logit_real_rate_pnb_reference.R` constructs the independent 800-draw oracle.
`logit_real_rate_pnb_replay.R` verifies the production path against that oracle
and archived successes. Both accept the local market-quality run root and a
new output directory; existing results are never overwritten. Raw fits, plans,
simulation objects, and logs remain in ignored local directories. The compact
four-row production table is `logit_rate_domain_replay_validation.csv`; the
27-row latent failure classification is `logit_rate_domain_latent_failures.csv`.

The complete seven-file coordination patch is retained at
`doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/logit-real-rate-01/coordination-real-rate-final.patch`.
Its SHA-256 is
`b23207b95dc410c19399232b363968c2bfe80977fec09183d75d6810d97bae84`;
its tested base is coordination main
`a7a928a981c42f4f57751daca858bb0d4081cd46`. `git apply --check` succeeds on
the original sibling checkout. The user applied the patch, and all seven
files were verified against the final tested snapshot. It increments
coordination to 0.0.0.9004. Install the updated packages before using the
public real-rate API with the normal R package library.

## Statistical and economic limits

No posterior density, prior, likelihood, bridge target, or HMC settings were
changed. This deterministic domain extension does not resolve the
separately recorded PNB sampler gate failure or authorize withheld BMA output.

Real rates require an application interpretation such as net depositor fees.
If nominal deposit rates are constrained to be nonnegative, a constrained game
requires its own derivation. Clipping the interior solution at zero is not a
valid substitute. CES, output, Bertrand, and ordinary simultaneous routes are
outside this initial implementation.
