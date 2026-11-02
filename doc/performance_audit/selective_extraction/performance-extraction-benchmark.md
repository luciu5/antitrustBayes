# Selective posterior extraction benchmark

Date: 2026-09-16. This is a read-only benchmark in the isolated review
clone. The archived fits are unchanged. Policy cells use 200 fixed posterior
draw positions from each 2,000-draw fit, three policies, and three warmed
sequential repeats. Full and selected runs use the same seeds and controls.

## Archived policy comparison

All four fixtures have 9/9 complete rows, 200/200 valid draws per row, and a
zero solver-failure rate in both modes. Across all 36 paired rows, draw
indices and failure vectors are identical. The maximum absolute difference
between full and selected results is zero for prices, shares, compensating
variation, producer surplus, and total welfare.

| fixture | full/selected extracted bytes | full/selected fields | extraction seconds |
|---|---:|---:|---:|
| ordinary Logit | 3,195,776 / 33,736 | 37 / 2 | 0.009 / 0.003 |
| exact BF | 5,273,120 / 130,536 | 54 / 3 | 0.020 / 0.004 |
| exact BS | 5,209,120 / 114,536 | 54 / 3 | 0.018 / 0.003 |
| ordinary CES | 3,162,200 / 33,736 | 35 / 2 | 0.013 / 0.003 |

Policy solve wall time has no consistent end-to-end improvement because the
solver dominates. Mean selected/full ratios by policy are:

| fixture | merger | cost +10% | merger + cost +10% |
|---|---:|---:|---:|
| ordinary Logit | 0.988 | 1.056 | 1.023 |
| exact BF | 0.982 | 1.011 | 0.994 |
| exact BS | 1.050 | 1.032 | 0.931 |
| ordinary CES | 0.992 | 0.969 | 1.059 |

## Fixed_param extraction stress

These cells are extraction and serialization stress checks, not posterior
inference and not HMC. Each uses the same `logit_supply` kernel, valid
collapsed-quality data, 200 Fixed_param rows, and `G=min(N,1000)` persistent
groups in separate full and selected processes. Canonical `a_market` and
`s0_prob` hashes match across modes and all retained values are finite.

| N | G | T | full/selected extract seconds | full/selected retained bytes | full/selected peak RSS KB |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 1,000 | 1 | 0.039 / 0.020 | 11,285,280 / 4,936 | 566,092 / 552,652 |
| 10,000 | 1,000 | 10 | 0.834 / 0.025 | 112,186,080 / 19,336 | 881,796 / 772,632 |

The Fixed_param setup itself took approximately 0.15 s (N=1,000) and 1.32 s
(N=10,000) in both modes; it dominates the small extraction timing at
N=1,000. Both stress processes stayed below the 1.5 GB review ceiling.

Scripts and compact result files are in this directory. The individual RDS
and `/usr/bin/time -v` logs are retained only to make the measurements
auditable; they are not production artifacts.
