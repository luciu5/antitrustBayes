# Optional bridge candidate: archived benchmark

## Executed evidence

Actual main is committed and its locally recorded origin/main matches ea29e96.
No bridge-dispatch implementation was added to main. The candidate remains in
an ignored isolated worktree; source, plans and compact evidence are retained.

Eight benchmark jobs (each evaluating two separate public comparisons,
two seeds, four execution modes) completed. Each fit retained 5,000
maximum bridge iterations and three repetitions: 48 bridge estimates in total.
There were no failed repetitions or bridge warnings in these calls. Every
paired model/repetition log marginal likelihood and iteration count was exactly
equal. Input fit hashes were unchanged. The compiled native DLL for every
mode was 417e584bd6fb96f7b3aa65860c953fe6e09c8fdf715a8daeca49a50ccc0a1aaa.

The archived HMC evidence is existing, not new: both BF and BS fits were marked
usable, with zero divergences and zero treedepth hits. Their maximum R-hats
were 1.00635 and 1.00498; minimum BFMIs 0.93593 and 0.77199. Minimum bulk ESS
was 981 and 1,032, and minimum tail ESS 715 and 778, respectively. These
diagnostics do not establish that the later bridge estimates themselves have
small Monte Carlo uncertainty or support a cross-model scientific conclusion.

| Seed | Serial seconds | Two-worker seconds | Speed ratio | PSS ratio |
| --- | ---: | ---: | ---: | ---: |
| 707101 | 9.332 | 8.211 | 1.137 | 1.784 |
| 707102 | 9.751 | 8.233 | 1.184 | 1.783 |

Serial public timing covers two separate full comparisons, including validated
context restoration and result assembly. Parallel execution runs each fit's
three repetitions concurrently with two workers. Process-tree PSS includes
startup and all children, sampled every 0.25 seconds. Peak PSS rose from
774–776 MiB to 1,380–1,384 MiB. This is not a measured large-market speedup
or a decomposition of bridge kernel time alone.

## Defects and retained failures

The original combined BF/BS comparison failed the same-observations contract:
the archived fits have different DGP datasets. The validation was correct.
Original plan, runner and failure log remain saved. The corrected performance
plan treats them separately; it does not compare their scientific posterior
model probabilities. A missing existing R library path also caused an initial
metadata-only load failure; fixing the path did not change any dependency.

A synthetic warning test exposed lost caller-visible fork warnings: serial
dispatch delivered two warnings, fork dispatch zero, with no returned warning
records. The candidate now captures worker warnings with model/repetition,
returns them, and re-emits them once in the parent. The independent new test
has 13 passing assertions, including warnings followed by errors, stream/result
equality and caller-visible parity. Its full output is retained in the ignored
recovery/logs/bridge-worker-warnings.log. The original candidate source is
preserved under bridge-before-warning/. Earlier focused recovery test results
predate this repair; they are explicitly not presented as current full-suite
evidence. Current-source R-only installation and focused warning, recovery,
migration and context-restoration tests now pass with the unchanged Helmert
DLL and model hashes. The private help update also passes parsing and installed
documentation checks. Initial missing-dependency and interrupted load checks
remain saved. The full installed suite was not rerun for this held candidate;
see candidate_package_verification.md for the exact verification scope.

Fresh Luna read-only review found no remaining substantive defect in the
repair, ordering, NULL-slot preservation, error records, RNG handling, source
migration safeguards or model-universe forwarding. No economic equation,
target helper, prior, default repetition count or iteration safeguard changed.

## Recommendation

Hold this candidate out of main. A 14–18% improvement on small fixtures with
approximately 78% more peak process-tree PSS does not yet justify the added
orchestration complexity, especially on this unstable 8 GiB WSL environment.
Keep serial production bridge execution. A larger valid saved-fit benchmark
may justify reconsideration, without refitting to manufacture a favorable
case or weakening evidence diagnostics. The isolated implementation and
warnings repair are preserved so this decision loses no research work.
