# Recovered serial vectorization gate

This is a computation-only candidate for the ordinary Logit and CES kernels.
The supplied patch replaces sums of normalized scalar Normal/Student-t
densities with their normalized vectorized equivalents. Predictors, scales,
degrees of freedom, conditional-rho formulas and missing-margin indices must
remain identical. The collapsed joint-Gaussian block, economic equations,
parameters, priors, exact-k kernels and generated quantities are unchanged.

Use the verified Helmert package in `recovery/Rlib/candidate` as the old
kernel. Reconstruct the vectorized source in its own ignored durable worktree
from `bayes-vectorization-unapproved-candidate.patch`; never edit the actual
main kernels at this stage. Match the compiler, flags and dependency manifests
of the recovered native builds. Compile one package at a time with `-j1`.
Keep build manifests, logs and DLL checksums in the ignored recovery directory.

Before timing or HMC, replay the eighteen fixed panel fixtures at identical
unconstrained points in separate processes. Compare adjusted and unadjusted
targets, complete gradients, raw schemas, constrained parameters, generated
quantities and representative finite differences. Include full large-N
gradient vectors, rather than only their norms. Exercise sampled Normal and
Student-t shares, rho on/off, missing and zero margins, smooth fringe,
collapsed shares with independent Student-t margins and the untouched
collapsed joint-Gaussian path. Retain every failure and warning.

For the warmed serial benchmark, freeze probes at N=100, 1,000 and 10,000,
with one-product and multiproduct/repeated-group data. Evaluate old/new probes
sequentially without other CPU-intensive tasks. Save every repetition; report
median target, gradient and GQ time separately, including setup and memory
where measured. A kernel speed improvement is not an HMC ESS/second result.

The original paired HMC arrays and initialization files were lost during the
environment reset. Any replacement is a new, explicitly labeled recovery
experiment. Freeze data, priors, raw initial states, seeds, source/DLL hashes
and all sampler settings before running either kernel. The preserved plan was
N=100, G=20, cyclic group IDs, prices evenly spaced from 1 to 2, quantities 1,
margins .2; Logit sampled Student-t shares/margins and CES collapsed Gaussian
shares/Student-t margins; rho and all optional pooling/cutoff/HMT flags off,
legacy coordinates. Seeds 61001/61002, four chains, 1,000 iterations/500
warmup, adapt_delta .95, depth 13, diag_e, cores 1. This artificial equal-share
fixture is a stress case, not a recovery DGP. Do not change its settings after
seeing diagnostics. Checkpoint each fit and capture warnings, sampler clocks,
divergences, depth hits, BFMI, structural bulk/tail ESS, ESS/second and economic
posterior summaries. Poor diagnostics remain evidence against acceptance.

No candidate enters main because a compiler succeeded or because a single
fit had no divergences. Root reviews the equations and numerical comparisons;
fresh independent review precedes integration. Reject or hold changes whose
small measured gain does not justify their complexity. Do not push remotely.
