# Additional native Helmert review

This directory contains a fresh, narrow review of the candidate sampled-quality
Helmert recurrence against the preserved baseline native kernels. It is
separate from `native_profile/run_helmert_panel_gate.R` so the original 18-case
profile remains unchanged.

The runner uses one candidate-generated unconstrained coordinate vector per
case and replays that vector in the baseline process. It records normalized and
raw targets, full normalized and raw gradients, representative finite
differences, constrained transformed parameters, complete fixed-draw generated
quantities, parameter schemas, finite-value patterns, and Helmert centering and
boundary values. The added fixtures are:

* Logit sampled quality, share-only, panel data;
* CES sampled quality, margin-only, panel data;
* Logit smooth cutoff with HMT enabled;
* Logit bargaining;
* CES second-score auction;
* Logit `G = 1000` with one product per group;
* CES `G = 1000` with the same 1000 groups repeated across two markets.

Run the candidate first, then the baseline with the same raw-coordinate file:

```sh
HELMERT_ROLE=candidate Rscript doc/performance_audit/recovery/native_review/run_extra_helmert_gate.R
HELMERT_ROLE=old Rscript doc/performance_audit/recovery/native_review/run_extra_helmert_gate.R
Rscript doc/performance_audit/recovery/native_review/compare_extra_helmert_gate.R
```

The runner resumes completed per-case checkpoints by default. Set
`HELMERT_FORCE=true` after changing either installed library or the fixture
source; a failed or interrupted case is then recomputed explicitly.

The installed package libraries are selected by the caller through `R_LIBS_USER`
or an equivalent library path. The runner does not install packages, modify
Stan source, or modify the active package source. Large RDS snapshots stay in
`.artifacts/`; compact summaries and failure ledgers are written beside this
README. The cache test exercises restart and replay guards with fake callbacks;
it does not establish a native equivalence result.
