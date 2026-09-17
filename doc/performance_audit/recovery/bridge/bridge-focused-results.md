# Focused bridge recovery results

The focused recovery file was run from the durable R-only candidate package on
2026-09-17. The package R code came from commit `1bf79e33c045f46616d590657780915b4dabef22`; focused tests used `50ad896c6f481d69a21a1fa0bf054811c271a27`;
the installed R-only package reused the unchanged baseline native DLL and did
not compile Stan or C++.

Command:

```text
R_LIBS_USER=Rlib/r-only:Rlib/baseline:Rlib:Rlib/candidate:<defaultuser> \
Rscript --vanilla -e 'library(antitrustBayes); library(testthat); test_env <-
new.env(parent = asNamespace("antitrustBayes")); test_dir(
"tests/testthat", filter = "bridge-performance-recovery", env = test_env,
load_helpers = FALSE, reporter = "summary")'
```

Result: exit status 0, no failed tests. The suite covered RNG presence and
restoration, seed and core type validation, deterministic L'Ecuyer streams,
serial worker errors, parallel worker errors with `mc.set.seed = FALSE`, the
serial `bridge_seed = NULL` path, model-universe control forwarding, and public
first/middle/last context restoration failures with two repetitions each.

The six displayed worker warnings are intentional synthetic failures: three
from serial dispatch and three from parallel dispatch after the first matching
warning was captured by `expect_warning`. The complete warning stream is in
`test-bridge-performance-recovery.log`; exercised failure rows are in
`bridge-focused-failures.csv`.

The installed function contract reports `bridge_maxiter = 5000L`,
`bridge_repetitions = 3L`, and universe formals including `bridge_cores` and
`bridge_seed`. The candidate and R-only installed Stan source hashes match;
the baseline and R-only native DLL SHA-256 is
`31d1beca03ac862872b1f96c9943c2b082c1d93b105cb64f4854aef51999afdc`.

No archived-fit bridge timing was run.
