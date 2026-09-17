# Dependency install checkpoint

The isolated dependency stage completed on 2026-09-17. Source copies are under
the ignored `doc/statistical_toolkit_runs/20260916/performance-audit/recovery/worktrees/`
directory and the library is under its `Rlib/` sibling.

- antitrust 0.99.40: antitrust refactor HEAD `684c009`, with
  `existing-antitrust-promotion-base.patch` applied after normalizing only the
  temporary copy's CRLF files.
- coordination 0.0.0.9003: coordination source HEAD `a7a928a`.
- installed dependency library: `Rlib/`.
- install logs: `install-antitrust.log` and `install-coordination.log`.

The first antitrust attempt used an isolated-only `R_LIBS_USER` and correctly
failed because it hid the already installed user-library dependencies. The
successful retry prepended the isolated library while retaining the user
library, with `MAKEFLAGS=-j1`.

The Bayes build harness is
`doc/performance_audit/recovery/run_native_install.sh`. It installs baseline
and Helmert candidate into separate `Rlib/baseline` and `Rlib/candidate`
libraries using the same flags, records compiler hashes and Stan source hashes,
and writes compact logs/manifests beside this file.
