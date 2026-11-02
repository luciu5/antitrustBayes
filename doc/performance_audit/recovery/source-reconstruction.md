# Isolated source reconstruction

The active Bayes worktree stayed unchanged. Temporary copies are under the
ignored statistical-toolkit recovery tree.

- Bayes baseline: `e67cf31` plus the reviewed current `R/bma.R`, `R/params.R`,
  and baseline consumer test/doc overlays.
- Bayes candidate: baseline copy plus
  `bayes-helmert-candidate.patch`, SHA-256
  `26b5e55a8fdcb0a31f1d5537361ee3f586fff5934f2275bca85ed9efaff15783`.
- antitrust dependency: refactor HEAD `684c009`, with
  `existing-antitrust-promotion-base.patch`, SHA-256
  `cf9a11577b4de7a1897977e5a00f85f37750a5e186750ad23ab554a6a3875f3b`,
  applied only in the temporary copy; installed as 0.99.40.
- coordination dependency: source HEAD `a7a928a`, installed as 0.0.0.9003.
- trade source was inspected at refactor HEAD `6609e28`; it was not needed for
  the native Bayes compile and its working tree was not modified.

The actual antitrust and trade branches remain untouched. The source and
dependency hashes used by the restart guard are in the ignored recovery tree;
the compact build and environment checkpoints are beside this file.
