# Working in antitrustBayes

## Start with the current task

- Routine implementation may be performed directly on `main`, the maintained default branch. The `refactor` migration is integrated; do not recreate that branch. Start any necessary isolated task branch from `main`.

- Read [the documentation index](doc/README.md), then only the topic documents
  needed for the requested change. It owns the current operational agenda;
  historical plans and session checkpoints are evidence, not instructions to resume.
- Check `git status --short` before editing and preserve unrelated work.
- Use `DESCRIPTION`, current R/Stan source, help pages, and tests to establish
  implemented behavior. Use dated result reports for executed evidence. A proposed
  interface or a passing algebra check does not establish statistical recovery.
- Update the index when a task changes a recorded status or validation route;
  link detailed evidence rather than copying it into another status document.

## Repository and package boundaries

- `R/` owns estimation, posterior reconstruction, evidence, and policy draw plans;
  `inst/stan/` owns the four Stan kernels; `man/` contains generated help.
- `tests/testthat/` holds routine regressions. `doc/` contains mathematical
  contracts, research runners, compact evidence, and historical checkpoints.
- `antitrust` owns structural models and common counterfactual machinery;
  `coordination` supplies fringe/Stackelberg models. Do not duplicate their
  equilibrium solvers here or add a dependency on downstream applications.
- Read dependency minima in `DESCRIPTION` and the development pin in
  `.github/workflows/qa.yml`; an installed package can differ from this checkout.
  Select a compatible R library explicitly before running integration checks.

## Choose proportionate verification

- Run commands from the repository root. Use `testthat::test_local(filter = ...)`
  for development-source checks of the behavior changed; see the index for
  commands and dependency preflight. Run the ordinary suite for broad changes.
- Keep opt-in Stan sampling, fixture checks, known-truth recovery, and large
  research benchmarks separate. Do not run HMC for documentation-only edits.
- The enterprise 124-task runner is archived RStan infrastructure; use the
  current CmdStanR validation entrypoint, and do not treat fixture checks as HMC.
- New tests should detect a named economic, numerical, API, or downstream failure.
  Prefer independent identities and small fixtures over implementation mirrors.

## Preserve scientific and execution evidence

- Keep demand/supply equations, likelihoods, priors, pooling, role support,
  normalization, and solver/sampler settings fixed during performance work
  unless the task explicitly changes them. Label separate experiments.
- Preserve failed/incomplete fits and unavailable evidence. Do not normalize
  model probabilities over successful positive-prior candidates only.
- Reuse compatible saved fits and fixed posterior draw plans when valid;
  verify source, prior, data, and draw provenance before avoiding a refit.
- Save long-run scripts, manifests, and completed-fit checkpoints in durable
  ignored workspace locations. `/tmp` session paths in old reports may be gone.
  Track compact evidence and reports, not raw fits, libraries, or native objects.
- Leave historical result bodies intact when clarifying status; use dated
  banners and links. Do not infer completion from an old running-session log.
