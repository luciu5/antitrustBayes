# Nutpie CI readiness, 2026-10-20

## Post-merge hosted run, 2026-10-20 EDT / 2026-10-21 UTC

The first [hosted run on main at `ab5878a`](https://github.com/luciu5/antitrustBayes/actions/runs/39532459744)
installed the pinned dependencies and passed the optional-free/provenance
checks, then failed package tests because the Stan parsing prerequisite called
`cmdstanr::cmdstan_path()` before handling an unconfigured CmdStan installation.
The ordinary job intentionally does not install CmdStan. Its installed tests
reported 4,927 passes, 21 existing warnings, 29 skips and this one error.

The follow-up catches path-discovery errors and uses the existing explicit
missing-toolchain skip. All source assertions and four actual parser checks
are unchanged. Focused local verification passed 44 assertions with CmdStan
2.39.0 and 40 assertions plus the intended single skip with path discovery
mocked to throw the hosted error; neither route had warnings or failures.
An independent reviewer found no blocker. Evidence is retained in
`.scratch/nutpie-closeout/ci-prerequisite-check.log` and its RDS summary.
No sampling was rerun and no production source was changed. A hosted rerun
after pushing this fix remains pending; clean native Nutpie CI remains the
separate future task below. The original readiness assessment follows.

This follow-up starts from `nutpie` commit `8d6daa1`. It fixes dependency-stack
reproducibility and adds explicit coverage of an installation without Nutpie.
It does not establish statistical certification or enable evidence/BMA.

## Verified dependency sources

The old CI antitrust pin `103c2491607d7f164f5e4e00fc21d6ea79331138` is
version 0.99.39. Coordination 0.0.0.9006 requires antitrust at least 0.99.41,
so the old pin was incompatible with the stack used for the backend validation.
The workflow now uses these independently checked sources:

| Package | Version | Git revision | DESCRIPTION Git blob |
| --- | --- | --- | --- |
| antitrust | 0.99.42 | `20f0923fc560295e04310aed5615c19dadf672cb` | `7cef63f6a951a91a4847544c6f9e37fec4340951` |
| coordination | 0.0.0.9006 | `6955c4845df0c1de98eb540423118b53732211e5` | `aa5e2ad582a9123dddcd93f79bb42e33df547927` |
| trade | 0.8.9.9002 | `2826df0fc3530c01066000725f643b060e71d136` | `f8f9c5cceb1d93a4b270477698f214981c17ade9` |
| cmdstanr | 0.9.0.9002 | `55342f4f6c9300f9a47046f107b4db9b5e468aa6` | `aa7dc6897009f5e9d3d3209a8f7851fa4903988c` |

The [antitrust](https://github.com/luciu5/antitrust/blob/20f0923fc560295e04310aed5615c19dadf672cb/DESCRIPTION),
[coordination](https://github.com/luciu5/coordination/blob/6955c4845df0c1de98eb540423118b53732211e5/DESCRIPTION),
[trade](https://github.com/luciu5/trade/blob/2826df0fc3530c01066000725f643b060e71d136/DESCRIPTION),
and [CmdStanR](https://github.com/stan-dev/cmdstanr/blob/55342f4f6c9300f9a47046f107b4db9b5e468aa6/DESCRIPTION)
descriptions were inspected at their exact revisions through read-only GitHub
access. Local structural checkouts have those tracked HEAD revisions; an
unrelated untracked trade patch was preserved. The installed validation
libraries have matching package versions. Directly reading the validation
library's installed antitrust, coordination and trade `DESCRIPTION` files
returns the exact `RemoteSha` values in the table. However,
`utils::packageDescription()` returns no `RemoteSha` for these three installed
packages; the frozen pilot runtime records that unavailable API metadata.
These two metadata views are recorded separately rather than rewriting the
saved provenance or inferring source identity from a version alone. CmdStanR's
installed `packageDescription()` `RemoteSha` matches the table. The new CI
installs from verified Git checkouts and records their actual SHAs directly.

Coordination's minimum is antitrust 0.99.41. Trade requires antitrust 0.99.41
and coordination 0.0.0.9005; the pinned stack satisfies both. Structural
packages are installed in dependency order from Git checkouts whose actual
HEAD must equal each expected full SHA. Installation success and the installed
version are checked before testing. The verified checkout SHAs are recorded
in the job summary. No later dependency resolver can replace them with the
moving `antitrust@refactor` reference in coordination's `Remotes` field.

## What the workflow now checks

The ordinary job uses R 4.6.1 and explicitly installs hard dependencies plus
the selected test/vignette tools. `setup-r-dependencies@v2` receives a complete
`packages` list, which replaces its default `deps::.` request, and
`dependencies = c("Depends", "Imports", "LinkingTo")`. Its
[action implementation](https://github.com/r-lib/actions/blob/v2/setup-r-dependencies/action.yaml)
passes those two inputs directly to the dependency lockfile calculation.
This avoids implicitly installing every package in `Suggests`. The explicit
list covers all current non-base antitrustBayes imports and all non-base
imports of the selected structural packages; hard transitive dependencies are
resolved normally. Future import changes must update this list.

The job asserts that neither `nutpieR` nor `bridgestan` is available, then
builds and checks the package with `_R_CHECK_FORCE_SUGGESTS_=false` and sampling
opted out. All other declared suggestions are installed explicitly. Expected
skips are the native Nutpie dependency check and opt-in sampling blocks;
these remain visible in test output. A separate source-based backend/runner
test filter also exercises research scripts excluded from the built package,
including the CES Stackelberg tariff regression. Trade is a pinned CI test
dependency, not a new package dependency.

The existing scheduled/manual CmdStan job retains its Logit/CES/pooling/
dispatch filter. It now installs and verifies CmdStan 2.39.0 before opting
into sampling, so a missing CmdStan installation cannot make that tier appear
successful solely through dependency skips. CmdStanR's
[documented installation interface](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)
provides the explicit version and core controls used here. These are the same
CmdStan/CmdStanR versions used for the executed validation. The fast and Stan
jobs have 45- and 90-minute limits respectively. Package-check artifacts are
uploaded from the fast job even after a failure.

Local verification parsed the YAML, all embedded R code and every shell block;
`bash -n` passed. The two dependency lists are identical, cover current imports,
and exclude both native optional packages. The optional-free installed-library
and CmdStanR provenance guards passed using the existing compatible library.
Independent review found no substantive CI correctness issue. No hosted job,
dependency installation or HMC was executed for this CI edit. A successful
hosted workflow run remains required evidence; syntax/provenance checks do not
stand in for that run. `actionlint` was unavailable locally.

## Remaining native Nutpie CI work

A native Nutpie job is deliberately not claimed as ready. The executed local
validation used nutpieR 1.9.1 at
`8e88a975d5e0b8ed3c76de84d22adb33de601a1e`, nuts-rs 0.18.2, and BridgeStan
2.8.0. NutpieR's exact
[Cargo lockfile](https://github.com/andytimm/nutpieR/blob/8e88a975d5e0b8ed3c76de84d22adb33de601a1e/src/rust/Cargo.lock)
and native binary paths independently identify the nuts-rs version. The local
BridgeStan source setup required recovered, hash-verified files and reuse of
the existing Stan toolchain after network failures. That recovery is not a
validated clean-runner installation recipe.

The next native CI task is a separately opted-in, clean-runner smoke job with
pinned nutpieR/BridgeStan artifacts, a verified BridgeStan source download,
recorded compiler/native ABI versions, and a bounded four-kernel sampling/
serialization/generated-quantity check. A prebuilt binary route can avoid a
Rust build where a compatible binary exists; artifact availability and hashes
must first be verified for the selected runner/R combination. NutpieR's
[installation instructions](https://github.com/andytimm/nutpieR/tree/8e88a975d5e0b8ed3c76de84d22adb33de601a1e)
distinguish binary installation from source installation with Rust.

Neither workflow currently supplies that native Nutpie tier. No claim is made
about Windows/macOS support, clean-runner native reproducibility, improved
sampling geometry, or public evidence certification from these CI changes.
