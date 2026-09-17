# Native build checkpoint

Recorded 2026-09-17 after the two matched native installs completed.

**September 17 recovery correction:** the checksum of the backup labeled
Helmert below was preserved, but loading it reports all four rotation-era
model-code digests. Its installed Stan source also contains the dense quality
matrix. That package is not accepted as a Helmert build. The historical record
is retained; `native-recovery-failures.csv` records the provenance failure.
A fresh source-verified build is required before the final native gate.

Both roles used R 4.6.1, StanHeaders 2.32.10, rstan 2.32.7, g++ 11.4.0,
`CXX17=g++`, the package's C++17 command, and `MAKEFLAGS=-j1`. The selected
compiler is `/usr/bin/x86_64-linux-gnu-g++-11` with SHA-256
`2360901d864cf10bfd6296e261cb2c14053552a80377761ab07146ec9ec9a2c0`; the
flags hash is `ca488275b6fe249adf16fa2d4dc6b7c9b08a2a1f851828482f423987a70b67cf`.

| role | source hash | installed DLL hash | package |
|---|---|---|---|
| baseline | `0b29fa1cf0de7ce0eb185dd1f73caa771e84e2b670d5fb998262b0b3fe071b8c` | `31d1beca03ac862872b1f96c9943c2b082c1d93b105cb64f4854aef51999afdc` | 0.2.0 |
| Helmert candidate | `fb2287a1fc50e832c9c61fe55c143b830ed25f945ca2c8af3c7526a11ed50176` | `08df2efd81d920c863459be6dfd053fcbdf1dc2eb7178a8b508535eee2fe44ce` | 0.2.0 |

The full 15 MB install logs are retained at
`doc/statistical_toolkit_runs/20260916/performance-audit/recovery/logs/`.
The durable libraries are in the ignored recovery `Rlib/baseline` and
`Rlib/candidate` directories. Atomic per-role records are
`stage-baseline.complete` and `stage-candidate.complete`; each records source,
compiler, flags, dependency-manifest, and installed-DLL hashes.

The original harness was
`doc/performance_audit/recovery/run_native_install.sh`. A subsequent dry
restart completed with both roles reported `action=SKIP reason=verified-complete`;
the guard checks all recorded hashes, the DLL, and `requireNamespace()` before
skipping a role. An interrupted role is represented by its `STARTED` record and
does not qualify for skipping.

The repaired harness excludes generated native objects from its source
fingerprint, checks the actual loaded model-code generation, freezes generated
source before compilation and propagates install errors. Its source-only
manifest uses version 2, distinct from the historical whole-tree hashes above.
Interrupted objects can be reused only with identical recorded build inputs
and explicit confirmation that no live writer remains. The small fingerprint
fixture passes 32 assertions; native compilation is a separate gate.

The corrected candidate install now loads all four expected Helmert digests.
Its DLL SHA-256 is
`417e584bd6fb96f7b3aa65860c953fe6e09c8fdf715a8daeca49a50ccc0a1aaa`,
with source-only fingerprint
`8c7746eb1ca3d2ecd7b4c496b9d6f162b7dba461a891e7bdf3068813686d5a87`.
All five compiler-command prefixes equal the preserved baseline exactly;
see `recovered_compiler_comparison.csv`. A stable snapshot of the repaired
harness exits zero and skips the verified completed candidate. The original
driver exited 2 after its successful completion record when it was edited
while executing; that orchestration failure is retained. Future long drivers
must run immutable script snapshots. Default recovery selects only candidate
so the preserved version-1 baseline is not unnecessarily rebuilt; choose
`PERF_RECOVERY_ROLES="baseline candidate"` explicitly for a fresh pair.

The ordinary installed-package suite exits zero with 13 skips and 13 warnings.
The extra skip is the public observed-baseline test because this recovery's
antitrust dependency has the promotion prerequisite but not the separate
lifecycle patch. Seven final fixtures resolve through the unchanged density,
full-gradient, schema, output and centering gates plus the documented
finite-difference cancellation diagnostic. The original 6/7 central gate and
all failures remain available under `native_review/`.
