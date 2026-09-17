# Private bridge review outcome

The isolated package source is based on main `ea29e96a443feab0436b3bbe47d1295e6bd43474`.
Candidate overlays and SHA-256:

| file | SHA-256 |
| --- | --- |
| `R/bma.R` | `eac3b32b21a2f04c4a2f1cc9afb155dc8163cc3e80e713dfa0dedee30811d7f0` |
| `R/model_universe.R` | `1aeeed9be05526856a40d3eb637656463391f83d5bcaed0b8d265064749046d4` |
| `tests/testthat/test-bridge-performance-recovery.R` | `ee670ac1aae841f8e8da07a6216578d5be3b6968b8f3814d70e3c76a7adda38b` |
| `tests/testthat/test-bridge-worker-warnings.R` | `789cc822e369b81705d8773c32ac7017d9e1821c4a0801266be6f613ca38ff53` |

The private `R CMD INSTALL --no-configure --no-libs --no-test-load` install
completed successfully after the recovery library root was added to
`R_LIBS_USER`. The earlier missing-`coordination` and interrupted temporary
load-check attempts are preserved in `logs/install.log` and
`logs/install-recovery-root.log`. The installed package was explicitly loaded
from `Rlib/bridge-reviewed`.

The DLL hash was unchanged at every checkpoint:

`417e584bd6fb96f7b3aa65860c953fe6e09c8fdf715a8daeca49a50ccc0a1aaa`

The loaded four Stan model code hashes equal the registered `helmert` hashes:

```text
ces_stackelberg 5906d0e9f7a00c6ec29da47cb9e870bddb65f8750f53855286b068aa24de73e6
ces_supply      e4a759841324154c27bdfc2508cf0babae3a0daaba0c942aefc1289deaae3592
logit_stackelberg 4e91beeba1f0155781d72a824bda92f9e1ad51b901d583f9700056b2bb3bcc76
logit_supply    59cdf553449bfe192aa8e13398fd54bcaf9b3e991289e2e9fd2afc0e93994ff3
```

The narrow `man/bayes_compare.Rd` update parses successfully and documents
`bridge_cores`, `bridge_seed`, and `bridge_warning_records`. Its before-image
and diff are in `logs/bayes_compare.Rd.before` and
`logs/bayes_compare.Rd.diff`; the installed help database contains all three
names, as recorded in `logs/check-installed-doc.log`.

Focused installed tests all passed:

- `test-bridge-worker-warnings.R`: pass; `logs/test-worker-warnings.log`.
- `test-bridge-performance-recovery.R`: pass; six expected synthetic failure
  warnings retained in `logs/test-performance-recovery.log`.
- `test-bridge-performance-migration.R`: pass; `logs/test-performance-migration.log`.
- `test-bridge-context-restoration.R`: pass; `logs/test-context-restoration.log`.

The full installed suite was outside this bounded review and was not run.
