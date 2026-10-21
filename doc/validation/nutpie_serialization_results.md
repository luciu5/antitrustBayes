# Nutpie saved-fit operation without native dependencies

Executed 2026-10-20 from branch `nutpie`, Git
`8d6daa115bf2fb7888a26d1c3ca4ba6666355c64`, package `0.2.0.9005`.
The no-HMC audit started at 08:39:13 EDT and passed for five actual saved
Nutpie fits from the completed paired benchmark. No fitting, compilation or
native binary loading was needed for posterior and policy operations.

The [audit runner](nutpie_serialization_validation.R) starts an independent
`Rscript --vanilla` process for each fit. Before process startup it sets
`R_LIBS`, `R_LIBS_USER` and `R_LIBS_SITE` to libraries excluding every directory
containing `nutpieR` or `bridgestan`. Each child confirms that both namespaces
are unavailable before package load, after package load and after all checks.
It loads the development antitrustBayes source and existing compatible
antitrust/coordination dependencies; it does not install or modify packages.

| Fixture | Kernel | Retained posterior states | Exact result fields | Merger solves | Tariff solves | Maximum physical-cost identity error | Result |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| B | logit_supply | 4000 | 10 | 2 | 2 | 0 | Passed |
| BS | logit_stackelberg | 4000 | 10 | 2 | 2 | 0 | Passed |
| CES_B | ces_supply | 4000 | 10 | 2 | 2 | 0 | Passed |
| CES_CS | ces_stackelberg | 4000 | 10 | 2 | 2 | 0 | Passed |
| CF_pooled | logit_stackelberg | 4000 | 10 | 2 | 2 | 0 | Passed |

The [compact CSV](nutpie_serialization_results.csv) records each input fit's
SHA-256 and every outcome. The ten exact comparisons are the complete posterior
array hash, public metadata hash, native metadata hash, indexed parameter
extraction, posterior summaries, common diagnostics, reconstructed quality,
posterior policy draw plan, reconstructed baseline class/accounts/provenance,
and solved policy prices. Quality and policy plans use two posterior states
with fixed seeds 721091 and 721092 respectively. All 4000 posterior states
are retained for array comparison, parameter summaries and diagnostics.

The merger joins the two largest observed firms. The tariff is 25% for the
other firms, represented by the existing tax-inclusive wedge `.25 / 1.25`.
Both use latent-quality reconstruction, the same fixed plan, inherited roles
and existing trade solvers. The audit compares native-enabled reference
results with optional-free child results exactly; it does not compare
independently sampled posterior distributions or establish statistical
recovery. Each child also repeats the posterior and policy calculations with
an absent native-library path on a copied fit. All results other than that
deliberately changed native metadata hash remain exactly equal.

## Serialization and native-artifact boundaries

Every non-call S4 slot is identical after `saveRDS()`/`readRDS()`. These fits
were constructed through `do.call()`, which retained a function closure in
the saved call. Its source-reference environments acquire new pointer
identities on reload, so whole-object `identical()` is too strict. The
serialized call digest is identical and `all.equal()` of the whole object
is true. This is an R call/source-reference identity detail; posterior,
data, scientific metadata and native draw records are bit-identical.

Evidence restoration retains its separate dependency and artifact contract:

| Condition | Observed rejection |
| --- | --- |
| Optional-free process | `Nutpie bridge target evaluation requires the optional 'bridgestan' R package.` |
| Native-enabled process, missing copied library path | `The saved Nutpie BridgeStan shared library is unavailable; normalized target evaluation cannot be restored.` |
| Native-enabled process, one byte appended to a workspace copy of the library | `The Nutpie BridgeStan shared library hash differs from its saved provenance.` |

The changed binary was never loaded. The original native library and all 38
sealed input-fit/production-R/audit-source files retained their exact SHA-256
values. No original fit, validation result or existing benchmark manifest
was rewritten. The public evidence/BMA gate remains closed.

## Reproduction and environment

Raw requests, references, child results, logs, seals and artifact rejections
are retained under `.scratch/nutpie-certification/serialization`. The runner
is explicitly opt-in and uses only existing saved fits:

```sh
AB_NUTPIE_SERIALIZATION_EXECUTE=true \
  Rscript --vanilla doc/validation/nutpie_serialization_validation.R
```

For each `<fixture>` in the table, that output root contains
`<fixture>/request.rds`, `reference.rds`, `reloaded-fit.rds`, `result.rds` and
`worker.log`. The child `result.rds` records its effective `.libPaths()`, loaded
namespace names, ten individual equality checks and complete numerical outputs.
The root's `manifest.rds` records all 38 input/source seals and the startup
library sets; `artifact-negatives.rds` records the two native-enabled rejection
messages and unchanged-source/library checks. Inputs remain in
`.scratch/nutpie-validation/<fixture>__nutpie_diag/fit.rds`.

Its dependency library can be selected with
`AB_NUTPIE_DEPENDENCY_LIBRARY`. This execution read the existing compatible
paper library as a dependency source only. Resolved versions were R 4.6.1,
antitrust 0.99.42, coordination 0.0.0.9006, trade 0.8.9.9002, posterior 1.7.0,
CmdStanR 0.9.0.9002 and pkgload 1.5.2. The reference/negative-artifact process
additionally had nutpieR 1.9.1 and bridgestan 2.8.0; the child processes had
neither available. Focused ordinary runner tests passed eight assertions
with zero failures or warnings.

This closes the bounded optional-free saved-fit functionality check across
all four kernels and one pooled configuration. It does not establish native
artifact portability between machines or operating systems. No Stan kernels,
priors, likelihoods, sampler controls or mergerBayes production analyses were
changed.
