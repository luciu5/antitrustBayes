# Pilot diagnostic availability: no-HMC precision audit

> Closeout update, 2026-10-20: the audit below is retained as sampler-quality
> evidence. It does not gate normalized evidence calculations. Public Nutpie
> and mixed-backend bridge support is implemented in version `0.2.0.9006`;
> see the [current closeout results](nutpie_closeout_results.md).

Executed 2026-10-20 after all nine pilot workers completed. This audit changes
the reference used to recognize exactly constant generated quantities; it
does not change posterior draws, numerical diagnostics or diagnostic limits.
All nine raw versus audited pass decisions are identical. No backend passes
all three pilot cases, and public Nutpie evidence/BMA remains disabled.

The frozen pilot compared every backend's constant values with the original
CmdStan reference, whose CSV retained eight significant figures. Nutpie
retains native doubles. That incorrectly marked nine fixed core-share fields
unavailable in each pooled-CF Nutpie fit and two in each CES-CS Nutpie fit.
The audit uses the original Nutpie-diagonal reference for both Nutpie
adaptations and retains the original CmdStan reference for CmdStan. The
existing exact-equality rule is unchanged: a field must be exactly constant,
match the reference value exactly and be finite. Nonconstant fields still
require finite R-hat, bulk ESS and tail ESS.

## Actual diagnostics

The fixed pilot used four chains, 2000 warmup plus 1000 retained draws per
chain, seed 419977, two cores, target acceptance .99 and maximum depth 13.
CmdStan used `dense_e` for pooled CF and its diagonal metric otherwise;
Nutpie adaptations remained explicit `diag` and `low_rank`.

| Case | Backend/adaptation | Divergences | Depth hits | BFMI range | Max R-hat | Min bulk ESS | Min tail ESS | Raw → audited unavailable fields | Raw/audited diagnostic pass |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | --- | --- |
| B | CmdStanR/diag_e | 0 | 0 | 0.843–0.957 | 1.006115 | 974.21 | 909.66 | 0 → 0 | Yes / yes |
| B | Nutpie/diag | 0 | 0 | 0.808–0.915 | 1.003215 | 720.29 | 597.64 | 0 → 0 | Yes / yes |
| B | Nutpie/low_rank | 0 | 0 | 0.713–0.892 | 1.009210 | 561.56 | 266.67 | 0 → 0 | No / no |
| CF_pooled | CmdStanR/dense_e | 0 | 0 | 0.729–0.838 | 1.009832 | 602.74 | 607.67 | 0 → 0 | Yes / yes |
| CF_pooled | Nutpie/diag | 0 | 0 | 0.730–0.874 | 1.019696 | 407.30 | 639.97 | 9 → 0 | No / no |
| CF_pooled | Nutpie/low_rank | 0 | 0 | 0.699–0.775 | 1.012634 | 319.40 | 679.28 | 9 → 0 | No / no |
| CES_CS | CmdStanR/diag_e | 69 | 0 | 0.751–0.832 | 1.018522 | 217.96 | 73.49 | 0 → 0 | No / no |
| CES_CS | Nutpie/diag | 4 | 0 | 0.756–0.944 | 1.004814 | 614.67 | 543.47 | 2 → 0 | No / no |
| CES_CS | Nutpie/low_rank | 3 | 0 | 0.755–0.876 | 1.007558 | 700.66 | 394.63 | 2 → 0 | No / no |

The frozen limits are zero divergences/depth hits, BFMI at least .3,
R-hat at most 1.01, bulk and tail ESS at least 400, and no unexplained
unavailable posterior diagnostics. B/low-rank fails tail ESS; pooled
CF/diagonal fails R-hat; pooled CF/low-rank fails R-hat and bulk ESS.
CES-CS fails divergences on every backend, with additional R-hat/ESS failures
as shown. Removing the 22 rounding-related false missingness entries cannot
change those conclusions. No genuine unavailable fields remain after the
audit, but that alone does not establish adequate sampling.

## Identity, preservation and reproduction

The [audit runner](nutpie_pilot_diagnostic_audit.R) sources the immutable
`pilot-accept99/frozen-runner.R` and verifies its hash against the original
pilot manifest. It does not depend on the subsequently corrected pilot
runner. It checks the complete production-R/base-helper runtime and all four
Stan source hashes, then verifies exact source/data/prior/model identity
between every pilot fit, its original CmdStan reference and its selected
same-backend constant reference. It reproduces the raw availability and pass
results exactly before changing reference selection.

All 67 sealed files remained unchanged, including the pilot/original
manifests, nine fits, nine raw results, original reference fits, frozen
runner, production sources and original diagnostic CSVs. The audit reads
existing posterior states and performs no compilation, HMC, evidence
integration or policy solving. It preserves the actual numerical diagnostic
values and only recomputes constant-field availability and its consequences
under the original limits.

Tracked detail is in the [36 chain rows](nutpie_pilot_diagnostic_audit_chains.csv)
and [per-variable availability/reasons](nutpie_pilot_diagnostic_audit_variables.csv).
The completed replay is retained at
`.scratch/nutpie-certification/pilot-accept99/diagnostic-audit-20261020T085046/`,
with its seal manifest and complete `result.rds`. The first replay at
`diagnostic-audit-20261020T084936/` emitted 18 harmless row-name recycling
warnings while constructing tables. Its exact executed source is preserved
as `frozen-audit.R`. Resetting the one-row task labels removed those warnings;
the second replay completed with identical diagnostic values and decisions
and no warnings. All original pilot artifacts remain intact.

```sh
AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE=true \
  Rscript --vanilla doc/validation/nutpie_pilot_diagnostic_audit.R
```

`AB_NUTPIE_PILOT_OUTPUT` selects the existing complete pilot root;
`AB_NUTPIE_DEPENDENCY_LIBRARY` selects a compatible dependency library. This
execution reused the existing nutpieR 1.9.1 / BridgeStan 2.8.0 library and the
read-only antitrust 0.99.42 / coordination 0.0.0.9006 dependency library.
Focused non-sampling tests passed nine assertions, including strict rejection
of a genuine constant change and unavailable nonconstant diagnostics.
