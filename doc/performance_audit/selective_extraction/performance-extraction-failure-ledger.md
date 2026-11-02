# Performance benchmark failure ledger

| cell | outcome | interpretation |
|---|---|---|
| Fixed_param N=1,000 first attempt | failed before model construction | The private runner passed `s0_pooling` together with legacy `use_s0_market_pool`/`use_s0_market_year_pool` spellings. The public API correctly rejected the contradictory argument contract. The runner was corrected to omit all s0 pooling spellings; no package source, economics, priors, or settings changed. |
| archived policy fixtures | 36/36 paired rows complete | No extraction, construction, solver, or policy-equivalence failure observed. |
| Fixed_param stress N=1,000 and N=10,000 | 4/4 mode-size cells complete | Same kernel, 200 rows, finite canonical fields, equal full/selected hashes. |
