# Current CmdStanR statistical validation

Run the known-truth recovery gate after installing the current package and
CmdStan:

```sh
AB_RECOVERY_OUTPUT=/tmp/antitrustBayes-recovery \
  Rscript --vanilla doc/validation/known_truth_recovery.R
```

Set `AB_RECOVERY_CMDSTAN` to a **writable** CmdStan installation for the
default full gate, and `AB_RECOVERY_CACHE` to a writable model cache if the
default temporary cache is unsuitable. The script uses fixed seeds, two HMC
chains per fit, 250 warmup and 250 retained iterations per chain. It writes
recovery, evidence, and diagnostic CSVs to `AB_RECOVERY_OUTPUT` and exits
nonzero if any required criterion fails. Use
`AB_RECOVERY_SCOPE=recovery` for the shorter parameter-recovery check.

The fixture is generated without antitrustBayes helpers. It has six markets,
three years per market, and five single-product firms per market. Conditional
Logit shares come from a softmax with known structural price coefficient
`alpha = 2.5`. The outside share is constant across years in each market.
Output margins follow the single-product Bertrand FOC
`1 / Lerner = price * alpha * (1 - (1 - s0) * inside_share)`
with independently generated margin noise. The resulting 90 observations
are fitted with the current public `bayes_est()` CmdStanR path.

The gate requires the true alpha inside its 95% posterior interval and its
posterior median within 8% of truth; at least five of six true outside
shares inside their 95% intervals with mean absolute median error at most
six percentage points; a plausible recovered margin-noise scale; key
rank-normalized R-hat at most 1.05, minimum bulk ESS at least 100, and zero
divergences. It also checks the economic invariants that a common alpha
does not vary across markets and a market-level outside share does not vary
across years in any posterior draw. These thresholds were chosen for a
bounded known-truth diagnostic, not as a claim of repeated-sample coverage.

The full gate also fits Bertrand and Cournot to a second, independently
generated **one-market**, eight-year Bertrand market. Its outside share is
constant across years and a free market intercept is inactive, so Cournot
cannot imitate Bertrand by raising an unconstrained market outside share.
The gate uses the public `bayes_compare()` bridge path with two repetitions.
It requires both fits to have sound chain diagnostics, both bridge
repetitions to be finite and consistent, the current Stan source hash to
match the fit, a Bertrand/Cournot Bayes factor of at least 10, and posterior
model probability for the true Bertrand model of at least 0.9 under equal
model priors. The completed full gate gave log marginal likelihoods 135.627
(Bertrand) and 82.129 (Cournot), with bridge repetition SDs 0.023 and 0.003
log units.

A first multi-market evidence pilot did **not** meet that standard: Cournot
could inflate its market outside shares, and its bridge repetitions differed
by about five log-evidence units. That more general model-selection problem
remains open; the successful one-market check establishes a bounded current
CmdStanR evidence path, not universal conduct identification.

The historical 124-task enterprise benchmark in
`doc/market_intercept_enterprise/` still requires an RStan-to-CmdStanR port.
These bounded checks establish current-backend parameter recovery and
model-evidence preference in two identified designs; they do not establish
broad repeated-sample coverage or replace the enterprise study.
