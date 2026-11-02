# Tariff counterfactuals from reusable posterior draws

Date: 2026-09-18. This extension adds deterministic policy bookkeeping and a
public posterior-realization API. It does not change Stan, demand, conduct,
priors, sampling defaults, ownership semantics or equilibrium tolerances.

## Analytically established

For known, action-independent seller retention, uniform within each active
firm, write \(r_f=1-\tau_f>0\). Then

\[
\pi_f^\tau=\sum_{j\in f}(r_fp_j-c_j)q_j
          =r_f\sum_{j\in f}(p_j-\kappa_j)q_j,
\qquad \kappa_j=c_j/r_f.
\]

Positive scaling preserves each firm's best-response set. Multiplication of
follower FOC rows by their respective retention factors cancels from the
implicit reaction Jacobian. Each noncooperative leader's reduced objective is
multiplied by its own retention factor. This argument includes multiproduct
firms and several simultaneous leaders; it is not a joint-profit objective.

An existing posterior is reusable when demand was fitted to consumer prices,
the margin likelihood uses net seller-revenue margins, recovered costs have
the effective-cost interpretation, and no incompatible physical-cost prior,
equation or linkage was imposed. Algebraically,

\[
m_j^{net}=\frac{r_fp_j-c_j}{r_fp_j}=\frac{p_j-\kappa_j}{p_j},\qquad
c_j^{pre}=r_j^{pre}\kappa_j^{pre},\qquad
\kappa_j^{post}=\kappa_j^{pre}(1+d_j)
                 \frac{r_j^{pre}}{r_j^{post}}.
\]

The margin **error model**, as well as its noiseless accounting formula, must
be compatible. Scaling posterior costs cannot repair a fit to a different
reporting scale. At zero baseline tariff the cost and margin distinctions
coincide. Consumer-price tariff fraction \(\tau\) equals \(t/(1+t)\) for a
producer-value tariff rate \(t\); the benchmark's \(\tau=.2\) corresponds to
\(t=.25\).

The full assumptions and deferred heterogeneous-tariff extension are in
[the implementation contract](tariff_posterior_reuse_plan.md). The production
cost state is immutable across scenarios. Quality and ordered-k realizations
are held fixed between each draw's pre/post equilibria.

## Architecture implemented

`antitrustBayes::posterior_policy_models()` realizes ordinary `AntitrustFit`
objects and native coordination core-fringe/Stackelberg objects from a reusable
`posterior_policy_draws()` plan. It preserves market-period and posterior draw
indices, local quality provenance, exact-k realizations and failure positions.
Default product labels are made unique without reordering product rows. There
is no tariff equation or hard trade dependency in antitrustBayes.

trade adds `TariffGameFit`, extending `TradeFit` through `StructuralFit`.
`as_trade_fit()` accepts the four native coordination classes and an explicit
basis contract for ordinary Logit/CES Bertrand, Cournot and MonCom models.
Promotion loads the existing baseline and adds zero equilibrium solves. Changed
policy delegates through public antitrust/coordination simulation methods;
there are no copied Stackelberg equations or private cross-package calls.
Historical ordinary promotion calls without a basis argument retain their
existing path.

`calcMC()` reports physical costs. `calcMargins(level=TRUE)` reports physical
unit profit \(rp-c\); the accounts also expose effective markup \(p-\kappa\).
`tariff_accounts()` separates physical producer profit and government revenue.
`tariff_welfare()` reports consumer, producer and government changes separately.
It retains the established antitrust CV convention: monetary consumer change
is \(-CV\). No exact national/global welfare or automatic revenue rebate is
asserted. Existing CES consumer-variation limitations remain.

Named tariff/cost vectors are aligned to product labels. Firm-uniformity is
checked using relative retention, including near-100% tariffs. Post checks use
active ownership; exited products contribute zero quantity and zero monetary
accounts. Native ownership changes require caller-supplied post role sets.
Empty core sets are valid. Repeating a scenario does not compound tariff or cost
adjustments. An omitted cost shock retains the current shock; explicit zero
restores original physical costs. Known contradictory cost/margin provenance
cannot be overridden by a caller declaration at nonzero baseline tariff.

For a compatible estimated **output** fit, select the desired market/period
through the existing selectors, realize a plan once, and reuse each returned
model across scenarios:

```r
plan <- antitrustBayes::posterior_policy_draws(fit, n_draws = 200, seed = 123)
prepared <- antitrustBayes::posterior_policy_models(
  fit, draw_plan = plan, market = chosen_market, year = chosen_year,
  demand_anchoring = "latent_quality"
)
stopifnot(nrow(prepared$failures) == 0)
model <- prepared$models[[1]][[1]]
baseline <- trade::as_trade_fit(
  model, tariffPre = tau_pre,
  cost_basis = "effective", margin_basis = "net_revenue"
)
scenario <- trade::simulate(baseline, tariffPost = tau_post, mcDelta = .10)
trade::tariff_accounts(scenario)
trade::tariff_welfare(scenario)
```

Here `tau_pre` and `tau_post` use the consumer-price fraction convention.
Counterfactual ownership changes additionally require the caller's appropriate
`leadersPost` or `corePost` choices for native models. Inspect failures and
preserve the unit/draw indices when aggregating multiple returned models.

## Numerically verified

The standalone standard-library reference in `tariff_uniform_oracle.py` passed
88/88 payoff, FOC, conditional best-response and cost/margin checks. Root
independently reran these checks. Maximum FOC scaling error was
\(1.11\times10^{-11}\); maximum conditional best-response action difference was
\(3.09\times10^{-8}\). A deliberately heterogeneous multiproduct tariff case
gave an action difference of 0.1201, confirming that productwise cost scaling
is insufficient there. These conditional response checks alone do not establish
a fully solved leader/follower equilibrium.

Separately, 16 analytic single-product follower FOCs for Logit/CES price and
quantity games were checked against independently differentiated actual
profits. All passed; maximum absolute error was
\(8.19\times10^{-10}\) against a scale-adjusted tolerance of \(10^{-7}\).

The production equivalence matrix exercises 14 games: ordinary B/C/MC and
native BF/CF/BS/CS, each under Logit and CES. It compares four policies using
nonzero baseline tariffs and a multiproduct firm: unchanged policy, ownership
change, tariff/cost change and combined intervention. Additional tests cover
multiproduct followers, named vector alignment, retained shocks, explicit cost
reset, source immutability, missing margins, exited products, typed rejection
of unsupported cases and zero-solve promotion/unchanged policy.

The final antitrustBayes suite passed 2,212 assertions with zero failures or
errors. It had 13 pre-existing opt-in/dependency skips and 17 fixture warnings;
wall time was 39.21 seconds. All 14 ordinary/native realization combinations
were covered with provided-parameter fixtures. These fixtures do not constitute
estimated CES/exact-k posterior validation.

The final trade suite passed 1,768 assertions with zero failures or errors in
105.43 seconds. It had 123 fixture warnings and one skipped constructor-based
fractional-ownership fixture. A separate explicit fractional-ownership rejection
test ran successfully, so this skip does not replace the rejection coverage.
The unchanged trade baseline had 1,189 passing assertions, zero failures/errors,
and 71 warnings. The added count is test coverage, not a performance comparison.

The corrected full nested reference in `tariff_stackelberg_oracle.py` passed
32/32 checks for Logit/CES Bertrand/Cournot: a multiproduct leader and two
single-product followers. At every leader trial it preserves the candidate
leader actions and resolves the follower equilibrium before evaluating the
leader's actual reduced payoff. Returned follower FOC residuals were at most
\(1.19\times10^{-11}\) (gate \(10^{-10}\)); finite-difference leader FOC residuals
were at most \(4.18\times10^{-6}\) (gate \(10^{-5}\)). Physical/effective leader
actions differed by at most \(1.34\times10^{-8}\), follower actions by
\(1.73\times10^{-10}\), common-profile reaction derivatives by
\(2.23\times10^{-12}\), and scaled leader gradients by
\(1.51\times10^{-13}\). Root reran immutable copies with both the optimized
reference inverse and the original published 180-step bisection inverse;
both passed all 32 checks. The committed script uses the latter published
primitive. These small stationary-equilibrium checks do not prove general
existence, uniqueness or global optimality.

## Saved-posterior policy evidence

Two existing usable joint-Gaussian, collapsed market-quality output-Bertrand
fits were reused, with seeds 314159 and 271828. Each has three markets, two
periods and six products per market-period (36 observations). Their earlier HMC
budgets were four chains with 2,000 iterations, including 1,000 warmup. Both
had zero divergences and zero treedepth hits; maximum R-hat was below 1.003 and
minimum BFMI exceeded 0.72. These are **pre-existing sampler results**, not new
HMC runs for this extension.

For each fit, the same 200 selected posterior draws were used for both
observed-share and latent-quality anchoring and four interventions: no-cost
merger, tariff introduction, +10% physical-cost shock, and merger plus tariff
plus cost shock. Each run made 9,600 actual policy evaluations and 9,600
equivalent-game comparisons. Across both runs, all 19,200 policy evaluations
succeeded. Maximum price difference was \(8.89\times10^{-16}\), share difference
\(5.00\times10^{-16}\), physical-profit accounting error
\(1.25\times10^{-16}\), and government-revenue accounting error
\(2.78\times10^{-17}\). The source fits, plans, realized qualities and draw
provenance remained unchanged.

The complete runs took 169.77 and 161.84 seconds, with peak process RSS of
836,848 and 873,872 KiB. Times include saved-fit loading, plan/structural
construction, promotion, both actual and reference simulations, accounting and
checkpoint output; they are not estimation times. Construction and promotion
are shared across four scenarios, so their repeated CSV columns must not be
summed four times.

These examples demonstrate posterior reuse and exact deterministic
equivalence for a validated output-Bertrand panel. They do not establish tariff
effects for PNB deposit inputs, empirical welfare rankings, or general
insensitivity to latent quality. Paired anchoring comparisons and per-draw
policy intervals are provided in compact CSVs. Market size defaults to one;
welfare is in those normalized units, not industry dollars. Pooled product
quantiles are distinct from posterior intervals for the average price change.

A further five-draw functional check reused an existing usable Gaussian
**sampled** market-quality fit (seed 314159) through both anchoring modes and
all four interventions: 240 policy cases and 240 reference comparisons, zero
failures, 13.05 seconds and 828,432 KiB peak RSS. This verifies the sampled-quality
realization path; it is not a new sampled-versus-collapsed HMC comparison.

## Performance evidence

A paired five-draw pilot reduced complete benchmark wall time from 19.05 to
14.23 seconds by promoting each baseline once and reusing it across scenarios
(1.34 times faster). Maximum RSS changed from 804,324 to 824,568 KiB. Identical
plans and results were checked. This is a small orchestration benchmark, not
an ESS/sec result or a general 200-draw speed estimate.

No likelihood, solver tolerance, prior, sampler setting or role support was
adjusted for performance. Underlying nonlinear iteration/evaluation counts are
not exposed consistently by the existing simulators and are not invented.
Promotion and unchanged-policy zero-solve behavior have direct regression
tests. Changed policy uses the normal public simulation APIs.

## Posterior targets and bridge sampling

There is no new sampling target or marginal-likelihood approximation. The Stan
sources were not edited, and testing reused the unchanged compiled Stan library
(SHA-256 `50b8a7ee6b3d32a82df3c5368962f219f614247d5606e4c6900bf1011adbad7a`).
The helper does not mutate `stanfit`, structural draws or cached-fit version
metadata. Fully normalized existing bridge densities remain the same. No new
bridge computation is claimed here. Physical costs are a deterministic
pushforward of compatible posterior effective costs; this does not assert that
an existing effective-cost prior equals a different prior on physical costs.

## Retained failures and limitations

The run ledger retains failed harness/install attempts, reporting/role/exit
defects subsequently fixed, a metadata partial-matching defect, and rejection
of a stale version-9 cached fit by the unchanged version-10 compatibility
gate. No fit metadata was modified to force reuse. All final suite failures
must be resolved before integration.

The independent full Stackelberg reference initially discarded candidate
leader actions, producing spurious passes; root rejected them. Its final
follower residual initially also evaluated a freshly optimized single-firm
response instead of the returned joint follower profile. That is an invalid
convergence diagnostic. Failed intermediate checks and a defect ledger are
retained. The original spurious-pass CSV was overwritten before it could be
archived; that loss is disclosed rather than reconstructed. Final reference
claims require corrected action propagation, actual returned-profile residuals
and independently re-solved reactions at unchanged strict tolerances.

No current-compatible estimated output CES/exact-role fit was available for
this run. Those families have deterministic and provided-parameter integration
tests; a new empirical posterior exercise remains separate work. Older caches
and input PNB fits are not substitutes.

Initial exclusions are input/deposit markets, auction, bargaining, smooth
fringe, quotas/capacity and within-firm heterogeneous tariffs. Nonzero baseline
tariffs require explicit effective-cost/net-margin declarations and compatible
actual estimation inputs. Declarations cannot certify the original data.
The initial percentage-cost adapter also requires finite, strictly positive
baseline and post costs. Zero/negative states fail explicitly, although the
payoff-scaling proof itself does not require positive costs. Preserve failed
draw identities and rates; do not silently censor them from policy summaries.

## Integration and recommendation

The Bayesian helper is in the antitrustBayes main workspace. The trade
extension is prepared against trade/refactor base
`6609e28e995d46f85c59ea91abae85df3fe164f9`. Automatic application after escalation
failed with `Read-only file system`; the sibling checkout was checked at that
point and had no source changes. The user subsequently applied the patch
locally. All 16 resulting files match the tested candidate's SHA-256 values
after CRLF/LF normalization. The pre-existing untracked promotion patch is
preserved. The active branches are trade/refactor and antitrustBayes/main.
The final 102,575-byte patch is
`doc/statistical_toolkit_runs/tariff-trade-reuse.patch`, SHA-256
`f62180498e8492267504928d1bc88fe201e580de147c4317863ecde377701fec`.
`git apply --check --whitespace=error-all` passed against the unchanged target.
After the final full suite, one trailing space was removed and help/documentation
clarified the existing positive-cost restriction. Identical parsed R syntax
trees and valid Rd syntax were verified. No remote push or new commit is claimed.

Manual command used to complete trade integration:

```bash
cd ~/Projects/trade
git apply ../antitrustBayes/doc/statistical_toolkit_runs/tariff-trade-reuse.patch
```

Source, documentation and compact CSVs are small. Posterior fits, local package
libraries, compiled binaries and full logs remain in the ignored local run
directory, outside the proposed Git changes.

Tests used the already validated installed antitrust 0.99.40 and coordination
0.0.0.9004. The sibling antitrust/refactor source is still 0.99.39, whereas the
pre-existing trade checkout already requires 0.99.40. A pending promotion patch
may belong to another session and is preserved. Resolve this dependency/source
mismatch before installing directly from those sibling sources; this task does
not silently apply an unrelated antitrust promotion.

Use the restricted effective-cost route now for compatible fitted output games.
Retain typed failures for unsupported economics. A generic tariff-payoff solver
for heterogeneous post tariffs is a later extension; it should remain in trade
and reuse coordination's strategic primitives.
