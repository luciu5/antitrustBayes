# Default-model history and promotion discipline

The historical specification used joint Gaussian share/margin disturbances and
explicit quality coordinates. Student-t residuals were subsequently adopted
amid divergence concerns. That history motivates investigation; it does not
prove that the Normal likelihood caused those divergences.

The controlled collapsed-Gaussian experiment dated 2026-09-16 held the Gaussian
statistical model fixed between sampled and analytically collapsed quality.
All seven collapsed Gaussian primary fits passed the frozen diagnostic gate;
none of the seven sampled Gaussian or seven sampled Student-t fits passed.
Poor BFMI and mixing, rather than many divergences, were the main reproduced
failure: each sampled group had only one divergence in total. Thus these
results support a sampled-quality geometry problem, not a universal claim
about the cause of historical divergences. The experiment and retained failures
are documented in [collapsed Gaussian results](collapsed_gaussian_results.md).

Analytical collapse preserves the historical Gaussian prior/likelihood exactly.
Market-specific quality is a subsequent substantive prior change: separate
persistent firm-market effects are centered within relevant markets and share
a common quality scale. It removes the historical global centering dependence
across markets. In one market the prior is exactly unchanged. It also permits
exact market-level ordered-k marginalization with active Gaussian rho after
quality integration.

The intended routine specification is Gaussian shares and margins, analytically
collapsed market-specific quality, and active rho. Student-t sampled-quality
and Gaussian-share/Student-t-margin hybrid specifications remain available.
Restoring Gaussian residuals is justified only if the new mathematical,
sampler, bridge, policy, and PNB gates succeed. Scope alone does not ensure
usable HMC or robust policy inference.

Status: the new quality architecture is under development. The public defaults
have **not** been promoted. See [validation contract](market_quality_validation.md)
for the gates. A future promotion must update function signatures, help pages,
tests, and this status together. Historical serialized fits without scope
metadata denote global quality and must never be reinterpreted as market
quality.

The first two-seed multi-market benchmark with active rho does not reproduce
the earlier collapsed-Gaussian usability success: both collapsed and sampled
Gaussian market-quality fits have divergences, and none of the six initial
Gaussian/Student-t fits passes the complete frozen gate. Collapsing improves
BFMI and effective counts but lowers minimum bulk ESS per wall second in this
comparison. These failures remain evidence and block promotion; the earlier
seven-fit result is not transferred to this different configuration. See
[current sampler evidence](market_quality_hmc_results.md).
