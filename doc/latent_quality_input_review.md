# E input-Logit review (read-only)

Reviewed commit `4c71473` in `/tmp/antitrust-lq-input`; no source edit,
build, install, or HMC.

The sign algebra is consistent with the installed antitrust API. `calcShares`
for Logit uses `meanval + alpha*(price - priceOutside)`; the patch constructs
`meanval = delta - alpha*(p - p0)`, so the latent `delta` is recovered at the
reference prices. Input `calcMargins(level=TRUE)` reports the positive
markdown `mc - p`; the patch uses that residual for both the reference and
solved pre-policy FOCs. It retains `mcPre` through the native and ordinary
routes, and the native route passes `mcDelta` into coordination's proportional
`mcPost = mcPre*(1+mcDelta)` path. This supports the stated value/procurement
primitive shock interpretation.

The dispatch guard correctly admits input Logit, output Logit, and output CES,
and excludes input CES and auction/bargaining. Both ordinary and exact native
routes apply the latent replacement before the corresponding simulator; exact
role plans retain the selected k. The independent policy-redteam input
reconnaissance (archived separately) supports the sign/FOC algebra for
nonzero local quality shifts across ordinary B/C and native Stack/core-fringe,
but it predates this public-test extension and uses a local helper copy.

The initial review identified a gap: the first public input tests called `.latent_policy_fixture()`,
which returns `BayesSpecifiedFit`; its quality mode is `none`, so the public
input tests exercise zero quality. The only new nonzero quality call is the
private `.lq_replace_demand()` output-oriented generic test. Add one public
input Logit fixture with a nonzero sampled or collapsed quality draw, verify
that latent pre shares differ from observed shares, solve both pre and post
FOCs, preserve the value primitive without a shock, and verify the documented
proportional value shock changes the post primitive while leaving pre fixed.

This gap is now resolved by the independent
`tests/testthat/test-latent-quality-input-nonzero.R` test. Its 76 assertions
passed for public input Logit Bertrand and Cournot calls using nonzero sampled
quality draws, reusable posterior indices, fixed pre-policy values, and a
10% post-policy value shock. It verifies solved pre/post FOCs and changed
latent baseline prices. These are executed interface/equilibrium tests with
mocked posterior draws, not input-market HMC evidence.

The review also identified a wording defect: `.lq_replace_demand()`'s unsupported-model error said
“input Logit and output CES only” although the same branch supports output
Logit. The error has been corrected to state Logit input/output and output CES.
