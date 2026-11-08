# Explicit Logit price scaling (2026-11-06)

`bayes_est(..., price_standardization = list(center = p_bar, scale = s_p))`
supports a fixed price center and positive scale for relative-share Logit.
The default remains within-sample mean/SD; CES and absolute-share fits reject
a supplied anchor. Metadata retains the actual center/scale and the transformed
prices and scale enter the Stan-data hash. Fit calls retain the explicit option.

For a comparison across information sets, compute the anchor once from the
target baseline snapshot and pass exactly that pair to each fit. Since
`a = alpha * s_p`, a common prior on log `a` then induces a common prior scale
for the economic slope. Centering cancels in relative-share contrasts. This
changes the economic prior relative to within-sample scaling; use new run/cache
identities and do not relabel existing inference. Pooling and parameter-count
differences across sets remain.

`tests/testthat/test-price-standardization.R` verifies invariant native slope
and inverse-wedge conversions, common snapshot coordinates, unchanged default
scaling, prior parameters, and rejection of invalid/unsupported anchors.
The 18 source-level assertions passed using the paper's existing shared library.
Full `testthat::test_local(filter = "price-standardization")` could not load
the development package: installed antitrust/coordination versions are below
its dependency minima. The paper's authorized source refresh failed GitHub
DNS resolution even after escalation. No inference was run, so these checks
do not establish posterior/evidence or policy robustness.
