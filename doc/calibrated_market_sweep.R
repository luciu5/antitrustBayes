# Reproducible J=7 calibrated-market sweep. Run from the package root with:
# R_LIBS=/path/to/dependencies Rscript --vanilla doc/calibrated_market_sweep.R
#
# Each conduct/seed cell is attempted exactly once. Failures are captured in
# the CSV; this script intentionally performs no redraw or replacement.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The sweep script requires devtools to load the package source.")
}
devtools::load_all(".", quiet = TRUE)

conducts <- c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")
seeds <- 9201:9210
J <- 7L
dirichlet_concentration <- 2.5
margin_range <- c(0.08, 0.16)

cells <- expand.grid(
  conduct = conducts,
  seed = seeds,
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

rows <- lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  result <- tryCatch(
    list(ok = TRUE, market = bayes_calibrated_market(
      conduct = cell$conduct,
      J = J,
      dirichlet_concentration = dirichlet_concentration,
      margin_range = margin_range,
      seed = cell$seed
    )),
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
  if (!result$ok) {
    return(data.frame(
      conduct = cell$conduct,
      seed = cell$seed,
      status = "failure",
      anchor_margin = NA_real_,
      realized_anchor_margin = NA_real_,
      alpha = NA_real_,
      anchor_index = NA_integer_,
      min_cost = NA_real_,
      max_price_error = NA_real_,
      max_share_error = NA_real_,
      baseline_foc = NA_real_,
      error = result$error,
      stringsAsFactors = FALSE
    ))
  }
  market <- result$market
  data.frame(
    conduct = cell$conduct,
    seed = cell$seed,
    status = "success",
    anchor_margin = market$design$anchor_margin,
    realized_anchor_margin = market$design$realized_anchor_margin,
    alpha = market$truth$alpha,
    anchor_index = market$design$anchor_index,
    min_cost = min(market$truth$costs),
    max_price_error = market$truth$equilibrium$max_relative_error,
    max_share_error = market$truth$equilibrium$max_share_error,
    baseline_foc = market$truth$equilibrium$max_foc_residual,
    error = "",
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, rows)
results$J <- J
results$dirichlet_concentration <- dirichlet_concentration
results$margin_low <- margin_range[[1L]]
results$margin_high <- margin_range[[2L]]
results <- results[c(
  "conduct", "seed", "J", "dirichlet_concentration", "margin_low",
  "margin_high", "status", "anchor_margin", "alpha", "anchor_index",
  "realized_anchor_margin", "min_cost", "max_price_error", "max_share_error",
  "baseline_foc", "error"
)]

write.csv(results, "doc/calibrated_market_sweep.csv", row.names = FALSE,
          na = "NA")
cat("J=", J, "; concentration=", dirichlet_concentration,
    "; anchor margin interval=", paste(margin_range, collapse = " to "),
    "; attempted=", nrow(results), "; failures=", sum(results$status == "failure"),
    "\n", sep = "")
print(with(results, table(conduct, status)))
