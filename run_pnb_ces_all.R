devtools::load_all()

data(pnb)
cat("PNB dataset:", nrow(pnb), "banks\n")
cat("Margin range:", range(pnb$margins), "\n\n")

supply_models <- c("bertrand", "monopolistic", "cournot", "auction")

for (sm in supply_models) {
  cat("=== CES +", sm, "===\n")
  fit <- tryCatch(
    bayes_est(
      prices       = pnb$prices,
      quantities   = pnb$quantities,
      margins      = pnb$margins,
      market       = rep("PNB", nrow(pnb)),
      group        = pnb$group,
      demand_model = "ces",
      supply_model = sm,
      stan_options = list(chains = 2, iter = 1000, warmup = 500, seed = 1960)
    ),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fit)) {
    print(fit)
    cat("\n")
  }
}
