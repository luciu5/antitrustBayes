devtools::load_all()

data(pnb)
# CES requires Lerner indices below one; use the common CES/Logit sample.
pnb <- pnb[pnb$margins < 1, , drop = FALSE]
cat("PNB common CES sample:", nrow(pnb), "banks\n")
cat("Margin range:", range(pnb$margins), "\n\n")

supply_models <- c("bertrand", "cournot", "auction", "monopolistic")
model_priors <- bank_model_priors(supply_models)

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
      prior_manifest = "baseline",
      model_prior = model_priors[[sm]],
      stan_options = list(chains = 2, iter = 1000, warmup = 500, seed = 1960)
    ),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fit)) {
    print(fit)
    cat("\n")
  }
}
