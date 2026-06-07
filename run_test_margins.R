devtools::load_all()
data(pnb)

cat("Margin range:", range(pnb$margins), "\n")
cat("L >= 1 count:", sum(pnb$margins >= 1), "\n\n")

# Logit should work fine with L > 1
cat("=== Logit + Bertrand (should work) ===\n")
fit <- tryCatch(
  bayes_est(
    prices = pnb$prices, quantities = pnb$quantities,
    margins = pnb$margins, market = rep("PNB", nrow(pnb)),
    group = pnb$group, demand_model = "logit", supply_model = "bertrand",
    stan_options = list(chains = 1, iter = 100, warmup = 50, seed = 1)
  ),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit)) cat("OK\n\n")

# CES should error with L > 1
cat("=== CES + Bertrand (should error) ===\n")
tryCatch(
  bayes_est(
    prices = pnb$prices, quantities = pnb$quantities,
    margins = pnb$margins, market = rep("PNB", nrow(pnb)),
    group = pnb$group, demand_model = "ces", supply_model = "bertrand",
    stan_options = list(chains = 1, iter = 100, warmup = 50, seed = 1)
  ),
  error = function(e) cat("ERROR (expected):", conditionMessage(e), "\n")
)
