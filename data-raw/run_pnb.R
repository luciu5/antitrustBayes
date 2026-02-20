# PNB merger example using antitrustBayes package
library(antitrustBayes)

data(pnb)
cat("PNB data:", nrow(pnb), "banks\n")
print(head(pnb))

# Estimate under 4 supply models
supply_models <- c("bertrand", "monopolistic", "cournot", "auction")

fits <- lapply(supply_models, function(sm) {
  cat("\n=== Estimating:", sm, "===\n")
  bayes_est(
    prices  = pnb$prices,
    shares  = pnb$shares,
    margins = pnb$margins,
    market  = rep("PHL", nrow(pnb)),
    group   = pnb$group,
    supply_model = sm,
    use_rho = FALSE,
    stan_options = list(
      chains = 4, iter = 2000, warmup = 1000,
      cores = 4, seed = 1960,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
  )
})
names(fits) <- c("Bertrand", "MonCom", "Cournot", "Auction")

# Show Bertrand fit
cat("\n=== Bertrand Fit ===\n")
fits$Bertrand

# Model comparison with bridge sampling
cat("\n=== Model Comparison ===\n")
comp <- bayes_compare(
  Bertrand = fits$Bertrand,
  MonCom   = fits$MonCom,
  Cournot  = fits$Cournot,
  Auction  = fits$Auction
)
print(comp)

# Weights
cat("\nBridge sampling weights:\n")
print(bayes_weights(comp, "bridge"))

# Simulate the PNB-Girard Trust merger
owner_post <- pnb$group
owner_post[owner_post == "889018"] <- "298816"

cat("\n=== Simulating merger under each model ===\n")
sims <- lapply(fits, function(f) {
  bayes_sim(f, ownerPost = owner_post, n_draws = 200)
})

cat("\n=== Bertrand simulation summary ===\n")
summary(sims$Bertrand)
summary(sims$Bertrand, market = TRUE)

# BMA
cat("\n=== Bayesian Model Averaging ===\n")
w <- bayes_weights(comp, "bridge")
avg_sim <- bayes_avg(sims, weights = w)
summary(avg_sim)
summary(avg_sim, market = TRUE)

cat("\nDone!\n")
