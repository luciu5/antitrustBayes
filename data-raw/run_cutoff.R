.libPaths(c("d:/Projects/antitrustBayes/tmp_lib", .libPaths()))
library(antitrustBayes)
data(pnb)

models <- c("bertrand", "monopolistic", "cournot", "auction")
model_names <- c("Bertrand", "MonCom", "Cournot", "Auction")

fits <- list()
for (i in seq_along(models)) {
  cat("\n=== Estimating:", model_names[i], "with use_cutoff=TRUE ===\n")
  fits[[model_names[i]]] <- bayes_est(
    prices = pnb$prices, shares = pnb$shares, margins = pnb$margins,
    market = rep("PHL", nrow(pnb)), group = pnb$group,
    supply_model = models[i], use_cutoff = TRUE, use_rho = FALSE,
    stan_options = list(chains = 4, iter = 2000, warmup = 1000, cores = 4,
      seed = 1960, control = list(adapt_delta = 0.95, max_treedepth = 15))
  )
  sp <- rstan::get_sampler_params(fits[[model_names[i]]]@stanfit, inc_warmup = FALSE)
  divs <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  td <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= 15)))
  cat(sprintf("  %s: %d divergences, %d max_treedepth\n", model_names[i], divs, td))
}

cat("\n=== Diagnostics Summary ===\n")
sapply(fits, function(f) {
  sp <- rstan::get_sampler_params(f@stanfit, inc_warmup = FALSE)
  c(divergences   = sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
    max_treedepth = sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= 15))))
})

cat("\n=== Model Comparison ===\n")
comp <- bayes_compare(
  Bertrand = fits$Bertrand, MonCom = fits$MonCom,
  Cournot = fits$Cournot, Auction = fits$Auction
)
print(comp)

cat("\n=== Bridge Sampling Weights ===\n")
print(bayes_weights(comp, "bridge"))

cat("\n=== Stacking Weights ===\n")
print(bayes_weights(comp, "stacking"))

cat("\n=== Key Parameters (cutoff_share) ===\n")
for (nm in model_names) {
  cs <- rstan::summary(fits[[nm]]@stanfit, pars = "cutoff_share")$summary
  cat(sprintf("  %s cutoff_share: mean=%.4f, sd=%.4f, 95%%CI=[%.4f, %.4f]\n",
    nm, cs[,"mean"], cs[,"sd"], cs[,"2.5%"], cs[,"97.5%"]))
}

# Simulate merger
owner_post <- pnb$group
owner_post[owner_post == "889018"] <- "298816"

cat("\n=== Simulating Merger (Bertrand) ===\n")
sim_b <- bayes_sim(fits$Bertrand, ownerPost = owner_post, n_draws = 200)
summary(sim_b)
summary(sim_b, market = TRUE)

cat("\n=== Simulating Merger (Cournot) ===\n")
sim_c <- bayes_sim(fits$Cournot, ownerPost = owner_post, n_draws = 200)
summary(sim_c)
summary(sim_c, market = TRUE)

cat("\n=== BMA ===\n")
sims_all <- lapply(fits, function(f) bayes_sim(f, ownerPost = owner_post, n_draws = 200))
w <- bayes_weights(comp, "bridge")
avg <- bayes_avg(sims_all, weights = w)
summary(avg)
summary(avg, market = TRUE)

cat("\nDone.\n")
