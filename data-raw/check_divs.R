.libPaths(c("d:/Projects/antitrustBayes/tmp_lib", .libPaths()))
library(antitrustBayes)
data(pnb)

models <- c("bertrand", "monopolistic", "cournot", "auction")
for (sm in models) {
  cat("=== Estimating:", sm, "===\n")
  fit <- bayes_est(
    prices = pnb$prices, shares = pnb$shares, margins = pnb$margins,
    market = rep("PHL", nrow(pnb)), group = pnb$group,
    supply_model = sm, use_rho = FALSE,
    stan_options = list(chains = 4, iter = 2000, warmup = 1000, cores = 4,
      seed = 1960, control = list(adapt_delta = 0.95, max_treedepth = 15))
  )
  sp <- rstan::get_sampler_params(fit@stanfit, inc_warmup = FALSE)
  divs <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  td <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= 15)))
  cat(sprintf("  %s: %d divergences, %d max_treedepth\n", sm, divs, td))
}
