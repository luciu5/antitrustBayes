devtools::load_all()

data(pnb)

# Both models on same sample: restrict to L < 1 (required by CES)
pnb_sub <- pnb[pnb$margins < 1, ]
cat("Common sample:", nrow(pnb_sub), "banks (margins < 1)\n\n")

stan_opts <- list(chains = 4, iter = 2000, warmup = 1000, seed = 1960,
                  cores = 1, control = list(adapt_delta = 0.95, max_treedepth = 15))

# --- Logit Bertrand ---
cat("=== Logit + Bertrand ===\n")
fit_logit <- bayes_est(
  prices       = pnb_sub$prices,
  quantities   = pnb_sub$quantities,
  margins      = pnb_sub$margins,
  market       = rep("PNB", nrow(pnb_sub)),
  group        = pnb_sub$group,
  demand_model = "logit",
  supply_model = "bertrand",
  use_rho      = FALSE,
  use_cutoff   = FALSE,
  stan_options = stan_opts
)
print(fit_logit)

# --- CES Bertrand ---
cat("\n=== CES + Bertrand ===\n")
fit_ces <- bayes_est(
  prices       = pnb_sub$prices,
  quantities   = pnb_sub$quantities,
  margins      = pnb_sub$margins,
  market       = rep("PNB", nrow(pnb_sub)),
  group        = pnb_sub$group,
  demand_model = "ces",
  supply_model = "bertrand",
  use_rho      = FALSE,
  use_cutoff   = FALSE,
  stan_options = stan_opts
)
print(fit_ces)

# --- WTP conversion ---
cat("\n=== WTP Conversion ===\n")
gamma_from_logit <- wtp_convert(fit_logit, to_gamma = TRUE)
alpha_from_ces   <- wtp_convert(fit_ces, to_gamma = FALSE)
cat("Logit alpha -> gamma equiv: mean =", round(mean(gamma_from_logit), 3),
    " sd =", round(sd(gamma_from_logit), 3), "\n")
cat("CES gamma -> alpha equiv:   mean =", round(mean(alpha_from_ces), 3),
    " sd =", round(sd(alpha_from_ces), 3), "\n")

# --- Model comparison ---
cat("\n=== Model Comparison ===\n")
comp <- bayes_compare(Logit = fit_logit, CES = fit_ces, bridge = FALSE)
print(comp)
