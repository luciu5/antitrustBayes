#!/usr/bin/env Rscript
# Independent small-N algebra reference; no Stan target or economic changes.
options(stringsAsFactors = FALSE)
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
out <- file.path(repo, "doc/statistical_toolkit_runs/20260916/performance-audit/recovery",
                 "market-quality/price-coordinate-reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
p <- c(1.7, 2.3, 3.4, 4.9)
p0 <- 1.3
b <- c(-.3, .2, -.1, .2)
B <- .7
eta <- qlogis(.24)
checks <- list()
check <- function(name, error, tolerance = 2e-12) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = is.finite(error) && error <= tolerance)
}
for (beta in c(-6, .4, 1.7)) {
  standard <- eta + B + b + beta * (p - mean(p))
  physical_delta <- B + b - beta * (mean(p) - p0)
  physical <- eta + physical_delta + beta * (p - p0)
  check(paste0("Logit_raw_centered_map_beta_", beta), max(abs(standard - physical)))
  for (factor in c(.01, 2, 100)) {
    scaled <- eta + B + b + (beta / factor) * (factor * p - mean(factor * p))
    check(paste0("Logit_common_price_units_", beta, "_", factor), max(abs(standard - scaled)))
  }
}
unit_rows <- list()
for (gamma in c(1.1, 3.5, 7)) {
  standard <- eta + B + b - gamma * (log(p) - mean(log(p))) + gamma * log(p0)
  log_A <- B + b + gamma * mean(log(p))
  physical <- eta + log_A + gamma * log(p0 / p)
  check(paste0("CES_raw_centered_map_gamma_", gamma), max(abs(standard - physical)))
  for (factor in c(.01, 2, 100)) {
    scaled <- eta + B + b - gamma * (log(factor * p) - mean(log(factor * p))) +
      gamma * log(factor * p0)
    expected_shift <- gamma * log(factor)
    check(paste0("CES_price_units_shift_", gamma, "_", factor),
          max(abs(scaled - standard - expected_shift)))
    remapped <- scaled - expected_shift
    check(paste0("CES_intercept_remapping_", gamma, "_", factor), max(abs(remapped - standard)))
    unit_rows[[length(unit_rows) + 1L]] <- data.frame(gamma = gamma, price_factor = factor,
      predictor_shift = expected_shift,
      physical_intercept_required_shift = -expected_shift,
      interpretation = "Equation identity; fixed single-market intercept cannot implement this remapping")
  }
}
result <- do.call(rbind, checks)
write.csv(result, file.path(out, "validation.csv"), row.names = FALSE)
write.csv(do.call(rbind, unit_rows), file.path(out, "CES-unit-shifts.csv"), row.names = FALSE)
print(result)
stopifnot(all(result$passed))
