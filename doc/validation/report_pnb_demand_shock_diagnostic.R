# Read-only inference collection: no sampling, compilation, bridges or averaging.
args <- commandArgs(TRUE)
output <- normalizePath(if (length(args)) args[[1L]] else ".scratch/demand_shock_diagnostic_20261109", mustWork = TRUE)
library(digest)
library(antitrustBayes)
manifest <- readRDS(file.path(output, "manifest.rds"))
inputs <- readRDS(file.path(output, "inputs.rds"))
cells <- c("snapshot_reference", "snapshot_candidate", "history_reference", "history_candidate")
tables <- parameters <- margins <- list()
status <- data.frame(Cell = cells, Status = "unavailable", Chains = NA_integer_,
  Retained_Draws = NA_integer_, Max_Rhat = NA_real_, Min_Bulk_ESS = NA_real_,
  Min_Tail_ESS = NA_real_, Divergences = NA_integer_, Treedepth_Hits = NA_integer_,
  Policy_Requested = 100L, Baseline_Valid = NA_integer_, Policy_Valid = NA_integer_,
  Max_FOC = NA_real_, Max_Noop = NA_real_, stringsAsFactors = FALSE)
summarize <- function(x) c(Mean = mean(x), SD = sd(x),
  Lower = quantile(x, .025, names = FALSE), Upper = quantile(x, .975, names = FALSE))
for (j in seq_along(cells)) {
  cell <- cells[j]; information <- sub("_.*", "", cell); candidate <- grepl("candidate$", cell)
  directory <- file.path(output, cell)
  if (!file.exists(file.path(directory, "policy_draws.rds"))) next
  policy <- readRDS(file.path(directory, "policy_draws.rds"))
  d <- read.csv(file.path(directory, "sampler_diagnostics.csv"))
  status[j, names(d)] <- d[1, ]
  stopifnot(length(policy$indices) == 100L, nrow(policy$effects) == 100L,
    !anyDuplicated(policy$indices))
  status$Baseline_Valid[j] <- sum(apply(is.finite(policy$share_solved), 1L, all))
  status$Policy_Valid[j] <- sum(is.na(policy$effects$Error))
  status$Max_FOC[j] <- max(policy$effects$FOC, na.rm = TRUE)
  status$Max_Noop[j] <- max(policy$effects$Noop, na.rm = TRUE)
  status$Status[j] <- if (j == 1L) "reused_sealed_reference_2x250" else "completed_new_4x500"
  x <- read.csv(file.path(directory, "comparison.csv"))
  x$Information <- information; x$Arm <- if (candidate) "candidate" else "reference"
  x$Chains <- d$Chains; x$Retained_Draws <- d$Retained_Draws
  x$Max_Rhat <- d$Max_Rhat; x$Min_Bulk_ESS <- d$Min_Bulk_ESS
  x$Convergence <- if (d$Max_Rhat <= 1.01 && d$Min_Bulk_ESS >= 400 &&
    d$Divergences == 0 && d$Treedepth_Hits == 0) "adequate_for_bounded_diagnostic" else "inconclusive"
  tables[[cell]] <- x
  if (j == 1L) {
    stopifnot(identical(digest(file = manifest$saved_reference_path, algo = "sha256"), manifest$saved_reference_sha256))
    native <- readRDS(manifest$saved_reference_path)@stanfit
  } else {
    saved <- readRDS(file.path(directory, "fit.rds"))
    input <- inputs[[information]]
    stopifnot(identical(saved$data_sha256, digest(input$data, algo = "sha256")),
      identical(saved$prototype, candidate), identical(saved$config, manifest$config))
    kernel <- file.path(output, "build", if (candidate) "period_demand_shock.stan" else "reference_logit_supply.stan")
    stopifnot(identical(saved$kernel_sha256, digest(file = kernel, algo = "sha256")))
    native <- saved$native_fit
  }
  variables <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho", "pred_marginNormInv")
  if (candidate) variables <- c(variables, "sigma_demand_shock", "rho_demand_shock")
  draws <- antitrustBayes:::.bayes_backend_extract(native, variables = variables)
  fields <- list(Alpha_Per_Rate_Unit = as.numeric(draws$a_market[, 1]) / inputs[[information]]$data$price_sd,
    Persistent_Quality_SD = as.numeric(draws$sigma_b_strat_economic),
    Total_Share_Residual_SD = as.numeric(draws$sigma_share),
    Margin_Residual_SD = as.numeric(draws$sigma_margin), Rho_Effective = as.numeric(draws$rho),
    Outside_Share = as.numeric(draws$s0_prob[, 1]))
  if (candidate) fields <- c(fields, list(Period_Shock_SD = as.numeric(draws$sigma_demand_shock),
    Rho_Shock = as.numeric(draws$rho_demand_shock), Fixed_Measurement_SD = rep(.1, nrow(as.matrix(draws$a_market)))))
  parameters[[cell]] <- do.call(rbind, lapply(names(fields), function(name)
    data.frame(Cell = cell, Parameter = name, t(summarize(fields[[name]])))))
  target <- which(inputs[[information]]$rows$year == "1960")
  observed <- as.logical(inputs[[information]]$data$margin_obs[target])
  idx <- target[observed]
  predicted <- log(draws$pred_marginNormInv[policy$indices, idx, drop = FALSE] + 1e-9)
  error <- sweep(predicted, 2L, inputs[[information]]$data$log_marginNormInv[idx], "-")
  sigma <- as.numeric(draws$sigma_margin[policy$indices])
  margins[[cell]] <- data.frame(Cell = cell, N_Target_Margins = length(idx),
    Mean_LogInv_Error = mean(error), Point_RMSE = sqrt(mean(error^2)),
    Marginal_Predictive_RMSE = sqrt(mean(error^2 + sigma^2)))
}
write.csv(do.call(rbind, tables), file.path(output, "comparison_completed.csv"), row.names = FALSE)
write.csv(status, file.path(output, "execution_status_completed.csv"), row.names = FALSE)
write.csv(do.call(rbind, parameters), file.path(output, "parameter_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, margins), file.path(output, "margin_fit.csv"), row.names = FALSE)
files <- c(manifest$saved_reference_path, file.path(output, "clock.json"),
  file.path(output, "manifest.rds"), file.path(output, "build", c("reference_logit_supply.stan", "period_demand_shock.stan")),
  unlist(lapply(cells, function(cell) list.files(file.path(output, cell),
    pattern = "[.](rds|csv)$", full.names = TRUE))))
write.csv(data.frame(Path = files, SHA256 = vapply(files, function(p)
  digest(file = p, algo = "sha256"), character(1))), file.path(output, "completed_artifact_seals.csv"), row.names = FALSE)
cat("Collected all available cells without inference. Original first-attempt records remain unchanged.\n")
print(status)
