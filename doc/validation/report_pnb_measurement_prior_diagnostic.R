# Collection and independent checks only; never fits, compiles or fills a cell.
args <- commandArgs(TRUE)
root <- normalizePath(".", mustWork = TRUE)
output <- normalizePath(if (length(args)) args[[1L]] else ".scratch/measurement_prior_diagnostic_20261109", mustWork = TRUE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
.libPaths(c(file.path(root, "../mergerBayes/.scratch/paper-library"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
pkgload::load_all(root, quiet = TRUE)
source("doc/validation/demand_shock_prototype.R")
source("doc/validation/measurement_prior_prototype.R")
manifest <- readRDS(file.path(output, "manifest.rds"))
inputs <- readRDS(file.path(output, "inputs.rds"))
plan <- read.csv(file.path(output, "execution_status.csv"))
hash <- function(path) digest::digest(file = path, algo = "sha256")
verify_source <- function(path) {
  expected <- manifest$source_hash[[path]]
  if (identical(hash(path), expected)) return(TRUE)
  # The initial summary had one complemented prior-tail probability reversed.
  # Accept only that exact reporting edit, with the executed source still
  # matching its original seal. No kernel, data, prior or inference exception.
  archived <- file.path(output, "reporting_correction/measurement_prior_prototype.executed.R")
  if (path != file.path(root, "doc/validation/measurement_prior_prototype.R") ||
      !file.exists(archived) || !identical(hash(archived), expected)) return(FALSE)
  old <- readLines(archived, warn = FALSE)
  repaired <- sub("Pr_Above_20 = plnorm(.20, mu, sdlog)",
    "Pr_Above_20 = 1 - plnorm(.20, mu, sdlog)", old, fixed = TRUE)
  identical(repaired, readLines(path, warn = FALSE))
}
stopifnot(all(vapply(names(manifest$source_hash), verify_source, logical(1))))
corrected <- file.exists(file.path(output, "reporting_correction/measurement_prior_prototype.executed.R"))
suffix <- if (corrected) "_corrected" else ""
metrics <- validity <- likelihood <- prior_posterior <- audit <- list()
finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
for (cell in plan$Cell) {
  folder <- file.path(output, cell)
  if (!file.exists(file.path(folder, "policy_draws.rds"))) next
  information <- sub("_.*", "", cell)
  input <- inputs[[information]]; data <- input$data
  policy <- readRDS(file.path(folder, "policy_draws.rds"))
  reference_cell <- if (grepl("observed$", cell)) paste0(information, "_reference") else cell
  if (reference_cell %in% names(manifest$reused_reference_paths)) {
    pointer <- manifest$reused_references[[reference_cell]]
    stopifnot(identical(hash(pointer$fit_path), pointer$fit_sha256))
    saved <- readRDS(pointer$fit_path)
  } else {
    saved <- readRDS(file.path(folder, "fit.rds"))
    stopifnot(identical(saved$config, manifest$config), identical(saved$data, data),
      identical(saved$data_sha256, digest::digest(data, algo = "sha256")),
      identical(saved$kernel_sha256, hash(file.path(output, "build/estimated_measurement_logit_supply.stan"))))
  }
  native <- saved$native_fit
  estimated <- grepl("estimated$", cell); shock <- estimated || grepl("fixed$", cell)
  variables <- c("relative_share_predictor", "pred_marginNormInv", "sigma_margin", "sigma_share",
    "rho", "sigma_b_strat_economic", "relative_share_loglik", "relative_joint_loglik")
  if (shock) variables <- c(variables, "sigma_demand_shock", "rho_demand_shock", "diagnostic_measurement_sd")
  d <- ds_private(".bayes_backend_extract")(native, variables = variables)
  target <- which(input$rows$year == manifest$config$target_year)
  truth <- input$rows$quantities[target] / sum(input$rows$quantities[target])
  idx <- policy$indices
  pred <- log(d$pred_marginNormInv[idx, target, drop = FALSE] + 1e-9)
  errors <- sweep(pred, 2L, data$log_marginNormInv[target], "-")
  paired <- as.logical(data$margin_obs[target])
  metrics[[cell]] <- data.frame(Cell = cell,
    Margin_Point_RMSE = sqrt(mean(errors[, paired, drop = FALSE]^2)),
    Margin_Predictive_RMSE = sqrt(mean(errors[, paired, drop = FALSE]^2 + as.numeric(d$sigma_margin[idx])^2)),
    Share_Solved_RMSE_PP = sqrt(mean(sweep(100 * policy$share_solved, 2L, 100 * truth, "-")^2)))
  likelihood[[cell]] <- data.frame(Cell = cell,
    Mean_Share_Log_Likelihood = mean(d$relative_share_loglik), Mean_Joint_Log_Likelihood = mean(d$relative_joint_loglik))
  effects <- policy$effects
  validity[[cell]] <- data.frame(Cell = cell, Requested = nrow(effects),
    Structural_Valid = sum(is.na(effects$Error) & is.finite(effects$FOC)),
    Recovery_or_Equilibrium_Failed = sum(!is.na(effects$Error)),
    Market_Index_Valid = sum(is.finite(effects$Market_Paasche)), Party_Index_Valid = sum(is.finite(effects$Party_Paasche)),
    Market_Index_Failed = sum(!is.na(effects$Market_Index_Error)), Party_Index_Failed = sum(!is.na(effects$Party_Index_Error)),
    Max_FOC = finite_max(effects$FOC), Max_Noop = finite_max(effects$Noop))
  stopifnot(length(idx) == manifest$config$policy_draws, identical(as.integer(idx), as.integer(effects$Posterior_Row)))
  max_share_error <- NA_real_; share_audited <- 0L
  for (j in seq_along(idx)) if (is.null(policy$states[[j]]$recovery_error)) {
    if (estimated) stopifnot(identical(policy$states[[j]]$measurement_sd,
      as.numeric(d$diagnostic_measurement_sd[idx[j]])))
    if (!all(is.finite(policy$share_at_observed[j, ]))) next
    score <- d$relative_share_predictor[idx[j], target] + policy$states[[j]]$structural_state[target]
    s <- exp(score - max(score)); s <- s / sum(s)
    max_share_error <- finite_max(c(max_share_error, abs(s - policy$share_at_observed[j, ])))
    share_audited <- share_audited + 1L
  }
  if (share_audited) stopifnot(max_share_error < 1e-12)
  if (grepl("observed$", cell)) {
    valid <- apply(is.finite(policy$share_solved), 1L, all)
    if (any(valid)) stopifnot(max(abs(sweep(policy$share_solved[valid, , drop = FALSE], 2, truth, "-"))) < 1e-10)
  }
  # Independent dense integration of persistent quality in Helmert contrasts.
  # This checks normalized native likelihoods, not the implementation's helper.
  n <- data$N; groups <- split(seq_len(n), data$market_year_idx)
  basis <- matrix(0, n, sum(lengths(groups) - 1L)); offset <- 0L
  for (ii in groups) if (length(ii) > 1L) {
    h <- qr.Q(qr(contr.helmert(length(ii))))
    basis[ii, offset + seq_len(ncol(h))] <- h; offset <- offset + ncol(h)
  }
  g <- data$N_quality; centered <- matrix(0, g, g)
  for (ii in split(seq_len(g), rep(seq_along(data$quality_block_n), data$quality_block_n)))
    centered[ii, ii] <- diag(length(ii)) - matrix(1 / length(ii), length(ii), length(ii))
  z <- diag(g)[as.integer(data$obs_to_quality), , drop = FALSE]
  covariance_quality_unit <- crossprod(basis, z %*% centered %*% t(z) %*% basis)
  gaussian <- function(residual, covariance) {
    r <- as.numeric(crossprod(basis, residual)); chol_c <- chol(covariance)
    -.5 * (length(r) * log(2 * pi) + 2 * sum(log(diag(chol_c))) +
      sum(backsolve(chol_c, r, transpose = TRUE)^2))
  }
  max_ll_error <- 0
  for (i in idx[seq_len(min(8L, length(idx)))]) {
    s <- as.numeric(d$sigma_share[i]); rho <- as.numeric(d$rho[i]); sm <- as.numeric(d$sigma_margin[i])
    tau <- as.numeric(d$sigma_b_strat_economic[i]); paired_all <- as.logical(data$margin_obs)
    residual <- as.numeric(data$log_shareIn) - d$relative_share_predictor[i, ]
    margin_residual <- as.numeric(data$log_marginNormInv) - log(d$pred_marginNormInv[i, ] + 1e-9)
    variance <- rep(s^2, n); conditional <- residual
    conditional[paired_all] <- conditional[paired_all] - rho * s / sm * margin_residual[paired_all]
    variance[paired_all] <- variance[paired_all] * (1 - rho^2)
    conditional_covariance <- tau^2 * covariance_quality_unit + crossprod(basis, variance * basis)
    joint <- gaussian(conditional, conditional_covariance) + sum(dnorm(margin_residual[paired_all], 0, sm, log = TRUE))
    share <- gaussian(residual, tau^2 * covariance_quality_unit + diag(s^2, ncol(basis)))
    max_ll_error <- max(max_ll_error, abs(joint - d$relative_joint_loglik[i]), abs(share - d$relative_share_loglik[i]))
    if (shock) {
      omega <- as.numeric(d$sigma_demand_shock[i]); m <- as.numeric(d$diagnostic_measurement_sd[i]); rd <- as.numeric(d$rho_demand_shock[i])
      stopifnot(abs(s^2 - omega^2 - m^2) < 1e-10, abs(rho * s - rd * omega) < 1e-10)
    }
  }
  # Older reference files used 8-digit CSV precision; retain their rounding bound.
  stopifnot(max_ll_error < if (grepl("reference$|observed$", cell)) 1e-4 else 1e-8)
  replay_exact <- NA
  if (estimated) {
    replay_folder <- file.path(folder, "replay_check"); dir.create(replay_folder, showWarnings = FALSE)
    if (!file.exists(file.path(replay_folder, "policy_draws.rds")))
      mp_make_evaluator()(native, input, TRUE, replay_folder,
        n_draws = manifest$config$policy_draws, seed = manifest$config$policy_seed)
    replay <- readRDS(file.path(replay_folder, "policy_draws.rds"))
    replay_exact <- identical(policy, replay); stopifnot(replay_exact)
    p <- read.csv(file.path(folder, "measurement_prior_posterior.csv")); p$Cell <- cell
    p$Pr_Above_20[p$Distribution == "prior"] <- plnorm(.20,
      manifest$config$measurement_prior_log_mean, manifest$config$measurement_prior_log_sd,
      lower.tail = FALSE)
    prior_posterior[[cell]] <- p
  }
  audit[[cell]] <- data.frame(Cell = cell, Independent_Likelihood_Rows = min(8L, length(idx)),
    Max_Normalized_Likelihood_Error = max_ll_error, Max_Independent_Share_Error = max_share_error,
    Share_Rows_Audited = share_audited, Same_Plan_Exact_Replay = replay_exact)
}
for (name in c("metrics", "validity", "likelihood", "prior_posterior", "audit")) {
  x <- get(name)
  if (length(x)) write.csv(do.call(rbind, x), file.path(output, paste0(name, "_summary", suffix, ".csv")), row.names = FALSE)
}
files <- list.files(output, recursive = TRUE, full.names = TRUE)
files <- files[!file.info(files)$isdir & grepl("\\.(rds|csv|stan|json|R)$", files) & !grepl("artifact_seals(_corrected)?[.]csv$", files)]
write.csv(data.frame(Path = files, SHA256 = vapply(files, hash, character(1))),
  file.path(output, paste0("artifact_seals", suffix, ".csv")), row.names = FALSE)
print(plan)
if (length(audit)) print(do.call(rbind, audit))
if (length(prior_posterior)) print(do.call(rbind, prior_posterior))
