# Saved-output collection only; never compiles, samples or fills missing cells.
args <- commandArgs(TRUE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
output <- normalizePath(if (length(args)) args[[1L]] else ".scratch/variance_time_trend_diagnostic_20261109", mustWork = TRUE)
manifest <- readRDS(file.path(output, "manifest.rds"))
plan <- read.csv(file.path(output, "execution_status.csv"), stringsAsFactors = FALSE)
inputs <- readRDS(file.path(output, "inputs.rds"))
hash <- function(path) digest::digest(file = path, algo = "sha256")
stopifnot(all(vapply(names(manifest$source_hash), function(path)
  identical(hash(path), manifest$source_hash[[path]]), logical(1))))
metrics <- parameters <- validity <- likelihood <- dispersion <- list()
for (cell in plan$Cell) {
  folder <- file.path(output, cell)
  if (!file.exists(file.path(folder, "policy_draws.rds"))) next
  input <- inputs[[sub("_.*", "", cell)]]
  policy <- readRDS(file.path(folder, "policy_draws.rds"))
  target <- which(input$rows$year == manifest$config$target_year)
  obs <- input$rows[target, ]
  truth <- obs$quantities / sum(obs$quantities)
  native <- if (cell %in% names(manifest$reused_reference_paths))
    readRDS(manifest$reused_reference_paths[[cell]])$native_fit else readRDS(file.path(folder, "fit.rds"))$native_fit
  ll <- posterior::as_draws_matrix(native$draws(variables = c("relative_share_loglik", "relative_joint_loglik")))
  likelihood[[cell]] <- data.frame(Cell = cell,
    Mean_Collapsed_Share_Log_Likelihood = mean(ll[, "relative_share_loglik"]),
    Mean_Collapsed_Joint_Log_Likelihood = mean(ll[, "relative_joint_loglik"]))
  draw <- posterior::as_draws_matrix(native$draws(variables = c("pred_marginNormInv", "sigma_margin")))
  pred_columns <- paste0("pred_marginNormInv[", target, "]")
  prediction <- log(draw[policy$indices, pred_columns, drop = FALSE] + 1e-9)
  sigma_m <- as.numeric(draw[policy$indices, "sigma_margin"])
  err_margin <- sweep(prediction, 2L, input$data$log_marginNormInv[target], "-")
  paired <- as.logical(input$data$margin_obs[target])
  err_obs <- sweep(100 * policy$share_at_observed, 2L, 100 * truth, "-")
  err_solved <- sweep(100 * policy$share_solved, 2L, 100 * truth, "-")
  metrics[[cell]] <- data.frame(Cell = cell,
    Margin_Point_RMSE = sqrt(mean(err_margin[, paired, drop = FALSE]^2)),
    Margin_Predictive_RMSE = sqrt(mean(err_margin[, paired, drop = FALSE]^2 + sigma_m^2)),
    All_Bank_At_Observed_RMSE_PP = sqrt(mean(err_obs^2)),
    All_Bank_Solved_RMSE_PP = sqrt(mean(err_solved^2)),
    All_Bank_Solved_MAE_PP = mean(abs(err_solved)),
    Max_Independent_FOC = max(policy$effects$FOC, na.rm = TRUE),
    Max_Noop = max(policy$effects$Noop, na.rm = TRUE))
  p <- read.csv(file.path(folder, "parameter_summary.csv")); p$Cell <- cell
  parameters[[cell]] <- p
  e <- policy$effects
  validity[[cell]] <- data.frame(Cell = cell, Requested = nrow(e),
    Structural_Valid = sum(is.na(e$Error) & is.finite(e$FOC)),
    Recovery_or_Equilibrium_Failed = sum(!is.na(e$Error)),
    Whole_Market_Valid = sum(is.finite(e$Market_Paasche)),
    Party_Valid = sum(is.finite(e$Party_Paasche)),
    Whole_Market_Index_Failed = sum(!is.na(e$Market_Index_Error)),
    Party_Index_Failed = sum(!is.na(e$Party_Index_Error)))
  stopifnot(identical(as.integer(policy$indices), as.integer(e$Posterior_Row)),
    nrow(e) == manifest$config$policy_draws)
  # Replay the state's demand shares directly, independent of the summary CSV.
  score <- posterior::as_draws_matrix(native$draws(variables = "relative_share_predictor"))
  for (j in seq_along(policy$indices)) if (is.null(policy$states[[j]]$recovery_error)) {
    predictor <- score[policy$indices[j], paste0("relative_share_predictor[", target, "]")]
    v <- predictor + policy$states[[j]]$structural_state[target]
    s <- exp(v - max(v)); s <- s / sum(s)
    stopifnot(max(abs(s - policy$share_at_observed[j, ])) < 1e-12)
  }
  scale_variables <- c("sigma_share", "sigma_margin", "rho", "pred_marginNormInv")
  heteroskedastic <- cell != "history_reference"
  if (heteroskedastic) scale_variables <- c(scale_variables, "sigma_share_observation")
  scale_draws <- posterior::as_draws_matrix(native$draws(variables = scale_variables))
  data <- input$data; n <- data$N; paired <- as.logical(data$margin_obs)
  by_year <- list()
  for (j in seq_along(policy$indices)) {
    if (!is.null(policy$states[[j]]$recovery_error)) next
    k <- policy$indices[j]
    predictor <- as.numeric(score[k, paste0("relative_share_predictor[", seq_len(n), "]")])
    sigma <- if (heteroskedastic)
      as.numeric(scale_draws[k, paste0("sigma_share_observation[", seq_len(n), "]")]) else
      rep(as.numeric(scale_draws[k, "sigma_share"]), n)
    rho <- as.numeric(scale_draws[k, "rho"]); sigma_m <- as.numeric(scale_draws[k, "sigma_margin"])
    pred_margin <- as.numeric(scale_draws[k, paste0("pred_marginNormInv[", seq_len(n), "]")])
    residual <- as.numeric(data$log_shareIn) - predictor - policy$states[[j]]$structural_state
    shift <- numeric(n)
    shift[paired] <- rho * sigma[paired] / sigma_m *
      (data$log_marginNormInv[paired] - log(pred_margin[paired] + 1e-9))
    variance <- sigma^2 * (1 - rho^2 * paired)
    for (year in sort(unique(input$rows$year))) {
      ii <- which(input$rows$year == year)
      basis <- qr.Q(qr(contr.helmert(length(ii))))
      raw <- as.numeric(crossprod(basis, residual[ii]))
      adjusted <- as.numeric(crossprod(basis, residual[ii] - shift[ii]))
      covariance <- crossprod(basis, variance[ii] * basis)
      standardized <- sqrt(as.numeric(crossprod(adjusted, solve(covariance, adjusted))) / length(raw))
      by_year[[length(by_year) + 1L]] <- data.frame(Cell = cell, Year = year,
        Draw = j, N = length(ii), Raw_RMS = sqrt(mean(raw^2)), Standardized_RMS = standardized)
    }
  }
  y <- do.call(rbind, by_year)
  dispersion[[cell]] <- do.call(rbind, lapply(split(y, y$Year), function(z) data.frame(
    Cell = cell, Year = z$Year[1L], N = z$N[1L], Raw_RMS_Mean = mean(z$Raw_RMS),
    Standardized_RMS_Mean = mean(z$Standardized_RMS),
    Standardized_RMS_Lower = unname(quantile(z$Standardized_RMS, .025)),
    Standardized_RMS_Upper = unname(quantile(z$Standardized_RMS, .975)))))
}
if (length(metrics)) write.csv(do.call(rbind, metrics), file.path(output, "fit_metrics.csv"), row.names = FALSE)
if (length(parameters)) write.csv(do.call(rbind, parameters), file.path(output, "parameter_comparison.csv"), row.names = FALSE)
if (length(validity)) write.csv(do.call(rbind, validity), file.path(output, "validity_summary.csv"), row.names = FALSE)
if (length(likelihood)) write.csv(do.call(rbind, likelihood), file.path(output, "normalized_likelihood_diagnostic.csv"), row.names = FALSE)
if (length(dispersion)) write.csv(do.call(rbind, dispersion), file.path(output, "residual_dispersion_by_year.csv"), row.names = FALSE)
files <- list.files(output, recursive = TRUE, full.names = TRUE)
files <- files[!file.info(files)$isdir & grepl("\\.(rds|csv|stan|json|R)$", files) &
  !grepl("artifact_seals.csv$", files)]
seals <- data.frame(Path = files, SHA256 = vapply(files, hash, character(1)))
write.csv(seals, file.path(output, "artifact_seals.csv"), row.names = FALSE)
print(plan)
if (length(metrics)) print(do.call(rbind, metrics))
if (length(validity)) print(do.call(rbind, validity))
