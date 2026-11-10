# Saved-output collection only; never compiles, samples or fills missing cells.
args <- commandArgs(TRUE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
output <- normalizePath(if (length(args)) args[[1L]] else ".scratch/asset_variance_diagnostic_20261109", mustWork = TRUE)
manifest <- readRDS(file.path(output, "manifest.rds"))
plan <- read.csv(file.path(output, "execution_status.csv"), stringsAsFactors = FALSE)
inputs <- readRDS(file.path(output, "inputs.rds"))
hash <- function(path) digest::digest(file = path, algo = "sha256")
stopifnot(all(vapply(names(manifest$source_hash), function(path)
  identical(hash(path), manifest$source_hash[[path]]), logical(1))))
metrics <- parameters <- validity <- list()
for (cell in plan$Cell) {
  folder <- file.path(output, cell)
  if (!file.exists(file.path(folder, "policy_draws.rds"))) next
  input <- inputs[[sub("_.*", "", cell)]]
  policy <- readRDS(file.path(folder, "policy_draws.rds"))
  target <- which(input$rows$year == manifest$config$target_year)
  obs <- input$rows[target, ]
  truth <- obs$quantities / sum(obs$quantities)
  native <- if (cell == "history_reference" && !file.exists(file.path(folder, "fit.rds")))
    readRDS(manifest$reused_history_path)$native_fit else readRDS(file.path(folder, "fit.rds"))$native_fit
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
}
if (length(metrics)) write.csv(do.call(rbind, metrics), file.path(output, "fit_metrics.csv"), row.names = FALSE)
if (length(parameters)) write.csv(do.call(rbind, parameters), file.path(output, "parameter_comparison.csv"), row.names = FALSE)
if (length(validity)) write.csv(do.call(rbind, validity), file.path(output, "validity_summary.csv"), row.names = FALSE)
files <- list.files(output, recursive = TRUE, full.names = TRUE)
files <- files[!file.info(files)$isdir & grepl("\\.(rds|csv|stan|json|R)$", files) &
  !grepl("artifact_seals.csv$", files)]
seals <- data.frame(Path = files, SHA256 = vapply(files, hash, character(1)))
write.csv(seals, file.path(output, "artifact_seals.csv"), row.names = FALSE)
print(plan)
if (length(metrics)) print(do.call(rbind, metrics))
if (length(validity)) print(do.call(rbind, validity))
