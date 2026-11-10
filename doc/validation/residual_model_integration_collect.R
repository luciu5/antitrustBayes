# Collect completed artifacts without fitting, changing seals, or hiding failures.
rmi_collect <- function(output, report_dir = "doc/validation") {
  manifest <- readRDS(file.path(output, "manifest.rds"))
  prefix <- file.path(report_dir, "residual_model_integration")
  write <- function(x, suffix) write.csv(x, paste0(prefix, "_", suffix, ".csv"), row.names = FALSE)
  tables <- setNames(vector("list", 7L), c("diagnostics", "parameters", "shares", "firm_shares",
    "margins", "policies", "status"))
  for (id in names(manifest$calls)) {
    directory <- file.path(output, id); attempt <- manifest$attempts[[id]]
    tables$status[[id]] <- data.frame(id, status = if (is.null(attempt)) "not_attempted" else attempt$status,
      fit_checkpoint = file.exists(file.path(directory, "fit.rds")),
      fit_seal = file.exists(file.path(directory, "seal.rds")),
      policy_report = file.exists(file.path(directory, "policy_all_rows.csv")),
      exit_status = if (is.null(attempt$exit_status)) NA_integer_ else attempt$exit_status)
    files <- c(diagnostics = "diagnostics.csv", parameters = "parameters.csv",
      shares = "share_prediction_means.csv", firm_shares = "firm_share_prediction_means.csv",
      margins = "margin_fit.csv")
    for (type in names(files)) {
      file <- file.path(directory, files[[type]])
      if (file.exists(file)) tables[[type]][[id]] <- cbind(id = id, read.csv(file))
    }
    file <- file.path(directory, "policy_all_rows.csv")
    if (file.exists(file)) {
      x <- read.csv(file); valid <- x$status == "valid"
      maximum <- function(name) if (any(valid)) max(x[[name]][valid]) else NA_real_
      minimum <- function(name) if (any(valid)) min(x[[name]][valid]) else NA_real_
      tables$policies[[id]] <- data.frame(id, attempted = nrow(x), valid = sum(valid), invalid = sum(!valid),
        max_pre_foc = maximum("pre_foc"), max_post_foc = maximum("post_foc"),
        max_no_op_error = maximum("no_op_error"),
        mean_market_paasche_pp = if (any(valid)) 100 * mean(x$market_paasche[valid]) else NA_real_,
        min_market_paasche_pp = 100 * minimum("market_paasche"),
        max_market_paasche_pp = 100 * maximum("market_paasche"),
        mean_party_paasche_pp = if (any(valid)) 100 * mean(x$party_paasche[valid]) else NA_real_)
    }
  }
  for (type in names(tables)) if (length(tables[[type]])) {
    tables[[type]] <- do.call(rbind, tables[[type]]); write(tables[[type]], type)
  }
  native <- file.path(output, "native/checks.csv")
  if (file.exists(native)) write(read.csv(native), "native_checks")
  write(read.csv(file.path(output, "fixture_validation.csv")), "fixture_checks")
  if (!is.null(manifest$source_hashes)) write(data.frame(path = names(manifest$source_hashes),
    sha256 = unname(manifest$source_hashes)), "source_seal")
  write(data.frame(case = names(manifest$fixture_data_hashes),
    data_sha256 = unname(manifest$fixture_data_hashes),
    truth_sha256 = unname(manifest$fixture_truth_hashes)), "fixture_seal")
  invisible(tables)
}
