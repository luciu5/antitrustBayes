# Separate asset-heteroskedasticity experiment; no production fit is migrated.
# Source demand_shock_prototype.R first for its unchanged input/policy helpers.
av_kernel <- function(source, destination) {
  text <- paste(readLines(source, warn = FALSE), collapse = "\n")
  once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  exactly_twice <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 2L)
    text <<- gsub(old, new, text, fixed = TRUE)
  }
  once("  vector[N] shareIn;", paste(
    "  int<lower=0, upper=1> asset_variance_active;",
    "  vector[N] institution_log_assets_z;",
    "  real<lower=0> asset_gamma_prior_sd;",
    "  vector[N] shareIn;", sep = "\n"))
  once("  real log_sigma_share;", paste(
    "  vector[asset_variance_active] asset_gamma;",
    "  real log_sigma_share;", sep = "\n"))
  once("  vector[N_market] sigma_share_market;", paste(
    "  vector[N] sigma_share_observation;",
    "  vector[N_market] sigma_share_market;", sep = "\n"))
  once("  if (rho_active == 1) rho = 0.95 * tanh(rho_raw[1]);", paste(
    "  for (i in 1:N) sigma_share_observation[i] = sigma_share_market[market[i]]",
    "    * exp(asset_variance_active == 1 ? asset_gamma[1] * institution_log_assets_z[i] : 0.0);",
    "  if (rho_active == 1) rho = 0.95 * tanh(rho_raw[1]);", sep = "\n"))
  once("  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", paste(
    "  target += normal_lpdf(asset_gamma | 0, asset_gamma_prior_sd);",
    "  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", sep = "\n"))
  exactly_twice("for (i in 1:N) relative_base_var[i] = square(sigma_share_market[market[i]]);",
    "for (i in 1:N) relative_base_var[i] = square(sigma_share_observation[i]);")
  exactly_twice("relative_cond_resid[i] -= rho * sigma_share_market[m]",
    "relative_cond_resid[i] -= rho * sigma_share_observation[i]")
  # This bounded prototype implements exactly the audited Gaussian ILR route.
  once("  rho_active = (use_rho == 1 && likelihood_mode == 1 && any_paired == 1) ? 1 : 0;",
    paste("  if (share_basis != 1 || quality_parameterization != 2 || N_market != 1",
      "      || share_likelihood != 2 || margin_likelihood != 2 || use_market_resid_scale != 0",
      "      || supply_model != 1 || use_cutoff != 0)",
      "    reject(\"Asset diagnostic requires one-market collapsed Gaussian relative-share Bertrand.\");",
      "  rho_active = (use_rho == 1 && likelihood_mode == 1 && any_paired == 1) ? 1 : 0;", sep = "\n"))
  writeLines(c("// Separate diagnostic: log sigma_i = log sigma + gamma z(log institution assets).", text), destination)
  invisible(destination)
}

av_asset_index <- function(rows) {
  required <- c("tophold", "year", "total_assets")
  stopifnot(all(required %in% names(rows)), all(is.finite(rows$total_assets)),
    all(rows$total_assets > 0))
  key <- paste(as.character(rows$tophold), as.character(rows$year), sep = ":")
  groups <- split(rows$total_assets, key)
  if (any(vapply(groups, function(x) max(x) != min(x), logical(1))))
    stop("Institution-year assets disagree across cases; do not sum repeated market records.")
  unique_rows <- rows[!duplicated(key), required, drop = FALSE]
  unique_rows$key <- key[!duplicated(key)]
  log_assets <- log(unique_rows$total_assets)
  center <- mean(log_assets); scale <- sd(log_assets)
  stopifnot(is.finite(center), is.finite(scale), scale > 0)
  unique_rows$institution_log_assets_z <- (log_assets - center) / scale
  list(rows = unique_rows, center = center, scale = scale,
    definition = "Full-institution assets per bank-year; unique observations across every banking case; no market normalization.")
}

av_inputs <- function(paper, template, config) {
  input <- ds_inputs(paper, template)
  e <- new.env(parent = emptyenv())
  load(file.path(paper, "data/analysis/supreme_data_operating_v3.RData"), e)
  assets <- av_asset_index(e$simdata)
  for (name in names(input)) {
    rows <- input[[name]]$rows
    positions <- match(paste(rows$group, rows$year, sep = ":"), assets$rows$key)
    stopifnot(!anyNA(positions))
    rows$total_assets <- assets$rows$total_assets[positions]
    rows$institution_log_assets_z <- assets$rows$institution_log_assets_z[positions]
    input[[name]]$rows <- rows
    input[[name]]$data$asset_variance_active <- 0L
    input[[name]]$data$institution_log_assets_z <- rows$institution_log_assets_z
    input[[name]]$data$asset_gamma_prior_sd <- config$gamma_prior_sd
  }
  list(inputs = input, assets = assets)
}

av_quality_state <- function(data, predictor, tau, sigma, sigma_margin, rho) {
  sigma <- rep_len(sigma, data$N)
  stopifnot(length(predictor) == data$N, all(is.finite(sigma)), all(sigma > 0))
  matched <- as.logical(data$margin_obs)
  residual <- as.numeric(data$log_shareIn) - predictor
  variance <- sigma^2
  residual[matched] <- residual[matched] - rho * sigma[matched] / sigma_margin *
    (data$log_marginNormInv[matched] - log(data$pred_margin[matched] + 1e-9))
  variance[matched] <- variance[matched] * (1 - rho^2)
  q <- as.numeric(ds_private(".relative_share_quality_draw")(
    residual, variance, as.integer(data$obs_to_quality),
    rep.int(seq_along(data$quality_block_n), data$quality_block_n),
    as.integer(data$market_year_idx), tau))
  list(persistent_quality = q, period_shock = rep(0, data$N),
    structural_state = q[data$obs_to_quality])
}

# Reuse the exact existing equilibrium, index, validity and table calculations.
# Only the extracted fields and conditional quality covariance differ. Inspect
# the generated evaluator in each run alongside the sealed source module.
av_make_evaluator <- function() {
  evaluator <- ds_evaluate
  changed <- c(extraction = 0L, state = 0L, label = 0L)
  walk <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("if")) &&
        length(x) == 3L && identical(x[[2L]], as.name("candidate")) &&
        is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("<-")) &&
        identical(x[[3L]][[2L]], as.name("variables"))) {
      changed[["extraction"]] <<- changed[["extraction"]] + 1L
      return(quote(if (candidate) variables <- c(variables, "sigma_share_observation", "asset_gamma")))
    }
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("state"))) {
      changed[["state"]] <<- changed[["state"]] + 1L
      return(quote(state <- tryCatch(av_quality_state(data,
        draws$relative_share_predictor[i, ], as.numeric(draws$sigma_b_strat_economic[i]),
        if (candidate) as.numeric(draws$sigma_share_observation[i, ]) else as.numeric(draws$sigma_share[i]),
        as.numeric(draws$sigma_margin[i]), as.numeric(draws$rho[i])),
        error = function(e) list(recovery_error = conditionMessage(e)))))
    }
    if (is.character(x) && identical(x, "diagnostic_joint_quality_period_shock")) {
      changed[["label"]] <<- changed[["label"]] + 1L
      return("diagnostic_asset_heteroskedastic_quality")
    }
    if (is.call(x) || is.pairlist(x)) for (j in seq_along(x)) {
      if (is.symbol(x[[j]]) && identical(x[[j]], quote(expr = ))) next
      x[[j]] <- walk(x[[j]])
    }
    x
  }
  body(evaluator) <- walk(body(evaluator))
  stopifnot(all(changed == 1L))
  environment(evaluator) <- environment()
  evaluator
}

av_diagnose <- function(fit, candidate) {
  variables <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho")
  if (candidate) variables <- c(variables, "asset_gamma")
  native <- fit$metadata()$model_params
  if (length(native)) variables <- unique(c(variables, native))
  s <- fit$summary(variables = variables); d <- fit$diagnostic_summary()
  # CmdStan exposes transformed/GQ names too in some versions. Constant
  # structural zeros have no defined R-hat and are reported separately.
  finite_rhat <- is.finite(s$rhat)
  shape <- dim(fit$draws(variables = "a_market", format = "draws_array"))
  data.frame(Chains = shape[[2L]], Retained_Draws = shape[[1L]] * shape[[2L]],
    Max_Rhat = max(s$rhat, na.rm = TRUE), Min_Bulk_ESS = min(s$ess_bulk, na.rm = TRUE),
    Min_Tail_ESS = min(s$ess_tail, na.rm = TRUE),
    Divergences = sum(d$num_divergent), Treedepth_Hits = sum(d$num_max_treedepth),
    Diagnostic_Variables = nrow(s), Undefined_Rhat = sum(!finite_rhat))
}

av_parameters <- function(fit, input, candidate, output) {
  names <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho")
  if (candidate) names <- c(names, "asset_gamma", "sigma_share_observation")
  draw <- ds_private(".bayes_backend_extract")(fit, variables = names)
  summary <- function(x) c(Mean = mean(x), Lower = unname(quantile(x, .025)), Upper = unname(quantile(x, .975)))
  p <- lapply(c("sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho"), function(n)
    data.frame(Parameter = n, t(summary(as.numeric(draw[[n]])))))
  p[[length(p) + 1L]] <- data.frame(Parameter = "alpha_original_units", t(summary(as.numeric(draw$a_market) / input$data$price_sd)))
  if (candidate) p[[length(p) + 1L]] <- data.frame(Parameter = "asset_gamma", t(summary(as.numeric(draw$asset_gamma))))
  write.csv(do.call(rbind, p), file.path(output, "parameter_summary.csv"), row.names = FALSE)
  if (candidate) {
    gamma <- as.numeric(draw$asset_gamma)
    target <- which(input$rows$year == "1960")
    size <- input$rows[target, c("group", "bank", "total_assets", "institution_log_assets_z")]
    size$Mean_Residual_SD <- colMeans(draw$sigma_share_observation[, target, drop = FALSE])
    write.csv(size, file.path(output, "asset_residual_scales.csv"), row.names = FALSE)
    write.csv(data.frame(Pr_Gamma_Negative = mean(gamma < 0)), file.path(output, "gamma_sign.csv"), row.names = FALSE)
  }
}
