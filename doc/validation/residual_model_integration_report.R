# Independent FOCs and draw-by-draw policy reporting; no refitting.
rmi_foc <- function(model, pre = TRUE) {
  p <- if (pre) model@pricePre else model@pricePost
  cost <- if (pre) model@mcPre else model@mcPost
  own <- if (pre) model@ownerPre else model@ownerPost
  n <- length(p)
  if (methods::is(model, "CES")) {
    gamma <- model@slopes$gamma; w <- model@slopes$meanval
    expenditure <- model@priceOutside^(1 - gamma) + sum(w * p^(1 - gamma))
    q <- w * p^-gamma / expenditure
    revenue <- w * p^(1 - gamma) / expenditure
    derivative <- -gamma * diag(q / p) +
      (gamma - 1) * outer(q, revenue / p)
  } else {
    alpha <- model@slopes$alpha
    score <- model@slopes$meanval + alpha * (p - model@priceOutside)
    # Stable outside-inclusive softmax, independent of package share methods.
    mx <- max(c(0, score)); w <- exp(score - mx)
    q <- w / (exp(-mx) + sum(w))
    derivative <- alpha * (diag(q) - tcrossprod(q))
  }
  as.numeric(q + t(derivative * own) %*% (p - cost))
}
rmi_parameter_report <- function(fit, fixture, config) {
  draws <- as.matrix(fit@stanfit$draws(format = "draws_matrix"))
  combined <- identical(fit@metadata$share_residual_model, "combined")
  total_sd <- sqrt(config$demand_shock_sd^2 + config$measurement_sd^2)
  effective_rho <- config$rho_demand * config$demand_shock_sd / total_sd
  stems <- c("alpha_structural_market", "gamma_structural_market", "s0_prob",
    "sigma_b_strat_economic", "sigma_demand_shock", "sigma_measurement", "sigma_share",
    "sigma_margin", "rho_demand_shock", "rho", "sigma_demand_shock_market", "rho_market")
  cols <- which(sub("\\[.*", "", colnames(draws)) %in% stems)
  out <- lapply(cols, function(j) {
    name <- colnames(draws)[j]; stem <- sub("\\[.*", "", name)
    truth <- switch(stem, alpha_structural_market = 3, gamma_structural_market = 6,
      s0_prob = .2, sigma_b_strat_economic = config$quality_sd,
      sigma_demand_shock = if (combined) total_sd else config$demand_shock_sd,
      sigma_measurement = if (combined) 0 else config$measurement_sd,
      sigma_share = total_sd, sigma_demand_shock_market = if (combined) total_sd else config$demand_shock_sd,
      sigma_margin = config$margin_sd,
      rho_demand_shock = if (combined) effective_rho else config$rho_demand,
      rho = effective_rho, rho_market = effective_rho, NA_real_)
    qq <- quantile(draws[, j], c(.025, .5, .975))
    data.frame(parameter = name, truth, q025 = qq[1L], median = qq[2L], q975 = qq[3L],
      interpretation = if (combined) "Combined-route equivalent covariance parameters; no fitted measurement component." else
        "Bounded recovery; variance components depend on priors and free shock correlation.")
  })
  # The identified observation variance and covariance are reported explicitly.
  omega <- draws[, "sigma_demand_shock"]; measurement <- draws[, "sigma_measurement"]
  derived <- list(total_variance = omega^2 + measurement^2,
    share_margin_covariance = draws[, "rho_demand_shock"] * omega * draws[, "sigma_margin"])
  for (name in names(derived)) {
    qq <- quantile(derived[[name]], c(.025, .5, .975))
    truth <- if (name == "total_variance") config$demand_shock_sd^2 + config$measurement_sd^2 else
      config$rho_demand * config$demand_shock_sd * config$margin_sd
    out[[length(out) + 1L]] <- data.frame(parameter = name, truth, q025 = qq[1L],
      median = qq[2L], q975 = qq[3L], interpretation = "Observation covariance target.")
  }
  do.call(rbind, out)
}

rmi_report_fit <- function(output, id) {
  manifest <- readRDS(file.path(output, "manifest.rds")); config <- manifest$config
  dest <- file.path(output, id); fit <- readRDS(file.path(dest, "fit.rds"))
  fixture <- readRDS(file.path(output, "fixtures.rds"))[[sub("__.*", "", id)]]
  write.csv(rmi_parameter_report(fit, fixture, config), file.path(dest, "parameters.csv"), row.names = FALSE)
  plan <- posterior_policy_draws(fit, n_draws = config$policy_draws, seed = config$policy_seed)
  saveRDS(plan, file.path(dest, "policy_plan.rds"))
  warnings <- character(); quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  prepared <- quiet(posterior_policy_models(fit, draw_plan = plan,
    year = unique(fixture$data$year), insideSize = 1000))
  saveRDS(prepared, file.path(dest, "baseline_models.rds"))
  results <- shares <- list()
  for (unit in seq_len(nrow(prepared$units))) {
    uid <- as.character(prepared$units$unit_id[unit])
    ix <- which(fixture$data$market == prepared$units$market[unit] &
      fixture$data$year == prepared$units$year[unit])
    dat <- fixture$data[ix, ]; truth <- fixture$truth[ix, ]; owner <- as.character(dat$group)
    post_owner <- owner; post_owner[post_owner == "b"] <- "a"
    parties <- owner %in% c("a", "b")
    for (j in seq_len(config$policy_draws)) {
      record <- data.frame(unit = uid, draw_position = j, posterior_row = plan$draw_indices[j],
        status = "failed", reason = "", pre_foc = NA_real_, post_foc = NA_real_,
        no_op_error = NA_real_, market_paasche = NA_real_, party_paasche = NA_real_,
        market_laspeyres = NA_real_, market_fisher = NA_real_)
      value <- tryCatch({
        source <- prepared$models[[uid]][[j]]
        if (is.null(source)) stop("Baseline unavailable; see preserved realization failures.")
        model <- if (methods::is(source, "AntitrustFit")) source@model else source
        record$pre_foc <- max(abs(rmi_foc(model, TRUE)))
        reference <- model; reference@pricePre <- dat$prices
        reference_share <- as.numeric(quiet(antitrust::calcShares(reference, preMerger = TRUE)))
        pre_share <- as.numeric(quiet(antitrust::calcShares(model, preMerger = TRUE)))
        no_op <- quiet(antitrust::simulate(source, antitrust::counterfactual(ownership = owner)))
        if (methods::is(no_op, "AntitrustFit")) no_op <- no_op@model
        record$no_op_error <- max(abs(no_op@pricePost - model@pricePre))
        merger <- quiet(antitrust::simulate(source, antitrust::counterfactual(ownership = post_owner)))
        if (methods::is(merger, "AntitrustFit")) merger <- merger@model
        record$post_foc <- max(abs(rmi_foc(merger, FALSE)))
        p0 <- model@pricePre; p1 <- merger@pricePost
        q1 <- as.numeric(quiet(antitrust::calcShares(merger, preMerger = FALSE)))
        record$market_paasche <- sum(p1 * q1) / sum(p0 * q1) - 1
        record$party_paasche <- sum(p1[parties] * q1[parties]) / sum(p0[parties] * q1[parties]) - 1
        record$market_laspeyres <- sum(p1 * pre_share) / sum(p0 * pre_share) - 1
        record$market_fisher <- sqrt((1 + record$market_paasche) * (1 + record$market_laspeyres)) - 1
        stopifnot(record$pre_foc < 1e-6, record$post_foc < 1e-6, record$no_op_error < 1e-6)
        shares[[length(shares) + 1L]] <- data.frame(unit = uid, draw_position = j,
          posterior_row = plan$draw_indices[j], product = seq_along(ix), group = owner,
          observed_share = dat$quantities / sum(dat$quantities),
          true_share = truth$conditional_share, at_observed_price = reference_share / sum(reference_share),
          solved_pre_share = pre_share / sum(pre_share), solved_pre_price = p0,
          observed_price = dat$prices)
        record$status <- "valid"; record
      }, error = function(e) {record$reason <- conditionMessage(e); record})
      results[[length(results) + 1L]] <- value
    }
  }
  results <- do.call(rbind, results)
  write.csv(results, file.path(dest, "policy_all_rows.csv"), row.names = FALSE)
  if (length(shares)) {
    shares <- do.call(rbind, shares)
    write.csv(shares, file.path(dest, "share_predictions_all_rows.csv"), row.names = FALSE)
    means <- aggregate(cbind(observed_share, true_share, at_observed_price, solved_pre_share) ~
      unit + product + group, shares, mean)
    means$at_observed_error_pp <- 100 * (means$at_observed_price - means$observed_share)
    means$solved_pre_error_pp <- 100 * (means$solved_pre_share - means$observed_share)
    write.csv(means, file.path(dest, "share_prediction_means.csv"), row.names = FALSE)
    firm <- aggregate(cbind(observed_share, true_share, at_observed_price, solved_pre_share) ~
      unit + group, means, sum)
    firm$at_observed_error_pp <- 100 * (firm$at_observed_price - firm$observed_share)
    firm$solved_pre_error_pp <- 100 * (firm$solved_pre_share - firm$observed_share)
    write.csv(firm, file.path(dest, "firm_share_prediction_means.csv"), row.names = FALSE)
  }
  pred <- 1 / as.matrix(fit@stanfit$draws("pred_marginNormInv", format = "draws_matrix"))
  write.csv(data.frame(row = seq_len(nrow(fixture$data)), observed_margin = fixture$data$margins,
    posterior_mean_margin = colMeans(pred), error = colMeans(pred) - fixture$data$margins),
    file.path(dest, "margin_fit.csv"), row.names = FALSE)
  write.csv(data.frame(unit = prepared$units$unit_id,
    attempted = as.integer(table(factor(results$unit, levels = prepared$units$unit_id))),
    invalid = as.integer(table(factor(results$unit[results$status != "valid"], levels = prepared$units$unit_id)))),
    file.path(dest, "invalid_counts.csv"), row.names = FALSE)
  saveRDS(list(warnings = unique(warnings), failures = prepared$failures,
    draw_indices = plan$draw_indices), file.path(dest, "policy_warnings_failures.rds"))
  invisible(results)
}
