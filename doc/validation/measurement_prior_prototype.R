# Source demand_shock_prototype.R first. Its fixed-.10 experiment is preserved.
mp_kernel <- function(source, destination) {
  stage <- tempfile(fileext = ".stan")
  on.exit(unlink(stage), add = TRUE)
  ds_kernel(source, stage)
  text <- paste(readLines(stage, warn = FALSE), collapse = "\n")
  once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  once("  vector[N] shareIn;", paste(
    "  real measurement_prior_log_mean;",
    "  real<lower=0> measurement_prior_log_sd;",
    "  vector[N] shareIn;", sep = "\n"))
  once("  real log_sigma_share;", paste(
    "  real log_sigma_measurement;",
    "  real log_sigma_share;", sep = "\n"))
  once("real diagnostic_measurement_sd = 0.10;",
    "real diagnostic_measurement_sd = exp(log_sigma_measurement);")
  once("model {", paste("model {",
    "  // Primitive log-SD: normalized Normal density; no exponential Jacobian.",
    "  target += normal_lpdf(log_sigma_measurement | measurement_prior_log_mean, measurement_prior_log_sd);",
    sep = "\n"))
  once("// Diagnostic prototype: independent measurement SD .10; period shock retains residual prior.",
    "// Diagnostic: estimate independent measurement SD with a proper lognormal prior; other fixed-.10 prototype priors unchanged.")
  writeLines(text, destination)
  invisible(destination)
}

mp_inputs <- function(paper, template, config) {
  inputs <- ds_inputs(paper, template)
  stopifnot(is.finite(config$measurement_prior_log_mean),
    is.finite(config$measurement_prior_log_sd), config$measurement_prior_log_sd > 0)
  for (name in names(inputs)) {
    inputs[[name]]$data$measurement_prior_log_mean <- config$measurement_prior_log_mean
    inputs[[name]]$data$measurement_prior_log_sd <- config$measurement_prior_log_sd
  }
  inputs
}

mp_base_data <- function(data) {
  data[c("measurement_prior_log_mean", "measurement_prior_log_sd")] <- NULL
  data
}

mp_joint_state <- function(data, predictor, tau, omega, sigma_margin,
                           rho_shock, pred_margin, measurement_sd) {
  stopifnot(length(measurement_sd) == 1L, is.finite(measurement_sd), measurement_sd >= 0,
    is.finite(omega), omega > 0, is.finite(rho_shock), abs(rho_shock) < 1)
  state <- ds_joint_state(data, predictor, tau, omega, sigma_margin,
    rho_shock, pred_margin, measurement_sd = measurement_sd)
  state$measurement_sd <- measurement_sd
  state$period_shock_sd <- omega
  state$rho_shock <- rho_shock
  state
}

# Only change extracted fields and the measurement SD passed to the existing
# state recovery. Ownership, costs, FOCs, signed rates, solved-pre denominators,
# invalid-draw retention and the economic solver are exactly the old evaluator.
mp_make_evaluator <- function(observed = FALSE) {
  evaluator <- ds_evaluate
  changed <- c(extraction = 0L, state = 0L, label = 0L)
  walk <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("if")) && length(x) == 3L &&
        identical(x[[2L]], as.name("candidate")) && is.call(x[[3L]]) &&
        identical(x[[3L]][[1L]], as.name("<-")) &&
        identical(x[[3L]][[2L]], as.name("variables"))) {
      changed[["extraction"]] <<- changed[["extraction"]] + 1L
      return(quote(if (candidate) variables <- c(variables,
        "sigma_demand_shock", "rho_demand_shock", "diagnostic_measurement_sd")))
    }
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("state"))) {
      changed[["state"]] <<- changed[["state"]] + 1L
      if (observed) return(quote(state <- tryCatch(mp_joint_state(data,
        draws$relative_share_predictor[i, ], as.numeric(draws$sigma_b_strat_economic[i]),
        as.numeric(draws$sigma_share[i]), as.numeric(draws$sigma_margin[i]),
        as.numeric(draws$rho[i]), pred, measurement_sd = 0),
        error = function(e) list(recovery_error = conditionMessage(e)))))
      return(quote(state <- tryCatch(mp_joint_state(data,
        draws$relative_share_predictor[i, ], as.numeric(draws$sigma_b_strat_economic[i]),
        as.numeric(draws$sigma_demand_shock[i]), as.numeric(draws$sigma_margin[i]),
        as.numeric(draws$rho_demand_shock[i]), pred,
        measurement_sd = as.numeric(draws$diagnostic_measurement_sd[i])),
        error = function(e) list(recovery_error = conditionMessage(e)))))
    }
    if (is.character(x) && identical(x, "diagnostic_joint_quality_period_shock")) {
      changed[["label"]] <<- changed[["label"]] + 1L
      return("diagnostic_estimated_measurement_joint_demand")
    }
    if (is.character(x) && identical(x, "conditional_persistent_quality") && observed)
      return("explicit_observed_full_residual_benchmark")
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

mp_parameters <- function(fit, input, output, config) {
  variables <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share",
    "sigma_margin", "rho", "sigma_demand_shock", "rho_demand_shock",
    "diagnostic_measurement_sd", "log_sigma_measurement")
  d <- ds_private(".bayes_backend_extract")(fit, variables = variables)
  summarize <- function(x) c(Mean = mean(x), SD = sd(x),
    Lower = unname(quantile(x, .025)), Median = median(x), Upper = unname(quantile(x, .975)))
  fields <- lapply(variables[-1L], function(v) as.numeric(d[[v]]))
  names(fields) <- variables[-1L]
  fields$alpha_original_units <- as.numeric(d$a_market) / input$data$price_sd
  fields$measurement_variance_fraction <- d$diagnostic_measurement_sd^2 / d$sigma_share^2
  fields$shock_margin_covariance_per_margin_sd <- d$rho_demand_shock * d$sigma_demand_shock
  p <- do.call(rbind, lapply(names(fields), function(v)
    data.frame(Parameter = v, t(summarize(as.numeric(fields[[v]]))))))
  write.csv(p, file.path(output, "parameter_summary.csv"), row.names = FALSE)
  mu <- config$measurement_prior_log_mean; sdlog <- config$measurement_prior_log_sd
  m <- as.numeric(d$diagnostic_measurement_sd)
  prior <- data.frame(Parameter = "measurement_sd", Mean = exp(mu + sdlog^2 / 2),
    SD = sqrt((exp(sdlog^2) - 1) * exp(2 * mu + sdlog^2)),
    Lower = qlnorm(.025, mu, sdlog), Median = exp(mu), Upper = qlnorm(.975, mu, sdlog),
    Pr_Below_05 = plnorm(.05, mu, sdlog), Pr_Below_10 = plnorm(.10, mu, sdlog),
    Pr_Above_20 = 1 - plnorm(.20, mu, sdlog), Distribution = "prior")
  post <- data.frame(Parameter = "measurement_sd", t(summarize(m)),
    Pr_Below_05 = mean(m < .05), Pr_Below_10 = mean(m < .10),
    Pr_Above_20 = mean(m > .20), Distribution = "posterior")
  write.csv(rbind(prior, post), file.path(output, "measurement_prior_posterior.csv"), row.names = FALSE)
  correlations <- cor(cbind(log_measurement = as.numeric(d$log_sigma_measurement),
    log_period_shock = log(as.numeric(d$sigma_demand_shock)),
    rho_shock = as.numeric(d$rho_demand_shock), log_total_sd = log(as.numeric(d$sigma_share)),
    quality_sd = as.numeric(d$sigma_b_strat_economic)))
  write.csv(correlations, file.path(output, "variance_split_correlations.csv"))
  write.csv(data.frame(Pr_Abs_Rho_Shock_Above_90 = mean(abs(d$rho_demand_shock) > .90)),
    file.path(output, "correlation_boundary_mass.csv"), row.names = FALSE)
  invisible(p)
}
