# Separate Bertrand-only diagnostic. No production kernel or fit is migrated.
ds_private <- function(name) get(name, envir = asNamespace("antitrustBayes"))

ds_kernel <- function(source, destination) {
  original <- readLines(source, warn = FALSE)
  text <- paste(original, collapse = "\n")
  replace_once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  replace_once("real sigma_share = exp(log_sigma_share_economic);", paste(
    "real sigma_demand_shock = exp(log_sigma_share_economic);",
    "real diagnostic_measurement_sd = 0.10;",
    "real sigma_share = sqrt(square(sigma_demand_shock) + square(diagnostic_measurement_sd));",
    "real rho_demand_shock = 0.0;", sep = "\n  "))
  replace_once("if (rho_active == 1) rho = 0.95 * tanh(rho_raw[1]);", paste(
    "if (rho_active == 1) {",
    "  rho_demand_shock = 0.95 * tanh(rho_raw[1]);",
    "  rho = rho_demand_shock * sigma_demand_shock / sigma_share;",
    "}", sep = "\n  "))
  replace_once("(prior_lkj - 1.0) * log1m(square(rho))",
    "(prior_lkj - 1.0) * log1m(square(rho_demand_shock))")
  writeLines(c("// Diagnostic prototype: independent measurement SD .10; period shock retains residual prior.", text), destination)
  invisible(destination)
}

ds_validate_data <- function(data) {
  stopifnot(data$share_basis == 1L, data$quality_parameterization == 2L,
    data$supply_model == 1L, data$use_cutoff == 0L,
    data$use_market_resid_scale == 0L, data$share_likelihood == 2L,
    data$margin_likelihood == 2L, data$N_market == 1L)
  invisible(data)
}

ds_saved_fit <- function(seal) {
  records <- read.csv(seal, stringsAsFactors = FALSE)
  record <- records[records$Kind == "fit", , drop = FALSE]
  stopifnot(nrow(record) == 1L, file.exists(record$Path),
    identical(digest::digest(file = record$Path, algo = "sha256"), record$SHA256))
  fit <- readRDS(record$Path)
  stopifnot(fit@supply_model == 1L, fit@market_side == "input",
    fit@metadata$model_spec_version == 14L,
    fit@metadata$quality_scale_parameterization_effective == "total_fraction")
  ds_validate_data(fit@stan_data)
  ds_private(".cmdstan_require_policy_inference")(fit, FALSE)
  ds_private(".bayes_cmdstan_evidence_contract")(fit)
  fields <- fit@metadata$stan_data_fields
  stopifnot(length(fields) > 0L, identical(digest::digest(
    fit@stan_data[fields], algo = "sha256"), fit@metadata$stan_data_hash))
  list(fit = fit, record = record)
}

ds_alias_check <- function(seal, output) {
  saved <- ds_saved_fit(seal)
  fit <- saved$fit
  raw <- ds_private(".bayes_backend_extract")(fit@stanfit,
    variables = c("sigma_b_strat", "sigma_b_strat_economic", "log_sigma_share"))
  extracted <- ds_private(".posterior_extract")(fit, "policy_plan")
  tau <- as.numeric(raw$sigma_b_strat_economic)
  reconstructed <- as.numeric(raw$sigma_b_strat) * sqrt(plogis(-as.numeric(raw$log_sigma_share)))
  stopifnot(isTRUE(all.equal(as.numeric(extracted$sigma_b_strat), tau)),
    isTRUE(all.equal(as.numeric(extracted$sigma_b_strat_economic), tau)),
    # Historical CmdStan CSVs retain eight significant digits. This check
    # compares two rounded parameterizations, not two in-memory aliases.
    max(abs((tau - reconstructed) / tau)) < 1e-6)
  plan <- antitrustBayes::posterior_policy_draws(fit, n_draws = 10L, seed = 10701L)
  run <- function() suppressWarnings(antitrustBayes::simulate(fit,
    as.character(fit@data$group), draw_plan = plan, mc.cores = 1L,
    price_domain = "real", demand_anchoring = "latent_quality"))
  a <- run(); b <- run()
  stopifnot(all(a@n_valid == 10L), identical(a@price_pre, b@price_pre),
    identical(a@price_post, b@price_post),
    max(abs(unlist(a@price_change))) < 1e-8)
  foc <- 0
  for (models in a@sims) for (model in models) for (pre in c(TRUE, FALSE)) {
    shares <- antitrust::calcShares(model, preMerger = pre)
    own <- if (pre) model@ownerPre else model@ownerPost
    prices <- if (pre) model@pricePre else model@pricePost
    values <- if (pre) model@mcPre else model@mcPost
    residual <- values - prices - 1 / (abs(model@slopes$alpha) *
      (1 - as.numeric(own %*% shares)))
    foc <- max(foc, max(abs(residual)))
  }
  stopifnot(foc < 1e-6, any(abs(plan$quality_draws) > 0),
    identical(digest::digest(file = saved$record$Path, algo = "sha256"),
      saved$record$SHA256))
  result <- list(status = "passed", n_draws = 10L, actual_tau_alias = TRUE,
    max_tau_transform_error = max(abs(tau - reconstructed)),
    max_noop_effect = max(abs(unlist(a@price_change))), max_independent_foc = foc,
    same_plan_exact_replay = TRUE, fit_sha256_unchanged = TRUE,
    fit_path = saved$record$Path, fit_sha256 = saved$record$SHA256,
    note = "Naming/recovery check only; does not explain share shrinkage.")
  saveRDS(list(result = result, draw_plan = plan), file.path(output, "alias_check.rds"))
  write.csv(as.data.frame(result), file.path(output, "alias_check.csv"), row.names = FALSE)
  result
}

# q | share contrasts,margins,theta first; x | SAME q,contrasts,margins,theta
# second. The measurement-error realization is excluded from demand state.
ds_joint_state <- function(data, predictor, tau, omega, sigma_margin,
                           rho_shock, pred_margin, measurement_sd = .1) {
  n <- data$N
  matched <- as.logical(data$margin_obs)
  margin_resid <- data$log_marginNormInv - log(pred_margin + 1e-9)
  shift <- numeric(n); variance <- rep(omega^2, n)
  shift[matched] <- rho_shock * omega / sigma_margin * margin_resid[matched]
  variance[matched] <- variance[matched] * (1 - rho_shock^2)
  residual <- as.numeric(data$log_shareIn) - predictor - shift
  block <- rep.int(seq_along(data$quality_block_n), data$quality_block_n)
  q <- as.numeric(ds_private(".relative_share_quality_draw")(
    residual, variance + measurement_sd^2, as.integer(data$obs_to_quality), block,
    as.integer(data$market_year_idx), tau))
  shock <- numeric(n)
  for (idx in split(seq_len(n), data$market_year_idx)) {
    if (length(idx) == 1L) next
    basis <- ds_private(".relative_share_basis")(rep(1L, length(idx)))
    x <- shift[idx] + sqrt(variance[idx]) * rnorm(length(idx))
    error <- measurement_sd * rnorm(length(idx))
    covariance <- crossprod(basis, (variance[idx] + measurement_sd^2) * basis)
    discrepancy <- crossprod(basis, data$log_shareIn[idx] - predictor[idx] -
      q[data$obs_to_quality[idx]] - x - error)
    x <- x + as.numeric(variance[idx] * basis %*% solve(covariance, discrepancy))
    shock[idx] <- as.numeric(basis %*% crossprod(basis, x))
  }
  list(persistent_quality = q, period_shock = shock,
    structural_state = q[data$obs_to_quality] + shock)
}

ds_reference_state <- function(data, predictor, tau, sigma, sigma_margin, rho) {
  matched <- as.logical(data$margin_obs)
  residual <- as.numeric(data$log_shareIn) - predictor
  variance <- rep(sigma^2, data$N)
  residual[matched] <- residual[matched] - rho * sigma / sigma_margin *
    (data$log_marginNormInv[matched] - log(data$pred_margin[matched] + 1e-9))
  variance[matched] <- variance[matched] * (1 - rho^2)
  q <- as.numeric(ds_private(".relative_share_quality_draw")(
    residual, variance, as.integer(data$obs_to_quality),
    rep.int(seq_along(data$quality_block_n), data$quality_block_n),
    as.integer(data$market_year_idx), tau))
  list(persistent_quality = q, period_shock = rep(0, data$N),
    structural_state = q[data$obs_to_quality])
}

ds_inputs <- function(paper_root, template) {
  # Read the versioned production input through its existing read-only loader;
  # never source the top-level seven-model estimation controller.
  e <- new.env(parent = globalenv())
  e$datadir <- file.path(paper_root, "data/analysis")
  e$rawdir <- file.path(paper_root, "data/raw")
  e$banking_ladder_design_version <- "banking_operating_v3"
  common <- new.env(parent = globalenv())
  sys.source(file.path(paper_root, "code/banking/banking_ladder_common.R"), common)
  e$banking_economic_market_id <- common$banking_economic_market_id
  e$market_slug <- function(event_id, year) sprintf("sc_e%02d_%d", as.integer(event_id), as.integer(year))
  for (x in parse(file.path(paper_root, "code/banking/supreme_antitrustbayes.R"))) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        as.character(x[[2L]]) %in% c("selected_ratio", "load_supreme_rows")) eval(x, e)
  }
  rows <- e$load_supreme_rows()
  rows <- rows[rows$event_id_num == 1L & rows$year_num %in% 1960:1965, ]
  rows <- rows[order(rows$year_num, -rows$quantities), ]
  target <- rows[rows$year_num == 1960L, ]
  stopifnot(nrow(target) == 39L, !anyDuplicated(paste(rows$year_num, rows$group)),
    identical(as.character(target$group), as.character(template@data$group)),
    isTRUE(all.equal(target$prices, template@data$prices)),
    isTRUE(all.equal(target$quantities, template@data$quantities)),
    isTRUE(all.equal(target$margins, template@data$margins)))
  history <- data.frame(prices = rows$prices, quantities = rows$quantities,
    margins = rows$margins, market = as.character(rows$economic_market_id),
    year = as.character(rows$year_num), group = as.character(rows$group),
    margin_observed = rows$margin_observed, outside_price = 1,
    bank = rows$bank, target = rows$target, acquirer = rows$acquirer)
  d <- ds_private("prepare_stan_data")(
    prices = history$prices, quantities = history$quantities, margins = history$margins,
    market = history$market, year = history$year, group = history$group,
    margin_observed = history$margin_observed, supply_model = 1L,
    market_side = "input", use_rho = TRUE, use_hmt = FALSE,
    use_alpha_market_pool = TRUE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = TRUE, priors = template@metadata$resolved_priors,
    quality_parameterization = "collapsed", quality_scope = "market",
    quality_scale_parameterization = "total_fraction", share_basis = "relative",
    share_likelihood = "gaussian", margin_likelihood = "gaussian")
  ds_validate_data(d)
  snapshot <- history[history$year == "1960", ]
  list(snapshot = list(data = template@stan_data, rows = snapshot),
    history = list(data = d, rows = history))
}

ds_evaluate <- function(stanfit, input, candidate, output, n_draws = 100L, seed = 10702L) {
  data <- input$data; rows <- input$rows
  ds_validate_data(data)
  variables <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share",
    "sigma_margin", "rho", "relative_share_predictor", "pred_marginNormInv")
  if (candidate) variables <- c(variables, "sigma_demand_shock", "rho_demand_shock")
  draws <- ds_private(".bayes_backend_extract")(stanfit, variables = variables)
  set.seed(seed)
  indices <- sample.int(nrow(as.matrix(draws$a_market)), n_draws)
  target <- which(rows$year == "1960")
  stopifnot(length(target) == 39L)
  observed <- rows[target, ]; observed_share <- observed$quantities / sum(observed$quantities)
  owner <- observed$group; post_owner <- owner
  post_owner[post_owner == observed$target[[1L]]] <- observed$acquirer[[1L]]
  parties <- c(observed$target[[1L]], observed$acquirer[[1L]])
  selected <- c(grep("First Pennsylvania", observed$bank, ignore.case = TRUE),
    grep("Philadelphia National", observed$bank, ignore.case = TRUE),
    grep("Girard", observed$bank, ignore.case = TRUE))
  stopifnot(length(selected) == 3L, !anyDuplicated(selected),
    sum(owner %in% parties) == 2L)
  states <- vector("list", n_draws); diagnostics <- vector("list", n_draws)
  share_observed <- share_solved <- matrix(NA_real_, n_draws, length(target))
  margin_mean <- margin_sd <- matrix(NA_real_, n_draws, length(target))
  effects <- data.frame(Draw = seq_len(n_draws), Posterior_Row = indices,
    Market_Paasche = NA_real_, Party_Paasche = NA_real_,
    FOC = NA_real_, Noop = NA_real_, Error = NA_character_,
    Market_Index_Error = NA_character_, Party_Index_Error = NA_character_)
  for (pos in seq_along(indices)) {
    i <- indices[[pos]]
    pred <- as.numeric(draws$pred_marginNormInv[i, ])
    data$pred_margin <- pred
    alpha <- as.numeric(draws$a_market[i, 1L]) / data$price_sd
    s0 <- as.numeric(draws$s0_prob[i, data$market_year_idx[target[[1L]]]])
    state <- tryCatch(if (candidate) ds_joint_state(data, draws$relative_share_predictor[i, ],
      as.numeric(draws$sigma_b_strat_economic[i]), as.numeric(draws$sigma_demand_shock[i]),
      as.numeric(draws$sigma_margin[i]), as.numeric(draws$rho_demand_shock[i]), pred) else
      ds_reference_state(data, draws$relative_share_predictor[i, ],
        as.numeric(draws$sigma_b_strat_economic[i]), as.numeric(draws$sigma_share[i]),
        as.numeric(draws$sigma_margin[i]), as.numeric(draws$rho[i])),
      error = function(e) list(recovery_error = conditionMessage(e)))
    states[[pos]] <- state
    margin_mean[pos, ] <- log(pred[target] + 1e-9)
    margin_sd[pos, ] <- as.numeric(draws$sigma_margin[i])
    if (!is.null(state$recovery_error)) {
      effects$Error[pos] <- paste("joint state recovery:", state$recovery_error)
      next
    }
    tryCatch({
      market <- antitrustBayes::bayes_market(observed$prices, observed$quantities,
        observed$margins, owner, outside_price = 1, market_side = "input",
        insideSize = sum(observed$quantities), labels = make.unique(observed$bank))
      params <- antitrustBayes::bayes_params(alpha = alpha, s0 = s0, supply_model = 1L)
      fit <- antitrustBayes::as_antitrust_fit(ds_private("build_conduct_object")(
        params, 1L, market, owner, owner))
      score <- list(score = draws$relative_share_predictor[i, target] +
        state$structural_state[target], s0 = s0,
        quality = state$structural_state[target], draw_index = i)
      realized <- ds_private(".lq_apply_to_fit")(fit, score,
        list(quality_source = if (candidate) "diagnostic_joint_quality_period_shock" else "conditional_persistent_quality",
          quality_conditioning = "relative_share_contrasts_and_observed_margins"), price_domain = "real")
      fit <- realized$fit
      share_observed[pos, ] <- realized$metadata$reference_quantity_shares /
        sum(realized$metadata$reference_quantity_shares)
      share_solved[pos, ] <- as.numeric(antitrust::calcShares(fit@model, TRUE)) /
        sum(antitrust::calcShares(fit@model, TRUE))
      diagnostics[[pos]] <- realized$metadata
      result <- suppressWarnings(antitrust::simulate(fit, post_owner, price_domain = "real"))
      noop <- suppressWarnings(antitrust::simulate(fit, owner, price_domain = "real"))
      residual <- 0
      for (pre in c(TRUE, FALSE)) {
        sh <- as.numeric(antitrust::calcShares(result, pre))
        own <- if (pre) result@ownerPre else result@ownerPost
        prices <- if (pre) result@pricePre else result@pricePost
        values <- if (pre) result@mcPre else result@mcPost
        residual <- max(residual, max(abs(values - prices -
          1 / (alpha * (1 - as.numeric(own %*% sh))))))
      }
      effects$FOC[pos] <- residual
      effects$Noop[pos] <- max(abs(noop@pricePost - noop@pricePre))
      stopifnot(residual < 1e-6, effects$Noop[pos] < 1e-6,
        identical(fit@model@slopes, result@slopes))
      index <- function(party) {
        use <- if (party) owner %in% parties else rep(TRUE, length(owner))
        weights <- as.numeric(antitrust::calcShares(result, FALSE))[use]
        pre <- result@pricePre[use]; post <- result@pricePost[use]
        denominator <- sum(weights * pre)
        if (any(!is.finite(c(weights, pre, post))) || any(weights < 0) || denominator <= 0)
          stop("Invalid Paasche basket or nonpositive solved-pre denominator")
        sum(weights * post) / denominator - 1
      }
      for (party in c(FALSE, TRUE)) tryCatch({
        effects[pos, if (party) "Party_Paasche" else "Market_Paasche"] <- index(party)
      }, error = function(e) effects[pos, if (party) "Party_Index_Error" else "Market_Index_Error"] <<- conditionMessage(e))
    }, error = function(e) effects$Error[pos] <<- conditionMessage(e))
  }
  summarize <- function(x) {
    finite <- x[is.finite(x)]
    c(Mean = if (length(finite)) mean(finite) else NA_real_,
      Lower = if (length(finite)) quantile(finite, .025, names = FALSE) else NA_real_,
      Upper = if (length(finite)) quantile(finite, .975, names = FALSE) else NA_real_,
      Valid = length(finite), Requested = length(x))
  }
  table <- do.call(rbind, lapply(seq_len(4L), function(j) {
    idx <- if (j <= 3L) selected[[j]] else selected
    weight <- observed_share[idx] / sum(observed_share[idx])
    obs <- summarize(100 * rowSums(share_observed[, idx, drop = FALSE]))
    solved <- summarize(100 * rowSums(share_solved[, idx, drop = FALSE]))
    err <- sweep(margin_mean[, idx, drop = FALSE], 2L,
      data$log_marginNormInv[target[idx]], "-")
    data.frame(Bank = if (j <= 3L) as.character(observed$bank[idx]) else "Three banks combined",
      Observed_Share_PP = 100 * sum(observed_share[idx]),
      At_Observed_Mean_PP = obs[["Mean"]], At_Observed_Lower_PP = obs[["Lower"]],
      At_Observed_Upper_PP = obs[["Upper"]], Solved_Mean_PP = solved[["Mean"]],
      Solved_Lower_PP = solved[["Lower"]], Solved_Upper_PP = solved[["Upper"]],
      At_Observed_Error_PP = obs[["Mean"]] - 100 * sum(observed_share[idx]),
      Solved_Error_PP = solved[["Mean"]] - 100 * sum(observed_share[idx]),
      Share_Valid = solved[["Valid"]], Share_Requested = n_draws,
      Margin_LogInv_Mean_Error = mean(as.numeric(err %*% weight)),
      Margin_LogInv_Point_RMSE = sqrt(mean(as.numeric((err^2) %*% weight))),
      Margin_LogInv_Predictive_RMSE = sqrt(mean(as.numeric((err^2 + margin_sd[, idx, drop = FALSE]^2) %*% weight))),
      Market_Paasche_Mean_PP = 100 * summarize(effects$Market_Paasche)[["Mean"]],
      Market_Paasche_Lower_PP = 100 * summarize(effects$Market_Paasche)[["Lower"]],
      Market_Paasche_Upper_PP = 100 * summarize(effects$Market_Paasche)[["Upper"]],
      Market_Index_Valid = sum(is.finite(effects$Market_Paasche)),
      Party_Paasche_Mean_PP = 100 * summarize(effects$Party_Paasche)[["Mean"]],
      Party_Paasche_Lower_PP = 100 * summarize(effects$Party_Paasche)[["Lower"]],
      Party_Paasche_Upper_PP = 100 * summarize(effects$Party_Paasche)[["Upper"]],
      Party_Index_Valid = sum(is.finite(effects$Party_Paasche)))
  }))
  saveRDS(list(indices = indices, states = states, diagnostics = diagnostics,
    share_at_observed = share_observed, share_solved = share_solved,
    effects = effects, summary = table), file.path(output, "policy_draws.rds"))
  write.csv(effects, file.path(output, "policy_validity.csv"), row.names = FALSE)
  write.csv(table, file.path(output, "comparison.csv"), row.names = FALSE)
  table
}
