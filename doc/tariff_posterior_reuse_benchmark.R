# Reuse saved posterior fits; this script never estimates a model.
# Run after installing the updated trade and antitrustBayes packages.
args <- commandArgs(trailingOnly = TRUE)
fixture <- if (length(args)) args[[1L]] else "logit_bertrand"
n_draws <- if (length(args) > 1L) as.integer(args[[2L]]) else 200L
base <- Sys.getenv("TARIFF_RUN_DIR", "doc/statistical_toolkit_runs/20260918/tariff-reuse-01")
old <- "doc/statistical_toolkit_runs/20260916/performance-audit/recovery"
if (dir.exists(file.path(base, "candidate-lib")))
  .libPaths(c(file.path(base, "candidate-lib"),
              file.path(old, "market-quality/logit-real-rate-01/Rlib"),
              file.path(old, "Rlib"), .libPaths()))
library(antitrustBayes)
library(trade)
paths <- c(
  logit_bertrand = file.path(old, "market-quality/rotation-rho-sensitivity-02-budget2k/B_collapsed_market_314159_budget2k-fit.rds"),
  rho_panel = file.path(old, "market-quality/rotation-rho-sensitivity-02-budget2k/B_collapsed_market_271828_budget2k-fit.rds")
)
# For another model family, supply a validated current-specification fit.
custom_fit <- Sys.getenv("TARIFF_FIT_PATH", "")
if (nzchar(custom_fit)) paths[fixture] <- custom_fit
stopifnot(fixture %in% names(paths), n_draws > 0L)
fit <- readRDS(paths[[fixture]])
stopifnot(is(fit, "BayesFit"))
original <- serialize(fit, NULL)
run_tag <- Sys.getenv("TARIFF_BENCHMARK_TAG", paste0("n", n_draws))
out <- file.path(base, "posterior-benchmark", fixture, run_tag)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[[3L]]
plan <- posterior_policy_draws(fit, n_draws = n_draws, seed = 20260918L)
plan_seconds <- proc.time()[[3L]] - t0
saveRDS(plan, file.path(out, "plan.rds"))
scenario <- list(
  merger = list(merger = TRUE, tau = 0, cost = 0),
  tariff = list(merger = FALSE, tau = .2, cost = 0),
  cost = list(merger = FALSE, tau = 0, cost = .1),
  combined = list(merger = TRUE, tau = .2, cost = .1)
)
rows <- list()
failures <- list()
models_raw <- list()
fail <- function(unit, index, anchoring, policy, stage, condition) {
  failures[[length(failures) + 1L]] <<- data.frame(
    fixture, unit, index, anchoring, policy, stage,
    category = class(condition)[[1L]], message = conditionMessage(condition))
}
for (anchoring in c("observed", "latent_quality")) {
  t0 <- proc.time()[[3L]]
  prepared <- suppressWarnings(posterior_policy_models(
    fit, draw_plan = plan, year = unique(fit@data$year),
    demand_anchoring = anchoring))
  construction_seconds <- proc.time()[[3L]] - t0
  saveRDS(prepared$failures, file.path(out, paste0(anchoring, "-construction-failures.rds")))
  # Promote each baseline once, then reuse it across all policy scenarios.
  promoted_cache <- prepared$models
  promotion_errors <- list()
  promotion_seconds <- 0
  promotion_count <- 0L
  for (uid in names(prepared$models)) for (pos in seq_along(plan$draw_indices)) {
    source <- prepared$models[[uid]][[pos]]
    if (is.null(source)) next
    key <- paste(uid, pos, sep = ":")
    t0 <- proc.time()[[3L]]
    promoted <- tryCatch({
      value <- as_trade_fit(source, tariffPre = 0,
                           cost_basis = "effective", margin_basis = "net_revenue")
      stopifnot(identical(attr(source, "bayes_policy_draw"),
                          attr(value, "bayes_policy_draw")),
                value@tariff_state$promotion_solves == 0L)
      value
    }, error = function(e) {
      promotion_errors[[key]] <<- e
      NULL
    })
    promoted_cache[[uid]][pos] <- list(promoted)
    promotion_seconds <- promotion_seconds + proc.time()[[3L]] - t0
    promotion_count <- promotion_count + !is.null(promoted)
  }
  models_raw[[anchoring]] <- list()
  for (policy in names(scenario)) {
    setting <- scenario[[policy]]
    products <- list()
    solves <- assemblies <- numeric()
    comparisons <- numeric()
    attempts <- 0L
    for (uid in names(prepared$models)) for (pos in seq_along(plan$draw_indices)) {
      index <- plan$draw_indices[[pos]]
      source <- prepared$models[[uid]][[pos]]
      attempts <- attempts + 1L
      if (is.null(source)) {
        fail(uid, index, anchoring, policy, "construction",
             simpleError("structural realization failed; see construction failures"))
        next
      }
      promoted <- promoted_cache[[uid]][[pos]]
      if (is.null(promoted)) {
        fail(uid, index, anchoring, policy, "promotion",
             promotion_errors[[paste(uid, pos, sep = ":")]])
        next
      }
      stage <- "promotion"
      tryCatch({
        owners <- promoted@tariff_state$ownerPre
        firms <- unique(owners)
        if (setting$merger && length(firms) > 1L) owners[owners == firms[[2L]]] <- firms[[1L]]
        remap <- function(roles) {
          if (setting$merger && length(firms) > 1L) roles[roles == firms[[2L]]] <- firms[[1L]]
          unique(roles)
        }
        native <- !is(source, "AntitrustFit")
        is_stack <- native && is(source, "StackelbergLogit") || native && is(source, "StackelbergCES")
        role_args <- if (!native) list() else if (is_stack)
          list(leadersPost = remap(promoted@tariff_state$leadersPre)) else
          list(corePost = remap(promoted@tariff_state$corePre))
        stage <- "policy_solve"
        t0 <- proc.time()[[3L]]
        result <- do.call(trade::simulate, c(list(
          promoted, tariffPost = setting$tau, ownerPost = owners, mcDelta = setting$cost), role_args))
        solves <- c(solves, proc.time()[[3L]] - t0)
        effective_delta <- (1 + setting$cost) / (1 - setting$tau) - 1
        stage <- "equivalent_game_solve"
        t0 <- proc.time()[[3L]]
        expected <- if (!native) antitrust::simulate(source, owners, mcDelta = rep(effective_delta, length(owners))) else if (is_stack)
          coordination::stackelberg_simulate(source, ownerPost = owners,
            leadersPost = role_args$leadersPost, mcDelta = rep(effective_delta, length(owners))) else
          coordination::core_fringe_simulate(source, ownerPost = owners,
            corePost = role_args$corePost, mcDelta = rep(effective_delta, length(owners)))
        comparisons <- c(comparisons, proc.time()[[3L]] - t0)
        expected_model <- if (is(expected, "AntitrustFit")) expected@model else expected
        stage <- "accounting"
        t0 <- proc.time()[[3L]]
        accounts <- tariff_accounts(result)
        before <- tariff_accounts(result, TRUE)
        welfare <- tariff_welfare(result)
        price_error <- max(abs(accounts$price - expected_model@pricePost))
        share_error <- max(abs(antitrust::calcShares(result, FALSE, revenue = FALSE) -
                               antitrust::calcShares(expected_model, FALSE, revenue = FALSE)))
        producer_error <- max(abs(accounts$producer_profit -
          ((1 - setting$tau) * accounts$price - promoted@tariff_state$cPre * (1 + setting$cost)) * accounts$q))
        government_error <- max(abs(accounts$government_revenue - setting$tau * accounts$price * accounts$q))
        stopifnot(price_error < 1e-7, share_error < 1e-7, producer_error < 1e-10,
                  government_error < 1e-10,
                  identical(attr(result, "bayes_policy_draw"), attr(source, "bayes_policy_draw")),
                  identical(result@tariff_state$source_model, promoted@tariff_state$source_model))
        products[[length(products) + 1L]] <- data.frame(
          fixture, anchoring, policy, uid, index, product = accounts$product,
          price_change = accounts$price / before$price - 1,
          price_error, share_error, producer_error, government_error,
          physical_ps_delta = welfare$physical_producer_surplus_delta,
          government_delta = welfare$importing_government_revenue_delta,
          consumer_change = welfare$monetary_consumer_change)
        assemblies <- c(assemblies, proc.time()[[3L]] - t0)
      }, error = function(e) fail(uid, index, anchoring, policy, stage, e))
    }
    product_data <- if (length(products)) do.call(rbind, products) else data.frame()
    models_raw[[anchoring]][[policy]] <- product_data
    unit_welfare <- if (nrow(product_data)) product_data[
      !duplicated(product_data[, c("uid", "index")]), , drop = FALSE] else data.frame()
    draw_price <- if (nrow(product_data)) stats::aggregate(
      price_change ~ index, product_data, mean)$price_change else numeric()
    draw_welfare <- if (nrow(unit_welfare)) stats::aggregate(
      cbind(physical_ps_delta, government_delta, consumer_change) ~ index,
      unit_welfare, mean) else data.frame()
    rows[[length(rows) + 1L]] <- data.frame(
      fixture, anchoring, policy, selected_draws = length(plan$draw_indices),
      units = nrow(prepared$units), attempts, usable = length(products),
      failures = attempts - length(products), plan_seconds, construction_seconds,
      promotion_seconds, promotion_count, promotion_shared = TRUE,
      solve_seconds = sum(solves),
      equivalent_game_seconds = sum(comparisons), assembly_seconds = sum(assemblies),
      price_error = if(nrow(product_data)) max(product_data$price_error) else NA_real_,
      share_error = if(nrow(product_data)) max(product_data$share_error) else NA_real_,
      producer_error = if(nrow(product_data)) max(product_data$producer_error) else NA_real_,
      government_error = if(nrow(product_data)) max(product_data$government_error) else NA_real_,
      physical_ps_delta_mean = if(nrow(unit_welfare)) mean(unit_welfare$physical_ps_delta) else NA_real_,
      physical_ps_delta_q025 = if(nrow(unit_welfare)) unname(quantile(unit_welfare$physical_ps_delta,.025)) else NA_real_,
      physical_ps_delta_q975 = if(nrow(unit_welfare)) unname(quantile(unit_welfare$physical_ps_delta,.975)) else NA_real_,
      government_delta_mean = if(nrow(unit_welfare)) mean(unit_welfare$government_delta) else NA_real_,
      consumer_monetary_change_mean = if(nrow(unit_welfare)) mean(unit_welfare$consumer_change) else NA_real_,
      price_change_mean = if(nrow(product_data)) mean(product_data$price_change) else NA_real_,
      price_change_q025 = if(nrow(product_data)) unname(quantile(product_data$price_change,.025)) else NA_real_,
      price_change_q975 = if(nrow(product_data)) unname(quantile(product_data$price_change,.975)) else NA_real_,
      draw_mean_price_change_q025 = if(length(draw_price)) unname(quantile(draw_price,.025)) else NA_real_,
      draw_mean_price_change_q975 = if(length(draw_price)) unname(quantile(draw_price,.975)) else NA_real_,
      draw_mean_physical_ps_delta_q025 = if(nrow(draw_welfare)) unname(quantile(draw_welfare$physical_ps_delta,.025)) else NA_real_,
      draw_mean_physical_ps_delta_q975 = if(nrow(draw_welfare)) unname(quantile(draw_welfare$physical_ps_delta,.975)) else NA_real_
    )
    write.csv(do.call(rbind, rows), file.path(out,"summary.csv"), row.names = FALSE)
    saveRDS(models_raw, file.path(out,"policy-draws.rds"))
    saveRDS(failures, file.path(out,"failures.rds"))
  }
}
stopifnot(identical(serialize(fit, NULL), original), identical(plan, readRDS(file.path(out,"plan.rds"))))
if(length(failures)) write.csv(do.call(rbind, failures),file.path(out,"failures.csv"),row.names=FALSE)
print(do.call(rbind,rows))
