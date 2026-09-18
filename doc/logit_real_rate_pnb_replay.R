# Replay the archived no-cost-change merger with saved posterior/quality/k
# plans; no sampling, bridge, or prior tuning. Separate unit tests cover cost
# and combined ownership/cost changes.
suppressPackageStartupMessages(library(antitrustBayes))
args <- commandArgs(trailingOnly = TRUE)
mq <- normalizePath(if (length(args)) args[[1L]] else
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality",
  mustWork = TRUE)
out <- if (length(args) >= 2L) args[[2L]] else file.path(mq,
  paste0("real-rate-replay-", format(Sys.time(), "%Y%m%dT%H%M%S")))
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
archive <- file.path(mq, "pnb-archived-retry-01")
policy <- file.path(archive, "postprocess-run-01/analysis/policy")
oracle_path <- if (length(args) >= 3L) args[[3L]] else file.path(
  mq, "logit-real-rate-01/oracle/pnb-oracle-matrices.rds")
oracle <- readRDS(oracle_path)$oracle_long
all_summary <- all_draws <- list()
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
maxabs <- function(x) max(abs(x))
for (case in c("CF", "CS")) {
  fit <- readRDS(file.path(archive, paste0(case, "-fit.rds")))
  plan <- readRDS(file.path(policy, paste0(case, "-draw-plan.rds")))
  status <- readRDS(file.path(policy, paste0(case, "-policy-status.rds")))
  labels <- as.character(fit@data$group)
  stopifnot(identical(status$acquirer, "298816"),
            identical(status$target, "889018"), length(labels) == 39L)
  ownership <- labels
  ownership[ownership == status$target] <- status$acquirer
  scenario <- antitrust::counterfactual(ownership = ownership)
  for (anchor in c("observed", "latent_quality")) {
    key <- paste(case, anchor, sep = "-")
    old <- readRDS(file.path(policy, paste0(key, "-sim.rds")))
    unit <- names(old@sims)[[1L]]
    first <- old@sims[[unit]][[which(is.na(old@failures[[unit]]))[[1L]]]]
    # Latent Logit replacement updates @shares while preserving @mktSize and
    # @insideSize.  Reconstruct from its original inside budget, not from the
    # realised latent inside share (which would rescale welfare quantities).
    inside <- first@insideSize
    start <- proc.time()[["elapsed"]]
    saved <- capture(simulate(fit, scenario, draw_plan = plan,
      n_draws = length(plan$draw_indices), demand_anchoring = anchor,
      price_domain = "real", mc.cores = 1L, labels = labels,
      insideSize = inside))
    elapsed <- proc.time()[["elapsed"]] - start
    sim <- saved$value
    saveRDS(sim, file.path(out, paste0(key, "-sim.rds")))
    saveRDS(saved$warnings, file.path(out, paste0(key, "-warnings.rds")))
    stopifnot(identical(sim@metadata$draw_indices, plan$draw_indices),
              identical(sim@metadata$k_cache, old@metadata$k_cache),
              identical(sim@product_labels, old@product_labels))
    rows <- lapply(seq_along(plan$draw_indices), function(i) {
      model <- sim@sims[[unit]][[i]]
      old_model <- old@sims[[unit]][[i]]
      original_failed <- !is.na(old@failures[[unit]][[i]])
      failed <- !is.na(sim@failures[[unit]][[i]])
      base <- data.frame(case = case, anchor = anchor, position = i,
        draw_index = plan$draw_indices[[i]],
        k = plan$k_cache[[as.character(plan$draw_indices[[i]])]][[unit]],
        original_failed = original_failed, real_failed = failed)
      if (failed) return(cbind(base, min_pre_rate = NA_real_, min_post_rate = NA_real_,
        max_pre_foc = NA_real_, max_post_foc = NA_real_,
        old_price_difference = NA_real_, old_share_difference = NA_real_,
        old_welfare_difference = NA_real_, quality_costs_preserved = NA,
        failure = sim@failures[[unit]][[i]]))
      pre_foc <- model@mcPre - model@pricePre -
        antitrust::calcMargins(model, preMerger = TRUE, level = TRUE)
      post_foc <- model@mcPost - model@pricePost -
        antitrust::calcMargins(model, preMerger = FALSE, level = TRUE)
      welfare <- c(as.numeric(antitrust::CV(model)),
                   as.numeric(antitrust::calcProducerSurplus(model, preMerger = TRUE)),
                   as.numeric(antitrust::calcProducerSurplus(model, preMerger = FALSE)))
      stopifnot(all(is.finite(welfare)),
                identical(as.numeric(model@mcPre), as.numeric(model@mcPost)))
      price_diff <- share_diff <- welfare_diff <- NA_real_
      if (!original_failed) {
        price_diff <- maxabs(c(model@pricePre - old_model@pricePre,
                               model@pricePost - old_model@pricePost))
        share_diff <- maxabs(c(sim@shares_pre[[unit]][i, ] - old@shares_pre[[unit]][i, ],
                               sim@shares_post[[unit]][i, ] - old@shares_post[[unit]][i, ]))
        old_welfare <- c(as.numeric(antitrust::CV(old_model)),
                        as.numeric(antitrust::calcProducerSurplus(old_model, preMerger = TRUE)),
                        as.numeric(antitrust::calcProducerSurplus(old_model, preMerger = FALSE)))
        welfare_diff <- maxabs(welfare - old_welfare)
        stopifnot(maxabs(model@mcPre - old_model@mcPre) < 1e-9)
        stopifnot(identical(as.numeric(model@insideSize),
                            as.numeric(old_model@insideSize)),
                  identical(as.numeric(model@mktSize),
                            as.numeric(old_model@mktSize)))
      }
      quality_preserved <- TRUE
      if (anchor == "latent_quality") {
        diag <- sim@metadata$latent_diagnostics[[unit]][[i]]
        quality_preserved <- identical(as.numeric(diag$quality),
                                        as.numeric(plan$quality_draws[i, ])) &&
          maxabs(diag$mc_pre - model@mcPre) < 1e-9
      }
      cbind(base, min_pre_rate = min(model@pricePre), min_post_rate = min(model@pricePost),
        max_pre_foc = maxabs(pre_foc), max_post_foc = maxabs(post_foc),
        old_price_difference = price_diff, old_share_difference = share_diff,
        old_welfare_difference = welfare_diff, quality_costs_preserved = quality_preserved,
        failure = NA_character_)
    })
    ledger <- do.call(rbind, rows)
    reference <- oracle[oracle$case == case & oracle$anchor == anchor, , drop = FALSE]
    stopifnot(identical(as.integer(reference$draw_position),
                        rep(seq_along(plan$draw_indices), each = length(labels))),
              identical(as.character(reference$product_label),
                        rep(labels, length(plan$draw_indices))),
              identical(as.integer(reference$draw_index),
                        rep(as.integer(plan$draw_indices), each = length(labels))))
    oracle_matrix <- function(field) matrix(reference[[field]],
      nrow = length(plan$draw_indices), ncol = length(labels), byrow = TRUE)
    oracle_price_error <- maxabs(c(sim@price_pre[[unit]] - oracle_matrix("pre_price"),
                                   sim@price_post[[unit]] - oracle_matrix("post_price")))
    oracle_share_error <- maxabs(c(sim@shares_pre[[unit]] - oracle_matrix("pre_share"),
                                   sim@shares_post[[unit]] - oracle_matrix("post_share")))
    oracle_value_error <- if (all(!ledger$real_failed)) maxabs(
      do.call(rbind, lapply(sim@sims[[unit]], function(x) as.numeric(x@mcPre))) -
        oracle_matrix("value_primitive")) else NA_real_
    write.csv(ledger, file.path(out, paste0(key, "-draw-validation.csv")), row.names = FALSE)
    summary <- data.frame(case = case, anchor = anchor, n_draws = nrow(ledger),
      positive_failures = sum(ledger$original_failed), real_failures = sum(ledger$real_failed),
      real_nonpositive_pre = sum(ledger$min_pre_rate <= 0, na.rm = TRUE),
      real_nonpositive_post = sum(ledger$min_post_rate <= 0, na.rm = TRUE),
      max_pre_foc = max(ledger$max_pre_foc, na.rm = TRUE),
      max_post_foc = max(ledger$max_post_foc, na.rm = TRUE),
      max_old_price_difference = max(ledger$old_price_difference, na.rm = TRUE),
      max_old_share_difference = max(ledger$old_share_difference, na.rm = TRUE),
      max_old_welfare_difference = max(ledger$old_welfare_difference, na.rm = TRUE),
      max_oracle_price_difference = oracle_price_error,
      max_oracle_share_difference = oracle_share_error,
      max_oracle_value_difference = oracle_value_error,
      elapsed_seconds = elapsed, all_quality_costs_preserved =
        all(ledger$quality_costs_preserved %in% TRUE))
    all_summary[[key]] <- summary
    all_draws[[key]] <- ledger
    write.csv(do.call(rbind, all_summary), file.path(out, "replay-summary.csv"), row.names = FALSE)
    print(summary)
    stopifnot(summary$real_failures == 0L, summary$all_quality_costs_preserved,
              summary$max_pre_foc < 1e-6, summary$max_post_foc < 1e-6,
              summary$max_old_price_difference < 1e-8,
              summary$max_old_share_difference < 1e-9,
              summary$max_old_welfare_difference < 1e-7,
              summary$max_oracle_price_difference < 1e-8,
              summary$max_oracle_share_difference < 1e-9,
              summary$max_oracle_value_difference < 1e-9)
    rm(old, sim, saved, first, rows); gc()
  }
  rm(fit); gc()
}
write.csv(do.call(rbind, all_draws), file.path(out, "replay-all-draws.csv"), row.names = FALSE)
