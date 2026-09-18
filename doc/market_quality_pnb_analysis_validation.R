#!/usr/bin/env Rscript
# Lightweight adversarial contracts; does not load fits or run Stan.
source_path <- "doc/market_quality_pnb_analysis.R"
reference <- new.env()
for (expression in parse(source_path))
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in%
      c("mq_stop", "mq_bind_rows", "mq_k_decomposition", "mq_sim_arrays",
        "mq_complete_evidence", "mq_validate_policy_pool"))
    eval(expression, reference)
joined <- reference$mq_bind_rows(list(
  data.frame(case = "BF", status = "available", k = 1L, posterior = .4),
  data.frame(case = "CS", status = "error", error = "retained failure")))
stopifnot(nrow(joined) == 2L, joined$posterior[[1L]] == .4,
          is.na(joined$posterior[[2L]]), joined$error[[2L]] == "retained failure")
plain <- list(long = data.frame(price_change = c(.1, .2),
  failure = NA_character_, k = NA_real_, product_index = 1L, product_label = "p1"))
rich <- list(long = data.frame(price_change = c(.1, .3, .2, .4),
  failure = NA_character_, k = c(1L, 1L, 2L, 2L), product_index = 1L, product_label = "p1"))
decomposition <- reference$mq_bind_rows(list(
  reference$mq_k_decomposition("B", plain), reference$mq_k_decomposition("BF", rich)))
stopifnot(nrow(decomposition) == 3L,
  decomposition$status[[1L]] == "no_realized_k_cells",
  abs(decomposition$within_k_variance[[2L]] - .01) < 1e-14,
  abs(decomposition$between_k_mean_variance[[2L]] - .0025) < 1e-14,
  abs(decomposition$total_variance[[2L]] - .0125) < 1e-14)
methods::setClass("MQReferenceSimulation", slots = c(price_change = "list",
  product_labels = "list", failures = "list", metadata = "list", n_draws = "integer"))
simulation <- methods::new("MQReferenceSimulation",
  price_change = list(m1_1960 = matrix(c(.1, .2, .3, .4), 2L)),
  product_labels = list(m1_1960 = c("p1", "p2")),
  failures = list(m1_1960 = c(NA_character_, "retained solver failure")),
  metadata = list(draw_indices = c(101L, 205L), k_cache = list(
    `101` = c(m1_1960 = 1L), `205` = c(m1_1960 = 2L))), n_draws = 2L)
reference$mq_groups <- c("p1", "p2")
plan <- list(metadata = list(), draw_indices = c(101L, 205L))
arrays <- reference$mq_sim_arrays(simulation, "BF", "observed", plan)
stopifnot(identical(arrays$k, c(1L, 2L)),
  all(arrays$long$k[arrays$long$draw_index == 101L] == 1L),
  all(arrays$long$failure[arrays$long$draw_index == 205L] == "retained solver failure"),
  all(arrays$product_summary$n_valid == 1L), all(arrays$product_summary$n_failed == 1L))
reference$mq_validate_policy_pool(arrays, "BF", expected_count = 2L)
malformed <- arrays
malformed$long <- rbind(malformed$long[-1L, ], malformed$long[2L, ])
error <- tryCatch(reference$mq_validate_policy_pool(malformed, "BF", expected_count = 2L),
                  error = identity)
stopifnot(inherits(error, "error"), grepl("one row per draw/product", conditionMessage(error)))
evidence <- data.frame(hypothesis = c("B", "BF"), logml_restricted = c(1, 2),
                      posterior_probability = c(.4, .6))
stopifnot(reference$mq_complete_evidence(evidence, c("BF", "B"), "complete"),
          !reference$mq_complete_evidence(evidence, c("BF", "B"), "incomplete_evidence"),
          !reference$mq_complete_evidence(evidence[c(1L, 1L), ], c("BF", "B"), "complete"))
bad <- simulation
bad@failures <- list(m1_1960 = NA_character_)
error <- tryCatch(reference$mq_sim_arrays(bad, "BF", "observed", plan), error = identity)
stopifnot(inherits(error, "error"), grepl("failure vector", conditionMessage(error)))
bad <- simulation
bad@metadata$draw_indices <- c(205L, 101L)
error <- tryCatch(reference$mq_sim_arrays(bad, "BF", "observed", plan), error = identity)
stopifnot(inherits(error, "error"), grepl("shared draw plan", conditionMessage(error)))
bad <- simulation
bad@product_labels <- list(m1_1960 = c("p2", "p1"))
error <- tryCatch(reference$mq_sim_arrays(bad, "BF", "observed", plan), error = identity)
stopifnot(inherits(error, "error"), grepl("ordered products", conditionMessage(error)))
bad <- simulation
bad@metadata$k_cache <- NULL
error <- tryCatch(reference$mq_sim_arrays(bad, "BF", "observed", plan), error = identity)
stopifnot(inherits(error, "error"), grepl("missing its realized role cache", conditionMessage(error)))
for (invalid_k in c(0L, 6L)) {
  bad <- simulation
  bad@metadata$k_cache[["101"]][["m1_1960"]] <- invalid_k
  error <- tryCatch(reference$mq_sim_arrays(bad, "BF", "observed", plan), error = identity)
  stopifnot(inherits(error, "error"), grepl("outside the frozen PNB support", conditionMessage(error)))
}
bad_pool <- arrays
bad_pool$failures[[1L]] <- ""
bad_pool$long$failure[bad_pool$long$draw_position == 1L] <- ""
error <- tryCatch(reference$mq_validate_policy_pool(bad_pool, "BF", expected_count = 2L),
                  error = identity)
stopifnot(inherits(error, "error"), grepl("failure labels", conditionMessage(error)))
simulation@metadata$k_cache[["205"]] <- c(other_unit = 2L)
error <- tryCatch(reference$mq_sim_arrays(simulation, "BF", "observed", plan),
                  error = identity)
stopifnot(inherits(error, "error"), grepl("missing posterior draw/unit", conditionMessage(error)))
cat("PNB analysis pure utility fixtures passed: mixed status, role cache, retained failures, variance decomposition.\n")
