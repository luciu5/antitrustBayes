mode <- Sys.getenv("BENCH_MODE", unset = "selected")
out_path <- Sys.getenv("BENCH_OUT", unset = "benchmark-policies.rds")
if (!mode %in% c("full", "selected")) stop("BENCH_MODE must be full or selected")

base <- "/mnt/d/projects/antitrustBayes/doc/statistical_toolkit_runs/20260916/hybrid/fits"
fixtures <- c(
  sim_logit = "logit_bertrand_noisy_271828__D.rds",
  exact_BF = "logit_BF_noisy_271828__D.rds",
  exact_BS = "logit_BS_noisy_271828__D.rds",
  sim_ces = "ces_bertrand_noisy_271828__D.rds"
)
fixture_filter <- Sys.getenv("BENCH_FIXTURES", unset = "")
if (nzchar(fixture_filter)) {
  keep <- intersect(names(fixtures), strsplit(fixture_filter, ",", fixed = TRUE)[[1L]])
  if (!length(keep)) stop("BENCH_FIXTURES selected no known fixture")
  fixtures <- fixtures[keep]
}

clock <- function(fun) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(fun(), error = function(e) e)
  list(value = value, elapsed = proc.time()[["elapsed"]] - start)
}

safe_hash <- function(x) tryCatch(digest::digest(x), error = function(e) NA_character_)
finite_max <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || !any(is.finite(x))) NA_real_ else max(abs(x[is.finite(x)]))
}
snapshot <- function(value, fit) {
  welfare <- lapply(value@sims[[1L]], function(model) {
    if (is.null(model)) return(c(cv = NA_real_, ps = NA_real_, welfare = NA_real_))
    model <- if (methods::is(model, "CounterfactualPath")) antitrust::final_result(model) else model
    tryCatch({
      cv <- -as.numeric(antitrust::CV(model))
      ps <- sum(antitrust::calcProducerSurplus(model, preMerger = FALSE)) -
        sum(antitrust::calcProducerSurplus(model, preMerger = TRUE))
      c(cv = cv, ps = ps, welfare = cv + ps)
    }, error = function(e) c(cv = NA_real_, ps = NA_real_, welfare = NA_real_))
  })
  welfare <- if (length(welfare)) do.call(rbind, welfare) else matrix(numeric(), 0L, 3L)
  list(
    price_pre = value@price_pre[[1L]], price_post = value@price_post[[1L]],
    price_change = value@price_change[[1L]], shares_pre = value@shares_pre[[1L]],
    shares_post = value@shares_post[[1L]], welfare = welfare,
    failures = value@failures[[1L]], n_valid = value@n_valid,
    draw_indices = value@metadata$draw_indices,
    hashes = c(price_pre = safe_hash(value@price_pre[[1L]]),
               price_post = safe_hash(value@price_post[[1L]]),
               shares_post = safe_hash(value@shares_post[[1L]]),
               failures = safe_hash(value@failures[[1L]]),
               draw_indices = safe_hash(value@metadata$draw_indices))
  )
}
metric_row <- function(fixture, policy, repeat_id, run, extraction, fit_load,
                       policy_time, snap, selected_bytes, full_bytes) {
  n <- length(snap$draw_indices)
  valid <- sum(as.integer(snap$n_valid))
  data.frame(
    mode = mode, fixture = fixture, policy = policy, `repeat` = repeat_id,
    status = "complete", error = "", requested = n,
    valid = valid, failure_rate = 1 - valid / n,
    extraction_seconds = extraction, fit_load_seconds = fit_load,
    policy_seconds = policy_time,
    selected_or_full_bytes = selected_bytes, full_bytes = full_bytes,
    price_post_mean = mean(snap$price_post, na.rm = TRUE),
    price_change_mean = mean(snap$price_change, na.rm = TRUE),
    inside_share_change_mean = mean(rowSums(snap$shares_post - snap$shares_pre), na.rm = TRUE),
    cv_mean = mean(snap$welfare[, "cv"], na.rm = TRUE),
    ps_mean = mean(snap$welfare[, "ps"], na.rm = TRUE),
    welfare_mean = mean(snap$welfare[, "welfare"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

library(antitrustBayes)
records <- list()
snapshots <- list()
direct <- list()
row_id <- 0L
for (fixture_id in names(fixtures)) {
  cat("fixture_start=", fixture_id, "\n", sep = "")
  flush.console()
  path <- file.path(base, fixtures[[fixture_id]])
  loaded <- clock(function() readRDS(path))
  if (inherits(loaded$value, "condition")) {
    row_id <- row_id + 1L
    records[[row_id]] <- data.frame(mode = mode, fixture = fixture_id,
      policy = NA_character_, `repeat` = NA_integer_, status = "fit_load_error",
      error = conditionMessage(loaded$value), requested = NA_integer_, valid = 0L,
      failure_rate = 1, extraction_seconds = NA_real_, fit_load_seconds = loaded$elapsed,
      policy_seconds = NA_real_, selected_or_full_bytes = NA_real_, full_bytes = NA_real_,
      price_post_mean = NA_real_, price_change_mean = NA_real_,
      inside_share_change_mean = NA_real_, cv_mean = NA_real_, ps_mean = NA_real_,
      welfare_mean = NA_real_, stringsAsFactors = FALSE)
    next
  }
  fit <- loaded$value
  extracted <- clock(function() {
    if (identical(mode, "full")) rstan::extract(fit@stanfit)
    else antitrustBayes:::.posterior_extract(fit, "simulation")
  })
  if (inherits(extracted$value, "condition")) {
    row_id <- row_id + 1L
    records[[row_id]] <- data.frame(mode = mode, fixture = fixture_id,
      policy = NA_character_, `repeat` = NA_integer_, status = "extraction_error",
      error = conditionMessage(extracted$value), requested = NA_integer_, valid = 0L,
      failure_rate = 1, extraction_seconds = extracted$elapsed, fit_load_seconds = loaded$elapsed,
      policy_seconds = NA_real_, selected_or_full_bytes = NA_real_, full_bytes = NA_real_,
      price_post_mean = NA_real_, price_change_mean = NA_real_,
      inside_share_change_mean = NA_real_, cv_mean = NA_real_, ps_mean = NA_real_,
      welfare_mean = NA_real_, stringsAsFactors = FALSE)
    next
  }
  draws <- extracted$value
  draw_count <- if (is.null(dim(draws$a_market))) length(draws$a_market) else dim(draws$a_market)[[1L]]
  selected_bytes <- as.numeric(object.size(draws))
  full_bytes <- if (identical(mode, "full")) selected_bytes else NA_real_
  direct[[fixture_id]] <- list(mode = mode, extraction_seconds = extracted$elapsed,
                               fit_load_seconds = loaded$elapsed,
                               n_total = draw_count, bytes = selected_bytes,
                               fields = names(draws))
  setup <- tryCatch(antitrustBayes:::.estimated_simulation_units(
    fit, market = "market_1", year = NULL), error = function(e) e)
  if (inherits(setup, "condition")) {
    row_id <- row_id + 1L
    records[[row_id]] <- data.frame(mode = mode, fixture = fixture_id,
      policy = NA_character_, `repeat` = NA_integer_, status = "unitprep_error",
      error = conditionMessage(setup), requested = NA_integer_, valid = 0L,
      failure_rate = 1, extraction_seconds = extracted$elapsed, fit_load_seconds = loaded$elapsed,
      policy_seconds = NA_real_, selected_or_full_bytes = selected_bytes, full_bytes = full_bytes,
      price_post_mean = NA_real_, price_change_mean = NA_real_,
      inside_share_change_mean = NA_real_, cv_mean = NA_real_, ps_mean = NA_real_,
      welfare_mean = NA_real_, stringsAsFactors = FALSE)
    next
  }
  owner <- as.character(fit@data$group)
  d <- fit@data
  selected <- which(as.character(d$market) == "market_1")
  size <- tapply(
    if (identical(as.character(fit@demand_model), "ces")) d$quantities * d$prices else d$quantities,
    as.character(d$group), sum
  )
  largest <- names(size)[order(-as.numeric(size), names(size))][seq_len(min(2L, length(size)))]
  merged <- owner
  if (length(largest) == 2L) merged[merged == largest[[2L]]] <- largest[[1L]]
  policies <- list(
    merger = list(owner = merged, delta = rep(0, nrow(d))),
    cost10 = list(owner = owner, delta = rep(.10, nrow(d))),
    merger_cost10 = list(owner = merged, delta = rep(.10, nrow(d)))
  )
  ## Warm the package and solver once outside the three timed paired repeats.
  set.seed(61000L + match(fixture_id, names(fixtures)))
  invisible(try(suppressWarnings(simulate(
    fit, ownerPost = policies$merger$owner, market = "market_1", n_draws = 5L,
    mcDelta = policies$merger$delta, mc.cores = 1L
  )), silent = TRUE))
  for (policy in names(policies)) {
    for (repeat_id in seq_len(3L)) {
      cat("policy_start=", fixture_id, "/", policy, "/", repeat_id, "\n", sep = "")
      flush.console()
      set.seed(62000L + 100L * match(fixture_id, names(fixtures)) +
        10L * match(policy, names(policies)) + repeat_id)
      run <- clock(function() suppressWarnings(simulate(
        fit, ownerPost = policies[[policy]]$owner, market = "market_1",
        n_draws = 200L, mcDelta = policies[[policy]]$delta, mc.cores = 1L
      )))
      if (inherits(run$value, "condition")) {
        row_id <- row_id + 1L
        records[[row_id]] <- data.frame(mode = mode, fixture = fixture_id,
          policy = policy, `repeat` = repeat_id, status = "simulation_error",
          error = conditionMessage(run$value), requested = 200L, valid = 0L,
          failure_rate = 1, extraction_seconds = extracted$elapsed,
          fit_load_seconds = loaded$elapsed, policy_seconds = run$elapsed,
          selected_or_full_bytes = selected_bytes, full_bytes = full_bytes,
          price_post_mean = NA_real_, price_change_mean = NA_real_,
          inside_share_change_mean = NA_real_, cv_mean = NA_real_, ps_mean = NA_real_,
          welfare_mean = NA_real_, stringsAsFactors = FALSE)
      } else {
        snap <- snapshot(run$value, fit)
        snapshots[[paste(fixture_id, policy, repeat_id, sep = "__")]] <- snap
        row_id <- row_id + 1L
        records[[row_id]] <- metric_row(
          fixture_id, policy, repeat_id, run$value, extracted$elapsed,
          loaded$elapsed, run$elapsed, snap, selected_bytes, full_bytes
        )
      }
    }
  }
  rm(fit, draws, extracted, loaded)
  gc()
}
result <- list(
  mode = mode,
  records = if (length(records)) do.call(rbind, records) else data.frame(),
  snapshots = snapshots, direct = direct
)
saveRDS(result, out_path)
write.csv(result$records, sub("\\.rds$", ".csv", out_path), row.names = FALSE, na = "")
cat("mode=", mode, " rows=", nrow(result$records),
    " complete=", sum(result$records$status == "complete"), "\n", sep = "")
