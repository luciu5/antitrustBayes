#!/usr/bin/env Rscript

# Base-R aggregation only.  This reads saved benchmark outputs and old CSV
# diagnostics; it never estimates, installs, or modifies source-fit objects.

run_base <- Sys.getenv(
  "TARIFF_RUN_BASE",
  "doc/statistical_toolkit_runs/20260918/tariff-reuse-01"
)
out_dir <- file.path(run_base, "evidence")
old_dir <- file.path(
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality",
  "rotation-rho-sensitivity-02-budget2k"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fixtures <- data.frame(
  fixture = c("logit_bertrand", "rho_panel"),
  seed = c(314159L, 271828L),
  stringsAsFactors = FALSE
)

run_path <- function(fixture, run) {
  file.path(run_base, "posterior-benchmark", fixture, run)
}

read_summary <- function(fixture, run) {
  path <- file.path(run_path(fixture, run), "summary.csv")
  if (!file.exists(path))
    return(data.frame())
  z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  z$fixture_seed <- fixtures$seed[match(fixture, fixtures$fixture)]
  z$run <- run
  z$source_summary <- path
  z
}

# Final rows retain policy-level timings and outcomes. Shared plan,
# construction, and promotion fields are copied only to the merger row for
# each fixture-anchor; they are repeated in the source summaries but are not
# additive across the four policies.
finals <- do.call(rbind, list(
  read_summary("logit_bertrand", "final-200"),
  read_summary("rho_panel", "final-200")
))
if (nrow(finals)) {
  shared <- finals$policy == "merger"
  finals$shared_timing_row <- shared
  finals$shared_plan_seconds <- ifelse(shared, finals$plan_seconds, NA_real_)
  finals$shared_construction_seconds <- ifelse(
    shared, finals$construction_seconds, NA_real_
  )
  finals$shared_promotion_seconds <- ifelse(
    shared, finals$promotion_seconds, NA_real_
  )
  finals$shared_promotion_count <- ifelse(
    shared, finals$promotion_count, NA_real_
  )
  finals$shared_promotion <- ifelse(shared, finals$promotion_shared, NA)
  finals$tariff_consumer_price_fraction <- ifelse(
    finals$policy %in% c("tariff", "combined"), 0.2, 0
  )
  finals$cost_fraction <- ifelse(
    finals$policy %in% c("cost", "combined"), 0.1, 0
  )
  keep <- c(
    "fixture", "fixture_seed", "run", "anchoring", "policy",
    "tariff_consumer_price_fraction", "cost_fraction",
    "selected_draws", "units", "attempts", "usable", "failures",
    "shared_timing_row", "shared_plan_seconds", "shared_construction_seconds",
    "shared_promotion_seconds", "shared_promotion_count", "shared_promotion",
    "solve_seconds", "equivalent_game_seconds", "assembly_seconds",
    "price_error", "share_error", "producer_error", "government_error",
    "physical_ps_delta_mean", "physical_ps_delta_q025",
    "physical_ps_delta_q975", "government_delta_mean",
    "consumer_monetary_change_mean", "price_change_mean",
    "price_change_q025", "price_change_q975", "draw_mean_price_change_q025",
    "draw_mean_price_change_q975", "draw_mean_physical_ps_delta_q025",
    "draw_mean_physical_ps_delta_q975", "source_summary"
  )
  keep <- keep[keep %in% names(finals)]
  write.csv(finals[, keep, drop = FALSE], file.path(out_dir, "benchmark.csv"),
            row.names = FALSE, na = "")
}

# Pair observed and latent-quality policy draws by the structural keys. Price
# is averaged over all products and units per draw. Welfare quantities are
# repeated across products in the source data, so first keep one row per
# uid/index and then average over units per draw.
paired_one <- function(fixture, run, seed) {
  path <- file.path(run_path(fixture, run), "policy-draws.rds")
  if (!file.exists(path)) return(data.frame())
  x <- readRDS(path)
  policies <- intersect(names(x$observed), names(x$latent_quality))
  out <- list()
  max_out <- list()
  for (policy in policies) {
    obs <- x$observed[[policy]]
    lat <- x$latent_quality[[policy]]
    keys <- c("uid", "index", "product")
    if (!all(keys %in% names(obs)) || !all(keys %in% names(lat))) next
    if (!all(c("price_change", "physical_ps_delta", "government_delta",
               "consumer_change") %in% names(obs))) next
    key_obs <- do.call(paste, c(obs[keys], sep = "\r"))
    key_lat <- do.call(paste, c(lat[keys], sep = "\r"))
    if (anyDuplicated(key_obs) || anyDuplicated(key_lat))
      stop("duplicate paired product key in ", fixture, "/", run, "/", policy)
    joined <- merge(
      obs[, c(keys, "price_change", "physical_ps_delta", "government_delta",
              "consumer_change")],
      lat[, c(keys, "price_change", "physical_ps_delta", "government_delta",
              "consumer_change")], by = keys, suffixes = c("_observed", "_latent"),
      sort = FALSE
    )
    if (!nrow(joined)) next
    # Product-level discrepancy is retained separately from draw-level means.
    max_out[[length(max_out) + 1L]] <- data.frame(
      fixture = fixture, fixture_seed = seed, run = run, policy = policy,
      matched_product_rows = nrow(joined),
      max_abs_price_change_discrepancy = max(abs(
        joined$price_change_latent - joined$price_change_observed
      )),
      max_abs_physical_ps_discrepancy = max(abs(
        joined$physical_ps_delta_latent - joined$physical_ps_delta_observed
      )),
      max_abs_government_discrepancy = max(abs(
        joined$government_delta_latent - joined$government_delta_observed
      )),
      max_abs_consumer_change_discrepancy = max(abs(
        joined$consumer_change_latent - joined$consumer_change_observed
      ))
    )
    draw_price_obs <- aggregate(price_change_observed ~ index, joined, mean)
    draw_price_lat <- aggregate(price_change_latent ~ index, joined, mean)
    names(draw_price_obs)[2] <- "price_observed"
    names(draw_price_lat)[2] <- "price_latent"
    draw <- merge(draw_price_obs, draw_price_lat, by = "index", sort = TRUE)
    unit_keys <- !duplicated(joined[, c("uid", "index")])
    unit <- joined[unit_keys, , drop = FALSE]
    draw_w_obs <- aggregate(
      cbind(physical_ps_delta_observed, government_delta_observed,
            consumer_change_observed) ~ index, unit, mean
    )
    draw_w_lat <- aggregate(
      cbind(physical_ps_delta_latent, government_delta_latent,
            consumer_change_latent) ~ index, unit, mean
    )
    draw_w <- merge(draw_w_obs, draw_w_lat, by = "index", sort = TRUE)
    draw <- merge(draw, draw_w, by = "index", sort = TRUE)
    out[[length(out) + 1L]] <- data.frame(
      fixture = fixture, fixture_seed = seed, run = run, policy = policy,
      paired_draws = nrow(draw),
      mean_price_change_latent_minus_observed = draw$price_latent - draw$price_observed,
      mean_physical_ps_delta_latent_minus_observed =
        draw$physical_ps_delta_latent - draw$physical_ps_delta_observed,
      mean_government_delta_latent_minus_observed =
        draw$government_delta_latent - draw$government_delta_observed,
      mean_consumer_change_latent_minus_observed =
        draw$consumer_change_latent - draw$consumer_change_observed
    )
  }
  list(draw = if (length(out)) do.call(rbind, out) else data.frame(),
       max = if (length(max_out)) do.call(rbind, max_out) else data.frame())
}

paired <- list(
  paired_one("logit_bertrand", "final-200", 314159L),
  paired_one("rho_panel", "final-200", 271828L)
)
paired_draws <- do.call(rbind, lapply(paired, `[[`, "draw"))
paired_max <- do.call(rbind, lapply(paired, `[[`, "max"))

if (nrow(paired_draws)) {
  groups <- unique(paired_draws[, c("fixture", "fixture_seed", "run", "policy")])
  pair_rows <- list()
  metrics <- c(
    price_change = "mean_price_change_latent_minus_observed",
    physical_ps_delta = "mean_physical_ps_delta_latent_minus_observed",
    government_delta = "mean_government_delta_latent_minus_observed",
    consumer_change = "mean_consumer_change_latent_minus_observed"
  )
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    z <- paired_draws[paired_draws$fixture == g$fixture &
                      paired_draws$policy == g$policy, , drop = FALSE]
    for (metric in names(metrics)) {
      v <- z[[metrics[[metric]]]]
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        fixture = g$fixture, fixture_seed = g$fixture_seed, run = g$run,
        policy = g$policy, metric = metric, paired_draws = nrow(z),
        mean = mean(v), q025 = unname(quantile(v, .025)),
        q975 = unname(quantile(v, .975)), stringsAsFactors = FALSE
      )
    }
  }
  write.csv(do.call(rbind, pair_rows), file.path(out_dir, "paired_anchor.csv"),
            row.names = FALSE)
  write.csv(paired_max, file.path(out_dir, "paired_anchor_max_product_discrepancy.csv"),
            row.names = FALSE)
}

# Parse /usr/bin/time blocks from the stored benchmark logs.
elapsed_seconds <- function(x) {
  x <- trimws(x)
  p <- strsplit(x, ":", fixed = TRUE)[[1L]]
  if (length(p) == 3L) return(as.numeric(p[1L]) * 3600 + as.numeric(p[2L]) * 60 + as.numeric(p[3L]))
  if (length(p) == 2L) return(as.numeric(p[1L]) * 60 + as.numeric(p[2L]))
  as.numeric(x)
}
time_log <- function(path) {
  if (!file.exists(path)) return(data.frame(
    elapsed_seconds = NA_real_, user_seconds = NA_real_, system_seconds = NA_real_,
    peak_rss_kb = NA_real_, exit_status = NA_real_, log_status = "MISSING"
  ))
  l <- readLines(path, warn = FALSE)
  get_num <- function(pattern) {
    i <- grep(pattern, l)
    if (!length(i)) return(NA_real_)
    as.numeric(sub("^.*: ", "", l[i[[length(i)]]]))
  }
  elapsed <- grep("Elapsed \\(wall clock\\) time", l, value = TRUE)
  elapsed <- if (length(elapsed)) elapsed_seconds(sub("^.*: ", "", elapsed[[length(elapsed)]])) else NA_real_
  data.frame(
    elapsed_seconds = elapsed,
    user_seconds = get_num("^\tUser time"),
    system_seconds = get_num("^\tSystem time"),
    peak_rss_kb = get_num("^\tMaximum resident set size"),
    exit_status = get_num("^\tExit status"),
    log_status = "READ"
  )
}

timing_spec <- data.frame(
  fixture = c("logit_bertrand", "logit_bertrand", "logit_bertrand",
              "rho_panel", "sampled_market"),
  seed = c(314159L, 314159L, 314159L, 271828L, 314159L),
  run = c("final-200", "pilot-5", "pilot-5-cached", "final-200", "sampled-5"),
  draws = c(200L, 5L, 5L, 200L, 5L),
  log = c(
    "logs/posterior-benchmark-logit-bertrand-200.log",
    "logs/posterior-benchmark-pilot-5.log",
    "logs/posterior-benchmark-pilot-5-cached.log",
    "logs/posterior-benchmark-rho-panel-200.log",
    "logs/posterior-benchmark-sampled-5.log"
  ), stringsAsFactors = FALSE
)
timing <- do.call(rbind, lapply(seq_len(nrow(timing_spec)), function(i) {
  z <- cbind(timing_spec[i, , drop = FALSE],
             time_log(file.path(run_base, timing_spec$log[[i]])))
  z
}))
timing$log <- file.path(run_base, timing$log)
timing$timing_scope <- ifelse(timing$run == "sampled-5",
                              "sampled quality validation; 5 draws",
                              ifelse(timing$run %in% c("pilot-5", "pilot-5-cached"),
                                     "pilot; 5 draws", "final fixture; 200 draws"))
write.csv(timing, file.path(out_dir, "timing.csv"), row.names = FALSE)

pilot <- timing[timing$run %in% c("pilot-5", "pilot-5-cached"), , drop = FALSE]
if (nrow(pilot) == 2L) {
  uncached <- pilot$elapsed_seconds[pilot$run == "pilot-5"]
  cached <- pilot$elapsed_seconds[pilot$run == "pilot-5-cached"]
  write.csv(data.frame(
    fixture = "logit_bertrand", seed = 314159L, draws = 5L,
    uncached_elapsed_seconds = uncached, cached_elapsed_seconds = cached,
    uncached_over_cached = uncached / cached,
    claim_scope = "pilot-5 only; no general runtime claim"
  ), file.path(out_dir, "pilot_cache_speed.csv"), row.names = FALSE)
}

# Optional sampled-quality functional check is kept separate from final/HMC
# evidence. Its source summary records 30/30 usable policy attempts and the
# stored time log exits successfully; it is not a posterior-HMC run.
sampled_path <- file.path(run_base, "posterior-benchmark/sampled_market/sampled-5/summary.csv")
if (file.exists(sampled_path)) {
  sampled <- read.csv(sampled_path, stringsAsFactors = FALSE, check.names = FALSE)
  sampled$tariff_consumer_price_fraction <- ifelse(
    sampled$policy %in% c("tariff", "combined"), 0.2, 0
  )
  sampled$cost_fraction <- ifelse(
    sampled$policy %in% c("cost", "combined"), 0.1, 0
  )
  sampled$fixture_seed <- 314159L
  sampled$run_class <- "sampled_quality_5draw_validation"
  sampled$hmc_status <- "NOT_HMC"
  sampled$source_summary <- sampled_path
  write.csv(sampled, file.path(out_dir, "sampled_quality_validation.csv"), row.names = FALSE)
}

# Copy old diagnostics for the two requested source-fit cases with explicit
# provenance. These are preexisting HMC outputs, not new HMC draws.
cases <- c("B_collapsed_market_314159_budget2k", "B_collapsed_market_271828_budget2k")
copy_filtered <- function(input, output, case_col = "case") {
  if (!file.exists(input)) return(0L)
  z <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  if (case_col %in% names(z)) z <- z[z[[case_col]] %in% cases, , drop = FALSE]
  z$provenance <- "PREEXISTING"
  z$source_file <- input
  write.csv(z, file.path(out_dir, output), row.names = FALSE)
  nrow(z)
}
diag_counts <- data.frame(
  diagnostic = c("summary", "chain", "parameter", "coordinate_correlation", "gradient", "warnings"),
  rows = c(
    copy_filtered(file.path(old_dir, "rotation-rho-sensitivity-summary.csv"), "sourcefit_hmc_summary.csv"),
    copy_filtered(file.path(old_dir, "rotation-rho-chain-diagnostics.csv"), "sourcefit_hmc_chain.csv"),
    copy_filtered(file.path(old_dir, "rotation-rho-parameter-diagnostics.csv"), "sourcefit_hmc_parameter.csv"),
    copy_filtered(file.path(old_dir, "rotation-rho-coordinate-correlations.csv"), "sourcefit_hmc_correlations.csv"),
    copy_filtered(file.path(old_dir, "native-target-gradient-check.csv"), "sourcefit_hmc_gradient.csv"),
    copy_filtered(file.path(old_dir, "sampler-warnings.csv"), "sourcefit_hmc_warnings.csv")
  )
)
write.csv(diag_counts, file.path(out_dir, "sourcefit_hmc_diagnostic_counts.csv"), row.names = FALSE)

writeLines(c(
  "Tariff posterior reuse evidence aggregation (2026-09-18)",
  "",
  paste0("Sources are under: ", run_base),
  "Final fixtures: logit_bertrand seed 314159 and rho_panel seed 271828; each has 200 saved posterior draws, 3 markets x 2 periods = 6 units, four policies, and observed/latent_quality anchors.",
  "Scenario columns identify tariff as 0.2 of consumer price and cost as 0.1; combined carries both values.",
  "benchmark.csv keeps policy outcomes row-wise and places shared plan/construction/promotion timing only on the merger row for each fixture-anchor. Do not sum those shared fields across policies.",
  "paired_anchor.csv matches observed and latent_quality by uid/index/product. Price is averaged over products and units per structural draw. Welfare is de-duplicated to uid/index before averaging over units. Reported intervals are quantiles of paired structural-draw means.",
  "paired_anchor_max_product_discrepancy.csv reports the maximum absolute matched product-level discrepancy separately; it is not a posterior interval.",
  "Columns named *_q025/*_q975 in benchmark.csv come from pooled source rows where present; they are pooled product/unit heterogeneity quantiles. draw_mean_* columns are the source draw-level quantiles.",
  "pilot_cache_speed.csv is a 5-draw cache comparison only; its ratio is not a general runtime claim.",
  "sampled_quality_validation.csv is a separately labeled 5-draw sampled-quality functional check and is NOT_HMC. It is not part of the final 200-draw benchmark.",
  "sourcefit_hmc_*.csv files are filtered preexisting diagnostics from rotation-rho-sensitivity-02-budget2k for the two B_collapsed_market cases. Their provenance column is PREEXISTING; no HMC was run by this aggregation.",
  "Missing source files would be recorded as missing timing/row counts; no success is inferred from code or theory."
), file.path(out_dir, "NOTES.txt"))

cat("Wrote evidence to", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("benchmark rows:", if (exists("finals")) nrow(finals) else 0L,
    "paired draw rows:", if (exists("paired_draws")) nrow(paired_draws) else 0L,
    "timing rows:", nrow(timing), "\n")
