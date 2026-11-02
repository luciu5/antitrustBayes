# Use the old matched native library for prepare, the Helmert library for verify.
# These small chains validate saved lp__ restoration, not inference or evidence.
library(antitrustBayes)
library(rstan)
mode <- Sys.getenv("MIGRATION_GATE_MODE")
out <- Sys.getenv("MIGRATION_GATE_DIR")
fixture_source <- Sys.getenv("MIGRATION_FIXTURE_SOURCE")
stopifnot(mode %in% c("prepare", "verify"), nzchar(out), nzchar(fixture_source))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(fixture_source)
ns <- asNamespace("antitrustBayes")
models <- get("stanmodels", ns)
e <- new.env(parent = globalenv())
for (n in ls(ns, all.names = TRUE)) assign(n, get(n, ns), e)
sys.source("R/bma.R", e)
rows <- list()
for (demand in c("logit", "ces")) for (structure in c("simultaneous", "stackelberg")) {
  for (coordinates in c("legacy", "rotated")) {
    key <- paste(demand, structure, coordinates, sep = "_")
    path <- file.path(out, paste0(key, ".rds"))
    cat(mode, key, "\n"); flush.console()
    if (mode == "prepare") {
      dat <- make_data(demand = demand, structure = structure,
        market_s0_parameterization = coordinates, use_s0_market_year_pool = TRUE,
        use_centered_s0_pool = FALSE,
        k_support = if (structure == "stackelberg") 1:3 else NULL,
        k_prior = if (structure == "stackelberg") rep(1/3, 3) else NULL)
      kernel <- paste0(demand, if (structure == "simultaneous") "_supply" else "_stackelberg")
      warnings <- character()
      sf <- withCallingHandlers(sampling(models[[kernel]], data = dat,
        chains = 1L, iter = 100L, warmup = 50L, seed = 93017L, init = 0,
        refresh = 0L, control = list(adapt_delta = .95, max_treedepth = 13L)),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        })
      fit <- methods::new("BayesLogit", stanfit = sf, stan_data = dat)
      sp <- get_sampler_params(sf, inc_warmup = FALSE)[[1L]]
      metadata <- list(kernel = kernel, coordinates = coordinates,
        divergences = sum(sp[, "divergent__"]), depth_hits = sum(sp[, "treedepth__"] >= 13),
        warnings = warnings)
      saveRDS(list(fit = fit, metadata = metadata), path)
      rows[[key]] <- data.frame(case = key, kernel = kernel,
        flag = dat$market_s0_parameterization, saved_rows = nrow(sp),
        divergences = metadata$divergences, depth_hits = metadata$depth_hits,
        warning_count = length(warnings), scope = "normalization fixture only; no inference claim")
    } else {
      saved <- readRDS(path); fit <- saved$fit
      before <- serialize(fit, NULL); set.seed(811L)
      before_rng <- .Random.seed
      ctx <- e$.bayes_bridge_restore_context(fit)
      stopifnot(identical(before_rng, .Random.seed),
                identical(serialize(fit, NULL), before), isTRUE(ctx$migration))
      checks <- ctx$target_check
      contract <- e$.bayes_bridge_code_contract(fit@stanfit)
      private <- e$.bayes_bridge_private_data(fit, contract)
      stopifnot(identical(private$market_s0_parameterization,
                          fit@stan_data$market_s0_parameterization),
                identical(contract$migration_geometry, "preserve"))
      rows[[key]] <- data.frame(case = key, kernel = ctx$kernel,
        saved_flag = fit@stan_data$market_s0_parameterization,
        restored_flag = private$market_s0_parameterization,
        checked_rows = nrow(checks), max_normalized_target_difference = max(checks$abs_difference),
        raw_schema = "PASS", rng_preserved = TRUE, fit_immutable = TRUE,
        divergences_in_saved_fixture = saved$metadata$divergences,
        depth_hits_in_saved_fixture = saved$metadata$depth_hits,
        scope = "saved normalized target gate; no inference or bridge precision claim")
    }
  }
}
write.csv(do.call(rbind, rows), file.path(out, paste0(mode, ".csv")), row.names = FALSE)
