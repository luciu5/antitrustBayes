# Uses installed candidate antitrust and Bayes libraries, unchanged native
# kernels, three warmed repeats and the existing 200-draw policy fixture runner.
# BASELINE_MODE=solve or observed; BENCH_OUT must end in .rds in an ignored dir.
main <- function() {
  library(antitrustBayes)
  baseline_mode <- Sys.getenv("BASELINE_MODE")
  stopifnot(baseline_mode %in% c("solve", "observed"),
            "baseline" %in% names(formals(antitrust::specify)))
  counter <- new.env(parent = emptyenv())
  counter$constructing <- FALSE
  counter$totals <- c(construction_calls = 0, construction_seconds = 0,
    baseline_price_calls = 0, policy_price_calls = 0,
    baseline_solver_calls = 0, baseline_solver_iterations = 0,
    baseline_solver_evaluations = 0, policy_solver_calls = 0,
    policy_solver_iterations = 0, policy_solver_evaluations = 0)
  counter$clock_rows <- list()
  assign(".baseline_consumer_benchmark_counter", counter, .GlobalEnv)
  ns <- asNamespace("antitrustBayes")
  old_realize <- get("as_antitrust_fit", ns)
  replacement <- function(object, spec = NULL, ...) {
    start <- proc.time()[["elapsed"]]
    counter$constructing <- TRUE
    on.exit({
      counter$constructing <- FALSE
      counter$totals["construction_calls"] <- counter$totals["construction_calls"] + 1
      counter$totals["construction_seconds"] <- counter$totals["construction_seconds"] +
        proc.time()[["elapsed"]] - start
    }, add = TRUE)
    old_realize(object, spec = spec, ..., baseline = baseline_mode)
  }
  unlockBinding("as_antitrust_fit", ns)
  assign("as_antitrust_fit", replacement, ns)
  lockBinding("as_antitrust_fit", ns)
  old_solver <- get("nleqslv", asNamespace("nleqslv"))
  assignInNamespace("nleqslv", function(...) {
    answer <- old_solver(...)
    prefix <- if (counter$constructing) "baseline" else "policy"
    keys <- paste0(prefix, c("_solver_calls", "_solver_iterations", "_solver_evaluations"))
    counter$totals[keys] <- counter$totals[keys] +
      c(1, if (is.null(answer$iter)) 0 else answer$iter,
        if (is.null(answer$nfcnt)) 0 else answer$nfcnt)
    answer
  }, "nleqslv")
  suppressMessages(trace("calcPrices", where = asNamespace("antitrust"),
    tracer = quote({
      cc <- get(".baseline_consumer_benchmark_counter", .GlobalEnv)
      key <- if (cc$constructing) "baseline_price_calls" else "policy_price_calls"
      cc$totals[key] <- cc$totals[key] + 1
    }), print = FALSE))
  on.exit({
    suppressMessages(untrace("calcPrices", where = asNamespace("antitrust")))
    assignInNamespace("nleqslv", old_solver, "nleqslv")
    unlockBinding("as_antitrust_fit", ns); assign("as_antitrust_fit", old_realize, ns)
    lockBinding("as_antitrust_fit", ns)
    rm(".baseline_consumer_benchmark_counter", envir = .GlobalEnv)
  }, add = TRUE)
  audit_clock <- function(fun) {
    gc()
    before <- counter$totals
    start <- proc.time()[["elapsed"]]
    value <- tryCatch(fun(), error = function(e) e)
    elapsed <- proc.time()[["elapsed"]] - start
    delta <- counter$totals - before
    if (delta[["construction_calls"]] > 0) {
      counter$clock_rows[[length(counter$clock_rows) + 1L]] <-
        as.data.frame(as.list(delta))
    }
    list(value = value, elapsed = elapsed)
  }
  expressions <- parse("doc/performance_audit/selective_extraction/benchmark-policies-archived.R")
  for (j in seq_along(expressions)) {
    expression <- expressions[[j]]
    if (is.call(expression) && identical(expression[[1L]], as.name("<-"))) {
      if (identical(expression[[2L]], as.name("clock"))) expressions[[j]][[3L]] <- as.name("audit_clock")
      if (identical(expression[[2L]], as.name("fixtures"))) {
        expressions[[j]][[3L]] <- quote(c(
          sim_logit = "logit_bertrand_noisy_271828__D.rds",
          sim_cournot = "logit_cournot_noisy_271828__D.rds",
          sim_ces = "ces_bertrand_noisy_271828__D.rds"))
      }
    }
  }
  eval(expressions, envir = environment())
  stopifnot(length(counter$clock_rows) == nrow(result$records))
  result$records <- cbind(baseline_mode = baseline_mode, result$records,
                         do.call(rbind, counter$clock_rows))
  saveRDS(result, out_path)
  write.csv(result$records, sub("\\.rds$", ".csv", out_path), row.names = FALSE, na = "")
}
main()
