# Run after the paired benchmark; heavy snapshots stay in the ignored run dir.
run_dir <- Sys.getenv("AUDIT_RUN_DIR", unset =
  "doc/statistical_toolkit_runs/20260916/performance-audit")
old <- readRDS(file.path(run_dir, "baseline-consumer-solve.rds"))
new <- readRDS(file.path(run_dir, "baseline-consumer-observed.rds"))
stopifnot(identical(names(old$snapshots), names(new$snapshots)),
          nrow(old$records) == 27L, nrow(new$records) == 27L)
numeric_difference <- function(a, b) {
  stopifnot(identical(dim(a), dim(b)), identical(length(a), length(b)),
            identical(is.na(a), is.na(b)),
            identical(is.finite(a), is.finite(b)))
  keep <- is.finite(a)
  if (any(keep)) max(abs(a[keep] - b[keep])) else 0
}
fields <- c("price_pre", "price_post", "price_change", "shares_pre",
            "shares_post", "welfare")
checks <- lapply(names(old$snapshots), function(key) {
  a <- old$snapshots[[key]]; b <- new$snapshots[[key]]
  differences <- vapply(fields, function(f) numeric_difference(a[[f]], b[[f]]), 0)
  stopifnot(identical(a$draw_indices, b$draw_indices),
            identical(a$failures, b$failures), identical(a$n_valid, b$n_valid),
            max(differences) <= 1e-6)
  data.frame(snapshot = key, as.list(differences), max_difference = max(differences))
})
checks <- do.call(rbind, checks)
a <- old$records; b <- new$records
keys <- c("fixture", "policy", "repeat.")
stopifnot(identical(a[keys], b[keys]), all(a$status == "complete"),
          all(b$status == "complete"), all(a$valid == 200L), all(b$valid == 200L),
          all(a$baseline_price_calls == 200L), all(b$baseline_price_calls == 0L),
          identical(a$policy_price_calls, b$policy_price_calls),
          identical(a$policy_solver_calls, b$policy_solver_calls),
          identical(a$policy_solver_iterations, b$policy_solver_iterations),
          identical(a$policy_solver_evaluations, b$policy_solver_evaluations))
paired <- data.frame(a[keys], old_wall_seconds = a$policy_seconds,
  new_wall_seconds = b$policy_seconds, speed_ratio = a$policy_seconds / b$policy_seconds,
  old_construction_seconds = a$construction_seconds,
  new_construction_seconds = b$construction_seconds,
  old_baseline_price_calls = a$baseline_price_calls,
  new_baseline_price_calls = b$baseline_price_calls,
  old_baseline_solver_calls = a$baseline_solver_calls,
  new_baseline_solver_calls = b$baseline_solver_calls,
  old_baseline_iterations = a$baseline_solver_iterations,
  new_baseline_iterations = b$baseline_solver_iterations,
  old_baseline_evaluations = a$baseline_solver_evaluations,
  new_baseline_evaluations = b$baseline_solver_evaluations,
  policy_solver_iterations = a$policy_solver_iterations,
  policy_solver_evaluations = a$policy_solver_evaluations,
  max_numerical_difference = checks$max_difference,
  old_failure_rate = a$failure_rate, new_failure_rate = b$failure_rate)
rss <- function(name) {
  lines <- readLines(file.path(run_dir, paste0("baseline-consumer-", name, ".time")))
  line <- lines[grepl("Maximum resident set size", lines, fixed = TRUE)]
  stopifnot(length(line) == 1L)
  as.numeric(sub(".*: *", "", line))
}
aggregate_rows <- lapply(split(seq_len(nrow(paired)),
  interaction(paired$fixture, paired$policy, drop = TRUE)), function(i) {
  data.frame(fixture = paired$fixture[i[1L]], policy = paired$policy[i[1L]],
    repetitions = length(i), posterior_draws_per_repetition = 200L,
    old_wall_seconds = median(paired$old_wall_seconds[i]),
    new_wall_seconds = median(paired$new_wall_seconds[i]),
    speed_ratio = median(paired$old_wall_seconds[i]) / median(paired$new_wall_seconds[i]),
    old_construction_seconds = median(paired$old_construction_seconds[i]),
    new_construction_seconds = median(paired$new_construction_seconds[i]),
    old_process_peak_rss_kib = rss("solve"), new_process_peak_rss_kib = rss("observed"),
    process_peak_memory_ratio = rss("observed") / rss("solve"),
    max_numerical_difference = max(paired$max_numerical_difference[i]),
    old_failures = sum(a$requested[i] - a$valid[i]),
    new_failures = sum(b$requested[i] - b$valid[i]))
})
aggregate_rows <- do.call(rbind, aggregate_rows)
out <- "doc/performance_audit"
write.csv(checks, file.path(out, "baseline_consumer_numerical_gate.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "baseline_consumer_paired.csv"), row.names = FALSE)
write.csv(aggregate_rows, file.path(out, "baseline_consumer_summary.csv"), row.names = FALSE)
print(aggregate_rows, row.names = FALSE)
