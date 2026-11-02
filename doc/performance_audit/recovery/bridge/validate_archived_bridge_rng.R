root <- "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/bridge-archived"
read_result <- function(mode, seed) {
  files <- list.files(file.path(root, paste(mode, seed, sep = "-")),
                      pattern = "finished.rds$", recursive = TRUE, full.names = TRUE)
  states <- lapply(files, readRDS)
  keep <- vapply(states, function(x) identical(x$status, "COMPLETE"), logical(1))
  stopifnot(sum(keep) == 1L)
  state <- states[[which(keep)]]
  payload <- file.path(dirname(files[[which(keep)]]), "result.rds")
  stopifnot(identical(digest::digest(file = payload, algo = "sha256"), state$result_sha256))
  readRDS(payload)
}
rows <- list()
for (seed in c(707101L, 707102L)) {
  old <- read_result("legacy_default", seed)
  new <- read_result("candidate_default", seed)
  serial <- read_result("candidate_seeded_serial", seed)
  fork <- read_result("candidate_seeded_parallel", seed)
  stopifnot(identical(old$rng_after, new$rng_after),
            identical(serial$rng_after, fork$rng_after))
  for (model in c("BF", "BS")) {
    a <- old$comparison[[model]]; b <- new$comparison[[model]]
    stopifnot(identical(names(a), names(b)))
    stopifnot(all(vapply(names(a), function(n) identical(a[[n]], b[[n]]), logical(1))))
    rows[[length(rows)+1L]] <- data.frame(seed = seed, model = model,
      default_public_columns_identical = TRUE, default_rng_identical = TRUE,
      seeded_serial_fork_rng_identical = TRUE,
      original_default_maxiter = 5000L, bridge_repetitions = 3L,
      scope = "separate-model public comparison; no cross-DGP Bayes factor")
  }
}
write.csv(do.call(rbind, rows), "doc/performance_audit/recovery/bridge/archived_rng_gate.csv", row.names = FALSE)
cat("Default columns and RNG gate PASS: four model/seed rows\n")
