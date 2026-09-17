# Benchmark-only checkpointing. One writer owns each job id; no production
# posterior, solver or sampler options are changed by these helpers.
audit_atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(".checkpoint-", tmpdir = dirname(path))
  saveRDS(value, staging)
  if (!file.rename(staging, path))
    stop("Checkpoint rename failed; staging file retained: ", staging)
  invisible(path)
}

audit_checkpoint_run <- function(id, specification, run, artifact_root,
                                 resume_incomplete = FALSE,
                                 retry_failed = FALSE) {
  stopifnot(is.function(run), is.list(specification))
  if (length(id) != 1L || is.na(id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", id))
    stop("A job id must be one safe directory name.")
  fingerprint <- digest::digest(list(schema = 1L, specification = specification,
    implementation = list(deparse(formals(run), width.cutoff = 500L),
                          deparse(body(run), width.cutoff = 500L))), algo = "sha256")
  location <- file.path(artifact_root, id, fingerprint)
  dir.create(location, recursive = TRUE, showWarnings = FALSE)
  audit_atomic_rds(specification, file.path(location, "specification.rds"))
  attempts <- sort(list.dirs(location, recursive = FALSE, full.names = TRUE))
  outcomes <- list()
  for (attempt in attempts) {
    finished <- file.path(attempt, "finished.rds")
    if (file.exists(finished)) {
      state <- readRDS(finished)
      if (!identical(state$fingerprint, fingerprint))
        stop("Checkpoint fingerprint mismatch: ", finished)
      if (identical(state$status, "COMPLETE")) {
        payload <- file.path(attempt, "result.rds")
        if (!file.exists(payload) || !identical(
            digest::digest(file = payload, algo = "sha256"), state$result_sha256))
          stop("Completed checkpoint payload is absent or corrupt: ", payload)
        state$result <- readRDS(payload)
        state$cached <- TRUE
        return(state)
      }
      outcomes[[length(outcomes) + 1L]] <- state
    } else if (file.exists(file.path(attempt, "started.rds"))) {
      if (!isTRUE(resume_incomplete))
        stop("Unfinished attempt exists; verify its writer has stopped before resuming: ", attempt)
      state <- readRDS(file.path(attempt, "started.rds"))
      if (!identical(state$fingerprint, fingerprint))
        stop("Unfinished checkpoint fingerprint mismatch: ", attempt)
      state$status <- "INTERRUPTED"
      state$observed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
      state$error <- "Previous attempt has no completion record; explicit recovery requested."
      audit_atomic_rds(state, finished)
      outcomes[[length(outcomes) + 1L]] <- state
    }
  }
  failed <- Filter(function(x) identical(x$status, "FAILED"), outcomes)
  if (length(failed) && !isTRUE(retry_failed)) {
    state <- failed[[length(failed)]]
    state$cached <- TRUE
    return(state)
  }
  attempt <- tempfile("attempt-", tmpdir = location)
  if (!dir.create(attempt)) stop("Cannot create attempt directory: ", attempt)
  state <- list(id = id, fingerprint = fingerprint, attempt = attempt,
    status = "RUNNING", started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string, warnings = character(), cached = FALSE)
  audit_atomic_rds(state, file.path(attempt, "started.rds"))
  warnings <- character()
  outcome <- NULL
  elapsed <- system.time(outcome <- withCallingHandlers(
    tryCatch(list(status = "COMPLETE", value = run(), error = NULL),
      error = function(e) list(status = "FAILED", value = NULL,
                               error = conditionMessage(e)),
      interrupt = function(e) list(status = "INTERRUPTED", value = NULL,
                                   error = conditionMessage(e))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))
  state$status <- outcome$status
  state$error <- outcome$error
  state$warnings <- warnings
  state$compute_time <- unclass(elapsed)
  state$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  if (identical(state$status, "COMPLETE")) {
    payload <- file.path(attempt, "result.rds")
    audit_atomic_rds(outcome$value, payload)
    state$result_sha256 <- digest::digest(file = payload, algo = "sha256")
  }
  audit_atomic_rds(state, file.path(attempt, "finished.rds"))
  if (identical(state$status, "COMPLETE")) state$result <- outcome$value
  state
}
