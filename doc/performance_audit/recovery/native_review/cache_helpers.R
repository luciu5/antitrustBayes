## Small restart-safe cache used only by the additional native review.

native_review_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."),
                  tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(object, tmp)
  if (!file.rename(tmp, path))
    stop("could not atomically replace ", path)
  invisible(path)
}

native_review_hash <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

native_review_cache <- function(cache_root, specification, run,
                                force = FALSE,
                                success = function(x) identical(x$status, "ok")) {
  stopifnot(is.list(specification), is.function(run))
  cache_hash <- native_review_hash(specification)
  generation <- file.path(cache_root, cache_hash)
  dir.create(generation, recursive = TRUE, showWarnings = FALSE)
  native_review_atomic_save(specification,
                            file.path(generation, "specification.rds"))

  attempt_dirs <- sort(list.dirs(generation, recursive = FALSE, full.names = TRUE))
  for (attempt in attempt_dirs) {
    finished_file <- file.path(attempt, "finished.rds")
    result_file <- file.path(attempt, "result.rds")
    if (!file.exists(finished_file)) {
      if (file.exists(file.path(attempt, "started.rds")) && !isTRUE(force)) {
        return(list(result = list(status = "error",
          case = specification$case,
          error = "unfinished native review attempt requires HELMERT_FORCE=true"),
          cache_hash = cache_hash, generation = generation, attempt = attempt,
          reused = TRUE, state = "INTERRUPTED"))
      }
      next
    }
    state <- readRDS(finished_file)
    if (!identical(state$cache_hash, cache_hash))
      stop("native review cache hash mismatch: ", finished_file)
    if (!file.exists(result_file) || is.null(state$result_sha256) ||
        !identical(digest::digest(file = result_file, algo = "sha256"),
                   state$result_sha256))
      stop("native review cache payload is absent or corrupt: ", result_file)
    result <- readRDS(result_file)
    if (!isTRUE(force) && (identical(state$status, "COMPLETE") ||
                           identical(state$status, "FAILED"))) {
      return(list(result = result, cache_hash = cache_hash,
                  generation = generation, attempt = attempt, reused = TRUE,
                  state = state$status))
    }
  }

  attempt_number <- length(attempt_dirs) + 1L
  attempt <- file.path(generation, sprintf("attempt-%04d", attempt_number))
  dir.create(attempt, recursive = TRUE, showWarnings = FALSE)
  started <- list(cache_hash = cache_hash, specification = specification,
                  status = "RUNNING",
                  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  native_review_atomic_save(started, file.path(attempt, "started.rds"))
  result <- tryCatch(run(), error = function(e)
    list(status = "error", error = conditionMessage(e)))
  ok <- isTRUE(success(result))
  state <- list(cache_hash = cache_hash,
                status = if (ok) "COMPLETE" else "FAILED",
                started_utc = started$started_utc,
                finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  result_file <- file.path(attempt, "result.rds")
  native_review_atomic_save(result, result_file)
  state$result_sha256 <- digest::digest(file = result_file, algo = "sha256")
  native_review_atomic_save(state, file.path(attempt, "finished.rds"))
  list(result = result, cache_hash = cache_hash, generation = generation,
       attempt = attempt, reused = FALSE, state = state$status)
}

native_review_validate_raw <- function(entry, case_name, input_hash,
                                       source_hashes, settings, seeds) {
  if (!is.list(entry) || !identical(entry$schema, 2L) ||
      !identical(entry$role, "candidate") ||
      !identical(entry$case, case_name) ||
      is.null(entry$u) || !is.character(entry$u_hash) ||
      !identical(native_review_hash(entry$u), entry$u_hash) ||
      !is.character(entry$input_hash) || length(entry$input_hash) != 1L ||
      !is.character(entry$kernel_code_hash) ||
      length(entry$kernel_code_hash) != 1L ||
      !identical(entry$input_hash, input_hash) ||
      !identical(entry$source_hashes, source_hashes) ||
      !identical(entry$settings, settings) ||
      !identical(entry$seeds, seeds)) {
    stop("candidate raw-coordinate provenance mismatch for ", case_name)
  }
  invisible(TRUE)
}
