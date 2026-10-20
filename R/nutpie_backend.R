## Native nutpieR uses the same installed Stan source through BridgeStan.
## Native model handles never enter the saved fit; arrays and provenance do.
.nutpie_require <- function() {
  if (!requireNamespace("nutpieR", quietly = TRUE) ||
      utils::packageVersion("nutpieR") < "1.9.1") {
    stop("backend = 'nutpie' requires optional package 'nutpieR' >= 1.9.1. ",
         "Install a matching release from https://github.com/andytimm/nutpieR ",
         "or https://andytimm.r-universe.dev.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.nutpie_cache_base <- function() {
  getOption("antitrustBayes.nutpie_cache",
            tools::R_user_dir("antitrustBayes", which = "cache"))
}

.nutpie_sample_native <- function(...) nutpieR::nutpie_sample(...)
.nutpie_diagnostics_native <- function(draws) nutpieR::nutpie_diagnostics(draws)
.nutpie_compile_native <- function(...) nutpieR::nutpie_compile_model(...)

## nutpieR's public cache root follows R_USER_CACHE_DIR. Scope it to the same
## explicit backend cache, restoring the caller's setting even after errors.
.nutpie_with_cache <- function(expr) {
  base <- file.path(.nutpie_cache_base(), "nutpie-native")
  if (!dir.exists(base) && !dir.create(base, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create native Nutpie cache: ", base, call. = FALSE)
  }
  old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else
    Sys.setenv(R_USER_CACHE_DIR = old), add = TRUE)
  Sys.setenv(R_USER_CACHE_DIR = normalizePath(base, mustWork = TRUE))
  force(expr)
}

## This deliberately accepts only the declaration grammar used in the four
## installed kernels. Unknown declarations fail instead of changing data hashes.
.nutpie_stan_data_fields <- function(kernel) {
  code <- paste(readLines(.cmdstan_kernel_file(kernel), warn = FALSE), collapse = "\n")
  code <- gsub("/\\*.*?\\*/", "", code, perl = TRUE)
  code <- gsub("//[^\n]*", "", code, perl = TRUE)
  start <- regexpr("(^|\\n)[[:space:]]*data[[:space:]]*\\{", code, perl = TRUE)
  if (start < 0L) stop("Cannot identify Stan data block.", call. = FALSE)
  body <- substring(code, start + attr(start, "match.length"))
  end <- regexpr("}", body, fixed = TRUE)
  if (end < 0L) stop("Cannot identify Stan data declarations.", call. = FALSE)
  declarations <- trimws(strsplit(substring(body, 1L, end - 1L), ";", fixed = TRUE)[[1L]])
  declarations <- declarations[nzchar(declarations)]
  if (any(!grepl("^(array\\s*\\[.*\\]\\s*)?(int|real|vector|matrix)\\b", declarations,
                 perl = TRUE)) ||
      any(!grepl("[[:alpha:]_][[:alnum:]_]*$", declarations))) {
    stop("Unsupported Stan data declaration; use compiler metadata.", call. = FALSE)
  }
  fields <- sub(".*?([[:alpha:]_][[:alnum:]_]*)$", "\\1", declarations, perl = TRUE)
  if (anyDuplicated(fields)) stop("Duplicate Stan data fields.", call. = FALSE)
  fields
}

.nutpie_model_for_kernel <- function(kernel) {
  .nutpie_require()
  source <- .cmdstan_kernel_file(kernel)
  source_hash <- digest::digest(file = source, algo = "sha256")
  version <- as.character(utils::packageVersion("nutpieR"))
  key <- digest::digest(list(source_hash, version, R.version$platform), algo = "sha256")
  ## Always let nutpieR check its content/BridgeStan/compiler cache. Retaining
  ## an in-memory entry alone could hide a changed compiler configuration.
  base <- .nutpie_cache_base()
  cache <- file.path(base, "nutpie", key)
  if (!dir.exists(cache) && !dir.create(cache, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create Nutpie source cache: ", cache, call. = FALSE)
  }
  staged <- file.path(cache, basename(source))
  if (!file.exists(staged) ||
      !identical(digest::digest(file = staged, algo = "sha256"), source_hash)) {
    if (!file.copy(source, staged, overwrite = TRUE)) {
      stop("Cannot stage unchanged Stan source for Nutpie.", call. = FALSE)
    }
  }
  started <- proc.time()[["elapsed"]]
  native <- .nutpie_with_cache(.nutpie_compile_native(stan_file = staged,
    compile_args = "STAN_THREADS=true"))
  elapsed <- proc.time()[["elapsed"]] - started
  lib <- normalizePath(native$lib_path, mustWork = TRUE)
  metadata <- list(
    kernel = kernel, source_code = paste(readLines(source, warn = FALSE), collapse = "\n"),
    source_hash = source_hash, source_file = source,
    lib_path = lib, lib_hash = digest::digest(file = lib, algo = "sha256"),
    nutpie_version = version,
    nutpie_remote_sha = utils::packageDescription("nutpieR")$RemoteSha %||% NA_character_,
    nuts_rs_version = NA_character_, bridgestan_version = NA_character_,
    engine_version_limitations = paste(
      "nutpieR's public API does not expose the linked nuts-rs and BridgeStan versions;",
      "the installed nutpieR version, source revision when available, and library hash are retained."),
    compile_time = elapsed,
    compile_args = "STAN_THREADS=true", stanc_args = character(),
    platform = R.version$platform
  )
  structure(list(native_model = native, metadata = metadata), class = "NutpieModel")
}

.nutpie_count <- function(value, label, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum || value > .Machine$integer.max ||
      value != floor(value)) {
    stop(label, " must be a finite integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.nutpie_sample_options <- function(options, nutpie_options = list()) {
  named_list <- function(value) is.list(value) && (!length(value) ||
    (!is.null(names(value)) && !anyNA(names(value)) &&
     all(nzchar(names(value))) && !anyDuplicated(names(value))))
  if (!named_list(options)) stop("stan_options must be a uniquely named list.", call. = FALSE)
  allowed <- c("chains", "iter", "warmup", "cores", "refresh", "control",
               "seed", "init", "thin", "save_warmup", "threads_per_chain")
  unknown <- setdiff(names(options), allowed)
  if (length(unknown)) stop("Unsupported stan_options for Nutpie: ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  control <- options$control %||% list()
  if (!named_list(control)) stop("Nutpie control must be a uniquely named list.", call. = FALSE)
  unknown <- setdiff(names(control), c("adapt_delta", "max_treedepth"))
  if (length(unknown)) stop("Unsupported Nutpie control option(s): ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  for (name in c("thin", "threads_per_chain")) {
    if (!is.null(options[[name]]) && !identical(.nutpie_count(options[[name]], name), 1L)) {
      stop("Nutpie supports only ", name, " = 1.", call. = FALSE)
    }
  }
  iter <- .nutpie_count(options$iter, "stan_options$iter")
  warmup <- .nutpie_count(options$warmup, "stan_options$warmup")
  if (iter <= warmup) stop("stan_options require iter > warmup > 0 for Nutpie.", call. = FALSE)
  result <- list(num_chains = .nutpie_count(options$chains, "stan_options$chains"),
    cores = .nutpie_count(options$cores, "stan_options$cores"),
    num_draws = iter - warmup, num_warmup = warmup,
    refresh = .nutpie_count(options$refresh %||% 100L, "stan_options$refresh", 0L),
    target_accept = control$adapt_delta, max_treedepth = control$max_treedepth)
  for (name in c("seed", "init", "save_warmup")) {
    if (!is.null(options[[name]])) result[[name]] <- options[[name]]
  }
  if (!is.null(result$seed)) result$seed <- .nutpie_count(result$seed, "stan_options$seed", 0L)
  if (!named_list(nutpie_options)) {
    stop("nutpie_options must be a uniquely named list.", call. = FALSE)
  }
  extra_allowed <- c("adaptation", "mindepth", "max_energy_error", "extra_doublings",
                     "mass_matrix_gamma", "mass_matrix_eigval_cutoff", "progress",
                     "store_divergences", "store_mass_matrix", "store_unconstrained",
                     "store_gradient")
  unknown <- setdiff(names(nutpie_options), extra_allowed)
  if (length(unknown)) stop("Unsupported nutpie_options: ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  if (is.null(nutpie_options$adaptation)) nutpie_options$adaptation <- "diag"
  if (!is.character(nutpie_options$adaptation) || length(nutpie_options$adaptation) != 1L ||
      !nutpie_options$adaptation %in% c("diag", "low_rank")) {
    stop("nutpie_options$adaptation must be 'diag' or 'low_rank'.", call. = FALSE)
  }
  if (nutpie_options$adaptation == "diag" &&
      any(c("mass_matrix_gamma", "mass_matrix_eigval_cutoff") %in% names(nutpie_options))) {
    stop("Nutpie low-rank mass-matrix controls require adaptation = 'low_rank'.", call. = FALSE)
  }
  if (!is.null(result$target_accept) &&
      (!is.numeric(result$target_accept) || length(result$target_accept) != 1L ||
       !is.finite(result$target_accept) || result$target_accept <= 0 || result$target_accept >= 1)) {
    stop("control$adapt_delta must be a finite number in (0, 1).", call. = FALSE)
  }
  if (!is.null(result$max_treedepth)) result$max_treedepth <-
    .nutpie_count(result$max_treedepth, "control$max_treedepth")
  for (name in intersect(c("mindepth", "extra_doublings"), names(nutpie_options))) {
    nutpie_options[[name]] <- .nutpie_count(nutpie_options[[name]], name, 0L)
  }
  if (!is.null(nutpie_options$mindepth) && !is.null(result$max_treedepth) &&
      nutpie_options$mindepth > result$max_treedepth) {
    stop("Nutpie mindepth cannot exceed control$max_treedepth.", call. = FALSE)
  }
  for (name in intersect(c("max_energy_error", "mass_matrix_gamma", "mass_matrix_eigval_cutoff"),
                         names(nutpie_options))) {
    value <- nutpie_options[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be a positive finite number.", call. = FALSE)
    }
  }
  for (name in intersect(c("store_divergences", "store_mass_matrix", "store_unconstrained",
                           "store_gradient"), names(nutpie_options))) {
    value <- nutpie_options[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.null(result$save_warmup) &&
      (!is.logical(result$save_warmup) || length(result$save_warmup) != 1L || is.na(result$save_warmup))) {
    stop("save_warmup must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(nutpie_options$progress) &&
      (!is.character(nutpie_options$progress) || length(nutpie_options$progress) != 1L ||
       !nutpie_options$progress %in% c("auto", "cli", "text", "none"))) {
    stop("progress must be 'auto', 'cli', 'text', or 'none'.", call. = FALSE)
  }
  c(result[!vapply(result, is.null, logical(1))], nutpie_options)
}

.nutpie_diagnostics_array <- function(diagnostics, num_draws, num_chains) {
  index_fields <- c("chain", "draw")
  total <- num_draws * num_chains
  if (!is.list(diagnostics) || !all(index_fields %in% names(diagnostics)) ||
      length(diagnostics$chain) != total || length(diagnostics$draw) != total) {
    stop("Nutpie returned incomplete diagnostic chain/draw indices.", call. = FALSE)
  }
  chain <- diagnostics$chain
  draw <- diagnostics$draw
  if (!is.numeric(chain) || !is.numeric(draw) ||
      any(!is.finite(chain)) || any(!is.finite(draw)) || any(chain != floor(chain)) ||
      any(draw != floor(draw)) || any(chain < 1L | chain > num_chains) ||
      any(draw < 1L | draw > num_draws)) {
    stop("Nutpie returned invalid diagnostic chain/draw indices.", call. = FALSE)
  }
  index <- draw + (chain - 1L) * num_draws
  if (anyDuplicated(index)) stop("Nutpie returned duplicate diagnostic draws.", call. = FALSE)
  aliases <- c(lp__ = "logp", divergent__ = "diverging", treedepth__ = "depth",
    n_leapfrog__ = "n_steps", energy__ = "energy", accept_stat__ = "mean_tree_accept",
    stepsize__ = "step_size", maxdepth_reached__ = "maxdepth_reached")
  present <- aliases %in% names(diagnostics)
  aliases <- aliases[present]
  result <- array(NA_real_, c(num_draws, num_chains, length(aliases)),
    dimnames = list(iteration = seq_len(num_draws), chain = seq_len(num_chains),
                    variable = names(aliases)))
  for (i in seq_along(aliases)) {
    value <- diagnostics[[aliases[[i]]]]
    if (!(is.numeric(value) || is.logical(value)) || length(value) != total) {
      stop("Invalid Nutpie diagnostic field: ", aliases[[i]], ".", call. = FALSE)
    }
    reordered <- numeric(total)
    reordered[index] <- as.numeric(value)
    result[, , i] <- matrix(reordered, num_draws, num_chains)
  }
  result
}

.nutpie_sample <- function(model, data, common_options, nutpie_options = list()) {
  .nutpie_require()
  options <- .nutpie_sample_options(common_options, nutpie_options)
  if (is.null(options$seed)) options$seed <- sample.int(.Machine$integer.max, 1L)
  if (!inherits(model, "NutpieModel")) stop("Expected a compiled NutpieModel.", call. = FALSE)
  json <- tempfile("antitrustBayes-nutpie-data-", fileext = ".json")
  on.exit(unlink(json), add = TRUE)
  cmdstanr::write_stan_json(data, json)
  data_json <- paste(readLines(json, warn = FALSE), collapse = "\n")
  old_threads <- Sys.getenv("STAN_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old_threads)) Sys.unsetenv("STAN_NUM_THREADS") else
    Sys.setenv(STAN_NUM_THREADS = old_threads), add = TRUE)
  Sys.setenv(STAN_NUM_THREADS = "1")
  draws <- .nutpie_with_cache(do.call(.nutpie_sample_native,
                   c(list(model = model$native_model, data = data_json), options)))
  shape <- dim(draws)
  if (!is.numeric(draws) || length(shape) != 3L || shape[[3L]] < 1L ||
      !identical(as.integer(shape[1:2]), c(options$num_draws, options$num_chains))) {
    stop("Nutpie did not return every requested chain and posterior draw.", call. = FALSE)
  }
  fields <- dimnames(draws)[[3L]]
  if (is.null(fields) || anyNA(fields) || any(!nzchar(fields)) ||
      anyDuplicated(fields)) {
    stop("Nutpie returned missing or duplicate posterior field names.", call. = FALSE)
  }
  raw_diagnostics <- .nutpie_diagnostics_native(draws)
  ## Drop the print helper's reference to draws: it duplicates all attributes
  ## in serialized fits and is unnecessary for numerical diagnostics.
  attr(raw_diagnostics, "draws") <- NULL
  diagnostics <- .nutpie_diagnostics_array(raw_diagnostics,
    options$num_draws, options$num_chains)
  lp <- match("lp__", dimnames(diagnostics)[[3L]])
  if (is.na(lp)) stop("Nutpie did not return log-density diagnostics.", call. = FALSE)
  if (!"lp__" %in% fields) {
    combined <- array(NA_real_, c(dim(draws)[1:2], dim(draws)[[3L]] + 1L),
      dimnames = list(iteration = dimnames(draws)[[1L]], chain = dimnames(draws)[[2L]],
        variable = c(dimnames(draws)[[3L]], "lp__")))
    combined[, , seq_len(dim(draws)[[3L]])] <- draws
    combined[, , dim(combined)[[3L]]] <- diagnostics[, , lp]
    combined <- posterior::as_draws_array(combined)
  } else {
    if (!identical(as.numeric(draws[, , "lp__"]), as.numeric(diagnostics[, , lp]))) {
      stop("Nutpie posterior and diagnostic log densities disagree.", call. = FALSE)
    }
    ## Keep one copy of native diagnostics and warmup arrays in the fit.
    combined <- posterior::as_draws_array(array(as.numeric(draws), shape,
                                               dimnames = dimnames(draws)))
  }
  metadata <- c(model$metadata, list(
    stan_variables = unique(sub("\\[.*$", "", fields[!grepl("__$", fields)])),
    inference = "hmc", backend = "nutpie", num_chains = options$num_chains,
    threads_per_chain = 1L,
    num_samples = options$num_draws, num_warmup = options$num_warmup,
    stan_data_hash = digest::digest(data, algo = "sha256"),
    data_json = data_json, data_json_hash = digest::digest(data_json, algo = "sha256"),
    seed = options$seed, resolved_options = options,
    sampler_config = attr(draws, "sampler_config", exact = TRUE),
    sampling_time_native = attr(draws, "sampling_time", exact = TRUE),
    warmup_time = NA_real_, post_warmup_sampling_time = NA_real_,
    sampling_time_scope = paste("Native warmup and sampling together, including callbacks and",
      "native extraction; excludes R setup/assembly. Separate phase timings unavailable."),
    provenance = list(nutpieR = model$metadata$nutpie_version,
      nutpie_remote_sha = model$metadata$nutpie_remote_sha,
      nuts_rs = model$metadata$nuts_rs_version,
      bridgestan = model$metadata$bridgestan_version,
      posterior = as.character(utils::packageVersion("posterior")),
      R = as.character(getRversion()), platform = R.version$platform)))
  fit <- structure(list(draws_array = combined, diagnostics_array = diagnostics,
    raw_diagnostics = raw_diagnostics,
    warmup_draws = attr(draws, "warmup_draws", exact = TRUE),
    warmup_diagnostics = attr(draws, "warmup_diagnostics", exact = TRUE),
    metadata = metadata), class = "NutpieFit")
  .nutpie_assert_fit(fit, "hmc", options$num_chains)
  fit
}

.nutpie_assert_fit <- function(fit, inference, requested_chains = NULL) {
  if (!inherits(fit, "NutpieFit") || !identical(inference, "hmc")) {
    stop("Nutpie supports HMC fits only.", call. = FALSE)
  }
  draws <- fit$draws_array
  dims <- dim(draws)
  requested_chains <- requested_chains %||% fit$metadata$num_chains
  if (!is.numeric(draws) || length(dims) != 3L || any(dims < 1L) ||
      !identical(as.integer(dims[[2L]]), as.integer(requested_chains)) ||
      !identical(as.integer(dims[[1L]]), as.integer(fit$metadata$num_samples))) {
    stop("Nutpie did not return every requested chain and posterior draw.", call. = FALSE)
  }
  fields <- dimnames(draws)[[3L]]
  if (is.null(fields) || anyNA(fields) || any(!nzchar(fields)) ||
      anyDuplicated(fields)) {
    stop("Nutpie returned missing or duplicate posterior field names.", call. = FALSE)
  }
  lp <- match("lp__", fields)
  if (is.na(lp) || any(!is.finite(draws[, , lp]))) {
    stop("Nutpie returned no usable posterior log densities.", call. = FALSE)
  }
  if (anyNA(draws)) stop("Nutpie returned unavailable parameter or generated-quantity draws.", call. = FALSE)
  ## The frozen Stackelberg kernels intentionally emit -Inf component_log_score
  ## for zero-prior partitions. A nonfinite structural coordinate or a failed
  ## positive-prior component must not be accepted as a completed usable fit.
  nonfinite <- which(!is.finite(draws), arr.ind = TRUE)
  if (nrow(nonfinite)) {
    for (i in seq_len(nrow(nonfinite))) {
      at <- nonfinite[i, ]
      field <- fields[[at[[3L]]]]
      permitted <- fit$metadata$kernel %in% c("logit_stackelberg", "ces_stackelberg") &&
        grepl("^component_log_score\\[[1-9][0-9]*\\]$", field) &&
        identical(as.numeric(draws[at[[1L]], at[[2L]], at[[3L]]]), -Inf)
      prior <- match(sub("^component_log_score", "prior_component_prob", field), fields)
      if (!isTRUE(permitted) || is.na(prior) ||
          !isTRUE(draws[at[[1L]], at[[2L]], prior] == 0)) {
        stop("Nutpie returned nonfinite posterior field: ", field, ".", call. = FALSE)
      }
    }
  }
  invisible(fit)
}
