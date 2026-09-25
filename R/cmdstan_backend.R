## CmdStanR owns compilation and inference.  Keep compiled executables outside
## the installed package, which may be on a read-only library path.
.cmdstan_model_cache <- new.env(parent = emptyenv())

.cmdstan_kernel_file <- function(kernel) {
  path <- system.file("stan", paste0(kernel, ".stan"),
                      package = "antitrustBayes")
  if (!nzchar(path)) path <- file.path("inst", "stan", paste0(kernel, ".stan"))
  if (!file.exists(path)) {
    stop("Cannot find the installed Stan kernel: ", kernel, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.cmdstan_model_for_kernel <- function(kernel) {
  source <- .cmdstan_kernel_file(kernel)
  version <- as.character(cmdstanr::cmdstan_version())
  hash <- digest::digest(list(
    source = digest::digest(file = source, algo = "sha256"),
    cmdstan = version,
    cmdstan_path = normalizePath(cmdstanr::cmdstan_path(), mustWork = TRUE),
    cmdstanr = as.character(utils::packageVersion("cmdstanr")),
    platform = R.version$platform
  ), algo = "sha256")
  key <- paste(kernel, hash, sep = "_")
  if (exists(key, envir = .cmdstan_model_cache, inherits = FALSE)) {
    return(get(key, envir = .cmdstan_model_cache, inherits = FALSE))
  }
  base <- getOption("antitrustBayes.cmdstan_cache",
                    tools::R_user_dir("antitrustBayes", which = "cache"))
  cache <- file.path(base, "stan", hash)
  if (!dir.exists(cache) && !dir.create(cache, recursive = TRUE,
                                        showWarnings = FALSE)) {
    stop("Cannot create CmdStan cache directory: ", cache, call. = FALSE)
  }
  cached_source <- file.path(cache, paste0(kernel, ".stan"))
  if (!file.exists(cached_source) && !file.copy(source, cached_source)) {
    stop("Cannot copy Stan kernel into the CmdStan cache.", call. = FALSE)
  }
  model_file <- file.path(cache, paste0(kernel, ".model.rds"))
  model <- if (file.exists(model_file)) tryCatch(readRDS(model_file),
                                                  error = function(e) NULL) else NULL
  if (is.null(model) || !inherits(model, "CmdStanModel") ||
      !file.exists(model$exe_file())) {
    ## Loading an executable alone loses CmdStanR's C++ source needed for
    ## unconstrain_draws()/bridge sampling. Preserve the compiled model object.
    model <- cmdstanr::cmdstan_model(
      cached_source, force_recompile = TRUE,
      cpp_options = list(PRECOMPILED_HEADERS = "false", stan_threads = TRUE)
    )
    saveRDS(model, model_file)
  }
  assign(key, model, envir = .cmdstan_model_cache)
  model
}

.cmdstan_sample_options <- function(options) {
  allowed <- c("chains", "iter", "warmup", "cores", "refresh", "control",
               "seed", "init", "thin", "save_warmup", "sig_figs",
               "output_dir", "threads_per_chain")
  unknown <- setdiff(names(options), allowed)
  if (length(unknown)) {
    stop("Unsupported stan_options for CmdStanR: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  iter <- as.integer(options$iter)
  warmup <- as.integer(options$warmup)
  if (length(iter) != 1L || length(warmup) != 1L || is.na(iter) ||
      is.na(warmup) || warmup < 0L || iter <= warmup) {
    stop("stan_options require iter > warmup >= 0.", call. = FALSE)
  }
  control <- options$control %||% list()
  unsupported_control <- setdiff(names(control),
                                 c("adapt_delta", "max_treedepth", "metric",
                                   "stepsize", "adapt_engaged"))
  if (length(unsupported_control)) {
    stop("Unsupported CmdStanR control option(s): ",
         paste(unsupported_control, collapse = ", "), ".", call. = FALSE)
  }
  result <- list(
    chains = as.integer(options$chains),
    parallel_chains = as.integer(options$cores),
    iter_warmup = warmup,
    iter_sampling = iter - warmup,
    refresh = options$refresh,
    adapt_delta = control$adapt_delta,
    max_treedepth = control$max_treedepth,
    metric = control$metric,
    step_size = control$stepsize,
    adapt_engaged = control$adapt_engaged
  )
  for (field in c("seed", "init", "thin", "save_warmup", "sig_figs",
                  "output_dir", "threads_per_chain")) {
    if (!is.null(options[[field]])) result[[field]] <- options[[field]]
  }
  result[!vapply(result, is.null, logical(1))]
}

.cmdstan_pathfinder_options <- function(options) {
  allowed <- c("seed", "refresh", "init", "output_dir", "draws",
               "num_paths", "single_path_draws", "psis_resample",
               "max_lbfgs_iters", "num_threads")
  unknown <- setdiff(names(options), allowed)
  if (length(unknown)) {
    stop("Unsupported pathfinder_options: ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  if (is.null(options$num_threads)) options$num_threads <- 1L
  options
}

.cmdstan_pathfinder_pareto_k <- function(fit) {
  output <- tryCatch(capture.output(fit$output()),
                     error = function(e) character())
  hit <- regmatches(output, regexec(
    "Pareto k value \\(([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)\\)", output,
    perl = TRUE
  ))
  value <- suppressWarnings(as.numeric(vapply(hit, function(x)
    if (length(x) >= 2L) x[[2L]] else NA_character_, character(1))))
  value <- value[is.finite(value)]
  if (length(value)) value[[length(value)]] else NA_real_
}

.cmdstan_assert_fit <- function(fit, inference, requested_chains = NULL) {
  expected <- if (identical(inference, "hmc")) "CmdStanMCMC" else "CmdStanPathfinder"
  if (!inherits(fit, expected)) {
    stop("CmdStanR returned an unexpected fit type: expected ", expected, ".",
         call. = FALSE)
  }
  codes <- fit$return_codes()
  if (!length(codes) || anyNA(codes) || any(codes != 0L)) {
    stop("CmdStanR did not complete every process successfully.", call. = FALSE)
  }
  if (identical(inference, "hmc") &&
      !identical(as.integer(fit$num_chains()), as.integer(requested_chains))) {
    stop("CmdStanR returned an unexpected chain count.", call. = FALSE)
  }
  draws <- fit$draws(variables = "lp__", format = "draws_array")
  if (length(dim(draws)) != 3L || dim(draws)[[1L]] < 1L ||
      any(!is.finite(as.numeric(draws)))) {
    stop("CmdStanR returned no usable posterior draws.", call. = FALSE)
  }
  invisible(fit)
}

.cmdstan_require_policy_inference <- function(fit, allow_experimental = FALSE) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(invisible(TRUE))
  if (!identical(fit@metadata$inference, "pathfinder")) return(invisible(TRUE))
  if (!isTRUE(allow_experimental)) {
    stop("Pathfinder policy draws are exploratory; set ",
         "allow_experimental = TRUE to use them explicitly.", call. = FALSE)
  }
  warning("Policy results use exploratory Pathfinder draws (Pareto k = ",
          if (is.finite(fit@metadata$pathfinder_pareto_k))
            format(fit@metadata$pathfinder_pareto_k) else "unknown",
          ").", call. = FALSE)
  invisible(TRUE)
}

## Convert CmdStanR's flattened Stan columns to the historical draw-first
## arrays used by policy code.  Assign by Stan indices rather than assuming
## Stan and R have the same matrix flattening order.
.cmdstan_extract <- function(fit, variables = NULL) {
  if (!inherits(fit, "CmdStanFit")) {
    stop("Expected a CmdStanR fit.", call. = FALSE)
  }
  if (is.null(variables)) variables <- fit$metadata()$stan_variables
  variables <- unique(as.character(variables))
  matrix <- as.matrix(fit$draws(variables = variables,
                                format = "draws_matrix"))
  columns <- colnames(matrix)
  result <- vector("list", length(variables))
  names(result) <- variables
  for (variable in variables) {
    prefix <- paste0(variable, "[")
    selected <- which(columns == variable | startsWith(columns, prefix))
    if (!length(selected)) {
      stop("CmdStanR fit is missing posterior field: ", variable, ".",
           call. = FALSE)
    }
    names_selected <- columns[selected]
    if (length(selected) == 1L && identical(names_selected, variable)) {
      result[[variable]] <- as.numeric(matrix[, selected])
      next
    }
    indexes <- lapply(names_selected, function(name) {
      text <- substring(name, nchar(variable) + 2L, nchar(name) - 1L)
      as.integer(strsplit(text, ",", fixed = TRUE)[[1L]])
    })
    rank <- unique(lengths(indexes))
    if (length(rank) != 1L || rank < 1L ||
        any(!is.finite(unlist(indexes)))) {
      stop("Invalid CmdStanR indexed field: ", variable, ".", call. = FALSE)
    }
    dimensions <- vapply(seq_len(rank), function(axis)
      max(vapply(indexes, `[[`, integer(1), axis)), integer(1))
    value <- array(NA_real_, dim = c(nrow(matrix), dimensions))
    for (i in seq_along(selected)) {
      index <- c(list(seq_len(nrow(matrix))), as.list(indexes[[i]]))
      value <- do.call(`[<-`, c(list(value), index,
                               list(value = as.numeric(matrix[, selected[[i]]]))))
    }
    result[[variable]] <- value
  }
  result
}
