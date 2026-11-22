## CmdStanR owns compilation and inference.  Keep compiled executables outside
## the installed package, which may be on a read-only library path.
.cmdstan_model_cache <- new.env(parent = emptyenv())

.cmdstan_kernel_file <- function(kernel, model_spec_version = NULL, stan_root = NULL) {
  legacy <- !is.null(model_spec_version) && identical(as.integer(model_spec_version), 14L)
  if (!legacy && (endsWith(kernel, "_relative_panel") || kernel == "logit_ordered_panel")) {
    base_kernel <- if (kernel == "logit_ordered_panel") "logit_stackelberg" else
      sub("_relative_panel$", "", kernel)
    base <- .cmdstan_kernel_file(base_kernel, stan_root = stan_root)
    source <- if (kernel == "logit_ordered_panel") .balanced_panel_kernel_file(base) else base
    return(.relative_panel_kernel_file(kernel, source, file.path(dirname(base), "include")))
  }
  if (!legacy && identical(kernel, "logit_stackelberg_panel"))
    return(.balanced_panel_kernel_file(.cmdstan_kernel_file("logit_stackelberg", stan_root = stan_root)))
  filename <- paste0(kernel, ".stan")
  if (legacy) filename <- file.path("legacy14", filename)
  if (is.null(stan_root)) {
    path <- system.file("stan", filename, package = "antitrustBayes")
    if (!nzchar(path)) path <- file.path("inst", "stan", filename)
  } else {
    path <- file.path(stan_root, filename)
  }
  if (!file.exists(path)) {
    stop("Cannot find the installed Stan kernel: ", kernel, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

# Compose a separately identified kernel, sharing the exact declarations, priors
# and generated quantities of the generic kernel. Its original source is also
# the evidence reference for saved fits, so it stays byte-for-byte intact.
.balanced_panel_kernel_file <- function(base_source) {
  base <- readLines(base_source, warn = FALSE)
  likelihood <- readLines(file.path(dirname(base_source), "include",
    "balanced_panel_likelihood.stan"), warn = FALSE)
  transformed_data <- which(base == "transformed data {")
  likelihood_start <- which(trimws(base) == "// The share likelihood is common across all ordered-k states.") - 1L
  quantities <- which(base == "generated quantities {")
  if (length(transformed_data) != 1L || length(likelihood_start) != 1L || length(quantities) != 1L ||
      !identical(trimws(base[likelihood_start]), "if (share_basis == 2) {") ||
      transformed_data >= likelihood_start || likelihood_start >= quantities)
    stop("Cannot locate the ordered-k kernel's shared blocks.")
  code <- c(base[seq_len(transformed_data)], "  array[N_firm_mky] int panel_firm_obs;",
    "  for (i in 1:N) panel_firm_obs[obs_to_firm_mky[i]] = i;",
    base[(transformed_data + 1L):(likelihood_start - 1L)], likelihood, "}", "",
    base[quantities:length(base)])
  bytes <- charToRaw(paste0(paste(code, collapse = "\n"), "\n"))
  hash <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  path <- file.path(tempdir(), paste0("logit_stackelberg_panel_", hash, ".stan"))
  if (!file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), hash))
    writeBin(bytes, path)
  path
}

# Keep previously sealed sources available. Panel calculations use the original local or
# global quality coordinates and preserve the period-only reporting contract.
.relative_panel_kernel_file <- function(kernel, source, includes) {
  lines <- readLines(source, warn = FALSE)
  insert <- function(anchor, addition, after = FALSE) {
    at <- which(lines == anchor)
    if (length(at) != 1L) stop("Missing or ambiguous panel kernel anchor: ", anchor)
    lines <<- append(lines, addition, after = at - as.integer(!after))
  }
  insert("  real bounded_inv_logit(real eta) {",
         readLines(file.path(includes, "relative_panel_functions.stan")))
  insert("transformed data {", c(
    "  array[N_market] int market_panel;", "  int global_panel = 0;"), after = TRUE)
  end <- max(which(lines == "}" & seq_along(lines) < match("parameters {", lines)))
  lines <- append(lines, readLines(file.path(includes, "relative_panel_design.stan")), after = end - 1L)
  lines <- sub("relative_local_gaussian_logpdf(relative_base_resid, relative_base_var,",
    "panel_local_gaussian_logpdf(market_panel[m] > 0, relative_base_resid, relative_base_var,", lines, fixed = TRUE)
  lines <- sub("relative_local_gaussian_logpdf(relative_cond_resid, relative_cond_var,",
    "panel_local_gaussian_logpdf(market_panel[m] == 2, relative_cond_resid, relative_cond_var,", lines, fixed = TRUE)
  lines <- sub("score[kk] += relative_gaussian_logpdf(",
    "score[kk] += panel_gaussian_logpdf(market_panel[m] == 2,", lines, fixed = TRUE)
  lines <- sub("relative_share_ll = relative_gaussian_logpdf(",
    "relative_share_ll = panel_gaussian_logpdf(global_panel > 0,", lines, fixed = TRUE)
  lines <- sub("relative_margin_ll + relative_gaussian_logpdf(",
    "relative_margin_ll + panel_gaussian_logpdf(global_panel == 2,", lines, fixed = TRUE)
  bytes <- charToRaw(paste0(paste(lines, collapse = "\n"), "\n"))
  hash <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  path <- file.path(tempdir(), paste0(kernel, "_", hash, ".stan"))
  if (!file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), hash))
    writeBin(bytes, path)
  path
}

.cmdstan_model_for_kernel <- function(kernel, model_spec_version = NULL) {
  source <- if (is.null(model_spec_version)) .cmdstan_kernel_file(kernel) else
    .cmdstan_kernel_file(kernel, model_spec_version)
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
  ## Concurrent fits can race to create the same directory or file; another
  ## process winning is success, so test existence after the attempt.
  if (!dir.exists(cache)) {
    suppressWarnings(dir.create(cache, recursive = TRUE, showWarnings = FALSE))
    if (!dir.exists(cache)) {
      stop("Cannot create CmdStan cache directory: ", cache, call. = FALSE)
    }
  }
  cached_source <- file.path(cache, paste0(kernel, ".stan"))
  if (!file.exists(cached_source)) {
    ## Rename is atomic, so a concurrent fit never compiles a half-written file.
    staging <- tempfile(paste0(kernel, "_"), tmpdir = cache, fileext = ".stan")
    copied <- suppressWarnings(file.copy(source, staging))
    if (copied) suppressWarnings(file.rename(staging, cached_source))
    unlink(staging)
    if (!file.exists(cached_source)) {
      stop("Cannot copy Stan kernel into the CmdStan cache.", call. = FALSE)
    }
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
    message <- "CmdStanR did not complete every process successfully."
    failed <- which(is.na(codes) | codes != 0L)
    if (length(failed)) {
      message <- paste0(message, " Failed process return codes: ",
        paste(paste0(failed, "=", codes[failed]), collapse = ", "), ".")
      if (is.function(fit$output)) {
        output <- vapply(failed, function(id) {
          lines <- tryCatch(capture.output(fit$output(id = id)),
                            error = function(e) character())
          lines <- lines[nzchar(trimws(lines))]
          if (!length(lines)) return("")
          tail_output <- paste(utils::tail(lines, 12L), collapse = "\n")
          paste0("Process ", id, " output:\n",
            substr(tail_output, max(1L, nchar(tail_output) - 1999L), nchar(tail_output)))
        }, character(1))
        output <- output[nzchar(output)]
        if (length(output)) message <- paste(message, paste(output, collapse = "\n"), sep = "\n")
      }
    }
    stop(message, call. = FALSE)
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
