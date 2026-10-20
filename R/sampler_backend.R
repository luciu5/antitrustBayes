## Sampler access ends here. Economics code consumes ordinary R draws and
## never needs a live native sampler object after estimation.
.bayes_backend_name <- function(fit) {
  if (inherits(fit, "NutpieFit")) return("nutpie")
  if (inherits(fit, "CmdStanFit")) return("cmdstanr")
  stop("Unsupported sampler fit; expected a CmdStanR or Nutpie fit.",
       call. = FALSE)
}

.bayes_backend_variable_names <- function(fit) {
  if (inherits(fit, "NutpieFit")) {
    ## Retained columns, rather than compiled declarations, are authoritative.
    columns <- dimnames(fit$draws_array)[[3L]]
    ## CmdStan's stan_variables excludes sampler variables such as lp__.
    columns <- columns[!grepl("__$", columns)]
    return(unique(sub("\\[.*$", "", columns)))
  }
  if (inherits(fit, "CmdStanFit")) {
    return(unique(as.character(fit$metadata()$stan_variables)))
  }
  ## Preserve discovery for historical serialized S4 fits and test doubles.
  slots <- methods::slotNames(fit)
  if ("sim" %in% slots && is.list(fit@sim) &&
      is.character(fit@sim$pars_oi) && length(fit@sim$pars_oi)) {
    return(unique(as.character(fit@sim$pars_oi)))
  }
  if ("model_pars" %in% slots && is.character(fit@model_pars)) {
    return(unique(as.character(fit@model_pars)))
  }
  character()
}

.bayes_backend_draws <- function(fit, variables = NULL,
                                 format = "draws_array") {
  format <- match.arg(format, c("draws_array", "draws_matrix", "draws_df",
                                "draws_list", "draws_rvars"))
  if (!inherits(fit, "NutpieFit")) {
    return(fit$draws(variables = variables, format = format))
  }
  draws <- posterior::as_draws_array(fit$draws_array)
  if (!is.null(variables)) {
    variables <- unique(as.character(variables))
    columns <- posterior::variables(draws)
    selected <- unique(unlist(lapply(variables, function(variable) {
      columns[columns == variable | startsWith(columns, paste0(variable, "["))]
    }), use.names = FALSE))
    missing <- variables[!vapply(variables, function(variable) {
      any(columns == variable | startsWith(columns, paste0(variable, "[")))
    }, logical(1))]
    if (length(missing)) {
      stop("Nutpie fit is missing posterior field(s): ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    draws <- posterior::subset_draws(draws, variable = selected)
  }
  switch(format,
    draws_array = draws,
    draws_matrix = posterior::as_draws_matrix(draws),
    draws_df = posterior::as_draws_df(draws),
    draws_list = posterior::as_draws_list(draws),
    draws_rvars = posterior::as_draws_rvars(draws)
  )
}

## Stan names carry their indices explicitly. Never infer matrix ordering
## from the position of a column in an upstream sampler's output.
.bayes_extract_indexed_draws <- function(matrix, variables) {
  columns <- colnames(matrix)
  if (is.null(columns) || anyNA(columns) || anyDuplicated(columns)) {
    stop("Posterior draw columns must have unique Stan names.", call. = FALSE)
  }
  result <- stats::setNames(vector("list", length(variables)), variables)
  for (variable in variables) {
    selected <- which(columns == variable |
                        startsWith(columns, paste0(variable, "[")))
    if (!length(selected)) {
      stop("Sampler fit is missing posterior field: ", variable, ".",
           call. = FALSE)
    }
    labels <- columns[selected]
    if (length(labels) == 1L && identical(labels, variable)) {
      result[[variable]] <- as.numeric(matrix[, selected])
      next
    }
    suffix <- substring(labels, nchar(variable) + 1L)
    if (any(!grepl("^\\[[1-9][0-9]*(,[1-9][0-9]*)*\\]$", suffix))) {
      stop("Invalid indexed posterior field: ", variable, ".", call. = FALSE)
    }
    indexes <- lapply(substring(suffix, 2L, nchar(suffix) - 1L), function(x)
      suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])))
    rank <- unique(lengths(indexes))
    if (length(rank) != 1L || anyNA(unlist(indexes))) {
      stop("Invalid indexed posterior field: ", variable, ".", call. = FALSE)
    }
    dimensions <- vapply(seq_len(rank), function(axis)
      max(vapply(indexes, `[[`, integer(1), axis)), integer(1))
    if (prod(dimensions) != length(selected)) {
      stop("Incomplete indexed posterior field: ", variable, ".", call. = FALSE)
    }
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

.bayes_backend_extract <- function(fit, variables = NULL) {
  ## Keep this compatibility entrypoint for existing CmdStan objects/mocks.
  if (!inherits(fit, "NutpieFit")) return(.cmdstan_extract(fit, variables))
  if (is.null(variables)) variables <- .bayes_backend_variable_names(fit)
  variables <- unique(as.character(variables))
  matrix <- as.matrix(.bayes_backend_draws(fit, variables, "draws_matrix"))
  .bayes_extract_indexed_draws(matrix, variables)
}

.bayes_backend_summary <- function(fit, variables = NULL) {
  if (!inherits(fit, "NutpieFit")) return(fit$summary(variables = variables))
  as.data.frame(posterior::summarise_draws(
    .bayes_backend_draws(fit, variables, "draws_array")
  ))
}

.bayes_backend_diagnostics <- function(fit, format = "draws_array") {
  format <- match.arg(format, c("draws_array", "draws_matrix"))
  if (!inherits(fit, "NutpieFit")) {
    return(fit$sampler_diagnostics(format = format))
  }
  if (is.null(fit$diagnostics_array)) {
    stop("Nutpie fit has no retained sampler diagnostics.", call. = FALSE)
  }
  if (format == "draws_array") {
    posterior::as_draws_array(fit$diagnostics_array)
  } else {
    posterior::as_draws_matrix(fit$diagnostics_array)
  }
}

.bayes_backend_sampling_contract <- function(fit, requested_chains,
                                              inference = "hmc") {
  if (inherits(fit, "NutpieFit")) {
    .nutpie_assert_fit(fit, inference, requested_chains)
  } else {
    .cmdstan_assert_fit(fit, inference, requested_chains)
  }
  invisible(fit)
}

.bayes_backend_provenance <- function(fit) {
  metadata <- if (inherits(fit, "NutpieFit")) fit$metadata else fit$metadata()
  list(backend = .bayes_backend_name(fit), metadata = metadata,
       provenance = metadata$provenance)
}

## Timing scopes differ between native samplers. Preserve that distinction;
## in particular a native warmup-plus-sampling total is not sampling-only time.
.bayes_backend_sampling_info <- function(fit) {
  backend <- .bayes_backend_name(fit)
  metadata <- .bayes_backend_provenance(fit)$metadata
  shape <- dim(.bayes_backend_draws(fit))
  if (backend == "nutpie") {
    timing <- list(compile_seconds = metadata$compile_time %||% NA_real_,
      warmup_seconds = NA_real_, sampling_seconds = NA_real_,
      total_seconds = metadata$sampling_time_native %||% NA_real_,
      scope = metadata$sampling_time_scope,
      unavailable = list(phase_times = "nutpieR does not expose separate warmup and sampling times."))
  } else {
    native <- if (is.function(fit$time)) fit$time() else list()
    timing <- list(compile_seconds = NA_real_,
      warmup_seconds = native$chains$warmup %||% NA_real_,
      sampling_seconds = native$chains$sampling %||% NA_real_,
      total_seconds = native$total %||% NA_real_,
      scope = "CmdStanR process elapsed time; phase times are per chain. Compilation is excluded.",
      unavailable = list(compile_time = "CmdStan fit objects do not retain compilation time."))
  }
  list(backend = backend, chains = shape[[2L]],
    retained_draws_per_chain = shape[[1L]],
    warmup_iterations = metadata$num_warmup %||% metadata$iter_warmup %||% NA_integer_,
    requested_controls = if (backend == "nutpie") metadata$resolved_options else NULL,
    native_config = if (backend == "nutpie") metadata$sampler_config else metadata,
    timing = timing)
}

## Threshold-free diagnostics. Missing fields stay unavailable, and native
## max-depth flags take precedence over a backend-specific depth convention.
.bayes_backend_diagnostic_summary <- function(fit, max_depth = NULL) {
  diagnostics <- .bayes_backend_diagnostics(fit)
  shape <- dim(diagnostics)
  fields <- dimnames(diagnostics)[[3L]]
  metadata <- .bayes_backend_provenance(fit)$metadata
  if (is.null(max_depth)) max_depth <-
    if (inherits(fit, "NutpieFit")) metadata$resolved_options$max_treedepth else
      metadata$max_treedepth
  depth_known <- is.numeric(max_depth) && length(max_depth) == 1L &&
    is.finite(max_depth) && max_depth >= 1
  values <- function(field, chain) {
    if (field %in% fields) as.numeric(diagnostics[, chain, field]) else numeric()
  }
  complete <- function(x) length(x) > 0L && all(is.finite(x))
  chains <- lapply(seq_len(shape[[2L]]), function(chain) {
    divergence <- values("divergent__", chain)
    hit <- values("maxdepth_reached__", chain)
    if (!length(hit) && depth_known) {
      depth <- values("treedepth__", chain)
      hit <- if (length(depth)) depth >= max_depth else numeric()
    }
    energy <- values("energy__", chain)
    stepsize <- values("stepsize__", chain)
    variance <- if (complete(energy) && length(energy) >= 3L) stats::var(energy) else NA_real_
    data.frame(chain = chain, retained_draws = shape[[1L]],
      divergences = if (complete(divergence)) sum(divergence) else NA_real_,
      depth_hits = if (complete(hit)) sum(hit) else NA_real_,
      bfmi = if (is.finite(variance) && variance > 0)
        sum(diff(energy)^2) / length(energy) / variance else NA_real_,
      stepsize_median = if (complete(stepsize)) stats::median(stepsize) else NA_real_)
  })
  chains <- do.call(rbind, chains)
  reasons <- list(
    divergences = "Divergence flags are absent or nonfinite.",
    depth_hits = "Depth-hit flags or a known maximum depth are unavailable, or depth values are nonfinite.",
    bfmi = "Hamiltonian energy is absent/nonfinite, has fewer than three retained draws, or has zero variance.",
    stepsize_median = "Per-transition step sizes are absent or nonfinite.")
  unavailable <- reasons[vapply(names(reasons), function(field)
    anyNA(chains[[field]]), logical(1))]
  list(backend = .bayes_backend_name(fit), chains = chains,
    posterior = .bayes_backend_summary(fit), unavailable = unavailable,
    bfmi_definition = "sum(diff(energy)^2) / length(energy) / var(energy), separately by retained chain; matches CmdStanR.",
    adaptation = metadata$resolved_options$adaptation %||% metadata$metric %||% NA_character_,
    sampling = .bayes_backend_sampling_info(fit))
}
