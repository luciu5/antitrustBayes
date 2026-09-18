.bayes_is_collapsed_quality <- function(fit) {
  metadata <- fit@metadata
  isTRUE(metadata$bayes_is_collapsed_quality) ||
    identical(metadata$quality_parameterization, "collapsed") ||
    isTRUE(fit@stan_data$quality_parameterization == 2L)
}

validate_comparison_contract <- function(fits) {
  reference <- fits[[1L]]
  columns <- c(
    "prices", "quantities", "margins", "market", "year", "group",
    "loo_cluster", "margin_observed", "outside_price"
  )
  missing <- setdiff(columns, names(reference@data))
  if (length(missing)) {
    stop("The reference fit lacks comparison data: ",
         paste(missing, collapse = ", "), ". Refit all models.")
  }

  for (i in seq_along(fits)[-1L]) {
    candidate <- fits[[i]]
    if (!identical(candidate@market_side, reference@market_side)) {
      stop("Models must use the same market_side for comparison.")
    }
    if (!all(columns %in% names(candidate@data)) ||
        !isTRUE(all.equal(
          candidate@data[columns], reference@data[columns],
          check.attributes = FALSE
        ))) {
      stop("Models must be fitted to the same ordered observations, margin ",
           "mask, identifiers, and outside prices.")
    }
    common_outcomes <- c(
      "shareIn", "clr_shareIn", "marginNormInv", "log_marginNormInv",
      "margin_obs", "market_year_idx", "cluster", "likelihood_mode"
    )
    if (!all(common_outcomes %in% names(candidate@stan_data)) ||
        !all(common_outcomes %in% names(reference@stan_data)) ||
        !isTRUE(all.equal(
          candidate@stan_data[common_outcomes],
          reference@stan_data[common_outcomes],
          check.attributes = FALSE
        ))) {
      stop("Models do not share the current quantity-space likelihood contract.")
    }
    shared_prior_names <- setdiff(
      names(reference@metadata$resolved_priors),
      grep("^bargaining_", names(reference@metadata$resolved_priors), value = TRUE)
    )
    if (!(.bayes_uses_cutoff_partition(reference) &&
          .bayes_uses_cutoff_partition(candidate))) {
      shared_prior_names <- setdiff(shared_prior_names,
        c("prior_cutoff_alpha", "prior_cutoff_beta"))
    }
    if (!(.bayes_is_smooth_fringe(reference) &&
          .bayes_is_smooth_fringe(candidate))) {
      shared_prior_names <- setdiff(shared_prior_names, "prior_sigma_meanval_fringe")
    }
    if (!identical(candidate@metadata$prior_manifest,
                   reference@metadata$prior_manifest) ||
        !isTRUE(all.equal(
          candidate@metadata$resolved_priors[shared_prior_names],
          reference@metadata$resolved_priors[shared_prior_names]
        ))) {
      stop("Models must use the same prior manifest on shared parameters.")
    }
  }
  invisible(TRUE)
}

.bayes_bridge_kernel_names <- c(
  "ces_stackelberg", "ces_supply", "logit_stackelberg", "logit_supply"
)

## These are the only source migrations for which an old saved posterior may
## be evaluated by a current package.  Keep this allow-list in the package so
## bridge behavior does not depend on an operator-maintained /tmp file.
.bayes_bridge_approved_hashes <- list(
  original = c(
    ces_stackelberg =
      "6e10cbea685736fabb13f325eaf5295a61ee19cfe33218af326c912b6b9d64d5",
    ces_supply =
      "44b9dd3944ee777c7ea6fc3f988da980dd12cba3bc0e3d030909f0f287b0e266",
    logit_stackelberg =
      "d224c0343a93762225f0fcdc2135dddc059dcd7b09020fe4c6763505ab093a22",
    logit_supply =
      "66103237317f96d4f14ba38ce9426f2c5991bba7fbcd6b437234242fb07ee7ce"
  ),
  rotation = c(
    ces_stackelberg =
      "d14c59ad142aa441ed51c0c9364f2602866f1de398f83746611fa87c73db1800",
    ces_supply =
      "b87f921d5f781245b57e7db90d12f30c315c188e37e6f3ff17e9bffd02065cf9",
    logit_stackelberg =
      "6b1c2ca570bcb62a6b3765979e9af873aac531960382568680306e1654e4c521",
    logit_supply =
      "1c759485e68a549fbe754623a4c615a3c36afed577c8861a8632819c02fc5313"
  ),
  helmert = c(
    ces_stackelberg =
      "5906d0e9f7a00c6ec29da47cb9e870bddb65f8750f53855286b068aa24de73e6",
    ces_supply =
      "e4a759841324154c27bdfc2508cf0babae3a0daaba0c942aefc1289deaae3592",
    logit_stackelberg =
      "4e91beeba1f0155781d72a824bda92f9e1ad51b901d583f9700056b2bb3bcc76",
    logit_supply =
      "59cdf553449bfe192aa8e13398fd54bcaf9b3e991289e2e9fd2afc0e93994ff3"
  )
)

.bayes_bridge_preserve_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

.bayes_bridge_read_hashes <- function(value) {
  ## Kept as a small validator for diagnostics.  Runtime restoration uses the
  ## package-owned allow-list above and never reads an operator-maintained file.
  if (!is.character(value) || is.null(names(value)) ||
      length(value) != length(.bayes_bridge_kernel_names) ||
      anyDuplicated(names(value)) ||
      !setequal(names(value), .bayes_bridge_kernel_names) ||
      any(!grepl("^[0-9a-f]{64}$", unname(value)))) {
    return(NULL)
  }
  value[.bayes_bridge_kernel_names]
}

.bayes_bridge_is_approved_migration <- function(saved_hash, current_hash,
                                                kernel, original, rotation,
                                                helmert = NULL) {
  !is.null(.bayes_bridge_migration_geometry(
    saved_hash, current_hash, kernel, original, rotation, helmert
  ))
}

.bayes_bridge_migration_geometry <- function(saved_hash, current_hash,
                                              kernel, original, rotation,
                                              helmert = NULL) {
  valid <- is.character(saved_hash) && length(saved_hash) == 1L &&
    is.character(current_hash) && length(current_hash) == 1L &&
    is.character(kernel) && length(kernel) == 1L &&
    kernel %in% .bayes_bridge_kernel_names &&
    is.character(original) && is.character(rotation) &&
    kernel %in% names(original) && kernel %in% names(rotation)
  if (!valid || identical(saved_hash, current_hash)) return(NULL)
  old_to_rotation <- identical(unname(original[[kernel]]), saved_hash) &&
    identical(unname(rotation[[kernel]]), current_hash)
  new_helmert <- is.character(helmert) && kernel %in% names(helmert) &&
    identical(unname(helmert[[kernel]]), current_hash)
  if (old_to_rotation || (new_helmert &&
      identical(unname(original[[kernel]]), saved_hash))) return("legacy")
  if (new_helmert && identical(unname(rotation[[kernel]]), saved_hash)) {
    return("preserve")
  }
  NULL
}

.bayes_bridge_code_contract <- function(stanfit) {
  if (!methods::is(stanfit, "stanfit")) {
    stop("Bridge context restoration requires a stanfit object.", call. = FALSE)
  }
  kernel <- as.character(stanfit@model_name)
  if (length(kernel) != 1L || !kernel %in% .bayes_bridge_kernel_names) {
    stop("Saved Stan model name is not a supported antitrustBayes kernel: ",
         kernel, call. = FALSE)
  }
  saved_code <- tryCatch(stanfit@stanmodel@model_code,
                          error = function(e) NULL)
  if (!is.character(saved_code) || length(saved_code) != 1L ||
      !nzchar(saved_code)) {
    stop("Saved stanfit has no recoverable Stan model code.", call. = FALSE)
  }
  models <- tryCatch(get("stanmodels", envir = asNamespace("antitrustBayes")),
                     error = function(e) NULL)
  current <- if (is.null(models)) NULL else models[[kernel]]
  if (is.null(current) || !methods::is(current, "stanmodel")) {
    stop("Installed antitrustBayes package has no matching kernel ", kernel,
         ".", call. = FALSE)
  }
  if (!identical(as.character(current@model_name), kernel)) {
    stop("Installed kernel name does not match the saved Stan model name.",
         call. = FALSE)
  }
  saved_hash <- digest::digest(saved_code, algo = "sha256", serialize = FALSE)
  current_hash <- digest::digest(current@model_code, algo = "sha256",
                                 serialize = FALSE)
  same_code <- identical(saved_hash, current_hash)
  migration <- FALSE
  migration_geometry <- NULL
  if (!same_code) {
    original <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$original)
    rotation <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$rotation)
    helmert <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$helmert)
    migration_geometry <- .bayes_bridge_migration_geometry(
      saved_hash, current_hash, kernel, original, rotation, helmert
    )
    migration <- !is.null(migration_geometry)
    if (!migration) {
      stop("Saved and installed Stan code hashes differ for ", kernel,
           "; no approved saved-to-current kernel pair matches.", call. = FALSE)
    }
  }
  list(kernel = kernel, current_model = current, saved_hash = saved_hash,
       current_hash = current_hash, same_code = same_code,
       migration = migration, migration_geometry = migration_geometry)
}

.bayes_bridge_canonicalize_data <- function(data) {
  ## Old serialized fits can contain an optional field twice (a NULL entry
  ## followed by its actual value).  Remove only absent NULL entries; every
  ## remaining duplicate must agree before the private context is built.
  nms <- names(data)
  if (length(nms) && anyDuplicated(nms)) {
    duplicate_names <- unique(nms[duplicated(nms)])
    agree <- function(values, name) {
      values <- values[!vapply(values, is.null, logical(1))]
      if (!length(values)) return(NULL)
      reference <- values[[1L]]
      same <- vapply(values[-1L], function(value) {
        if (is.numeric(reference) && is.numeric(value) &&
            length(reference) == length(value)) {
          if (!length(reference)) return(TRUE)
          delta <- abs(as.numeric(value) - as.numeric(reference))
          tol <- 1e-12 + 128 * .Machine$double.eps *
            pmax(1, abs(as.numeric(reference)), abs(as.numeric(value)))
          all(is.finite(delta)) && all(delta <= tol)
        } else {
          identical(value, reference)
        }
      }, logical(1))
      if (length(same) && !all(same)) {
        stop("Private bridge data has conflicting duplicate field '", name,
             "'.", call. = FALSE)
      }
      reference
    }
    for (name in duplicate_names) {
      value <- agree(data[which(nms == name)], name)
      data <- data[nms != name]
      if (!is.null(value)) data[[name]] <- value
      nms <- names(data)
    }
  }
  data
}

.bayes_bridge_private_data <- function(fit, contract) {
  data <- .bayes_bridge_canonicalize_data(fit@stan_data)
  if (isTRUE(contract$migration)) {
    if (identical(contract$migration_geometry, "legacy")) {
      flag <- data$market_s0_parameterization
      if ("market_s0_parameterization" %in% names(data) &&
          (!is.numeric(flag) || length(flag) != 1L || !is.finite(flag) ||
           flag != 1)) {
        stop("Approved original-kernel migration requires a legacy private data flag of 1.",
             call. = FALSE)
      }
      ## Only a pre-rotation source has this implicit legacy coordinate system.
      data$market_s0_parameterization <- 1L
    } else if (identical(contract$migration_geometry, "preserve")) {
      flag <- data$market_s0_parameterization
      if (!is.numeric(flag) || length(flag) != 1L || !is.finite(flag) ||
          !flag %in% c(1, 2)) {
        stop("Approved performance migration requires its saved private data flag of 1 or 2.",
             call. = FALSE)
      }
    } else {
      stop("Approved bridge migration has no validated coordinate policy.",
           call. = FALSE)
    }
  }
  data
}

.bayes_bridge_is_zero_width <- function(dims) {
  is.numeric(dims) && length(dims) == 1L &&
    isTRUE(as.numeric(dims) == 0)
}

.bayes_bridge_dim_width <- function(dims) {
  if (!is.numeric(dims) || any(!is.finite(dims)) || any(dims < 0) ||
      any(dims != as.integer(dims))) {
    stop("Stan parameter dimensions must be finite non-negative integers.",
         call. = FALSE)
  }
  if (!length(dims)) return(1)
  prod(as.numeric(dims))
}

.bayes_bridge_raw_parameter_names <- function(stanfit) {
  ## model_pars is the canonical order saved by rstan, but it also contains
  ## transformed parameters, generated quantities, and lp__.  The latter are
  ## not initialization fields and must not participate in the bridge schema
  ## contract.  When the model code is available, use its parameters block to
  ## select the raw subset while retaining model_pars' order.  Small test
  ## doubles may not carry a stanmodel slot; their model_pars is already the
  ## raw subset, so retain that explicit schema as a conservative fallback.
  model_pars <- stanfit@model_pars
  if (!is.character(model_pars) || !length(model_pars) ||
      anyDuplicated(model_pars)) {
    stop("Stan fit has no unique raw parameter schema.", call. = FALSE)
  }

  code <- NULL
  if ("stanmodel" %in% methods::slotNames(stanfit)) {
    code <- tryCatch(stanfit@stanmodel@model_code,
                     error = function(e) NULL)
  }
  if (is.null(code) && "stan_code" %in% methods::slotNames(stanfit)) {
    code <- tryCatch(stanfit@stan_code, error = function(e) NULL)
  }
  declared <- character(0L)
  if (is.character(code) && length(code) == 1L && nzchar(code)) {
    lines <- strsplit(code, "\\n", fixed = FALSE)[[1L]]
    parameter_start <- grep("^\\s*parameters\\s*\\{", lines,
                            perl = TRUE)
    transformed_start <- grep("^\\s*transformed[[:space:]]+parameters\\s*\\{",
                              lines, perl = TRUE)
    parameter_start <- parameter_start[1L]
    transformed_start <- transformed_start[transformed_start > parameter_start][1L]
    if (is.finite(parameter_start) && is.finite(transformed_start) &&
        transformed_start > parameter_start + 1L) {
      block <- paste(lines[(parameter_start + 1L):(transformed_start - 1L)],
                     collapse = "\n")
      block <- gsub("(?m)//.*$", "", block, perl = TRUE)
      statements <- strsplit(block, ";", fixed = TRUE)[[1L]]
      declared <- vapply(statements, function(statement) {
        statement <- gsub("\\s+", " ", trimws(statement), perl = TRUE)
        if (!nzchar(statement)) return(NA_character_)
        token <- regmatches(statement,
                            regexpr("[A-Za-z][A-Za-z0-9_]*$", statement,
                                    perl = TRUE))
        if (!length(token)) NA_character_ else token
      }, character(1L))
      declared <- unique(declared[!is.na(declared)])
    }
  }

  if (length(declared)) {
    raw_names <- model_pars[model_pars %in% declared]
    if (!length(raw_names)) {
      stop("Stan model code has no raw parameters in model_pars.",
           call. = FALSE)
    }
    return(raw_names)
  }

  ## Explicit mock objects used by package tests have no model code.  Keep
  ## their declared model_pars contract, while excluding lp__ if supplied.
  raw_names <- model_pars[model_pars != "lp__"]
  if (!length(raw_names)) {
    stop("Stan fit has no raw parameter declarations.", call. = FALSE)
  }
  raw_names
}

.bayes_bridge_positive_init_schema <- function(stanfit) {
  if (!length(stanfit@inits) || !is.list(stanfit@inits[[1L]]) ||
      is.null(names(stanfit@inits[[1L]]))) {
    stop("Stan fit has no named initialization schema.", call. = FALSE)
  }
  raw_names <- .bayes_bridge_raw_parameter_names(stanfit)
  init_names <- names(stanfit@inits[[1L]])
  dims <- stanfit@par_dims
  if (is.null(names(dims)) || any(!raw_names %in% names(dims))) {
    stop("Stan raw parameter schema contains an undeclared parameter.",
         call. = FALSE)
  }
  widths <- vapply(raw_names, function(name) {
    .bayes_bridge_dim_width(dims[[name]])
  }, numeric(1))
  positive <- raw_names[widths > 0]
  if (any(!positive %in% init_names)) {
    stop("Stan initialization is missing a declared raw parameter.",
         call. = FALSE)
  }
  list(
    names = positive,
    dims = lapply(dims[positive], function(value) as.integer(value))
  )
}

.bayes_bridge_assert_positive_schema <- function(saved_stanfit, context) {
  saved <- .bayes_bridge_positive_init_schema(saved_stanfit)
  current <- .bayes_bridge_positive_init_schema(context)
  if (!identical(saved$names, current$names) ||
      !identical(saved$dims, current$dims)) {
    stop("Saved and restored Stan contexts have different nonzero parameter ",
         "names or dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

.bayes_bridge_clone_sample_view <- function(saved_stanfit, context) {
  if (!methods::is(saved_stanfit, "stanfit") ||
      !methods::is(context, "stanfit")) {
    stop("Bridge sample schema adaptation requires stanfit objects.",
         call. = FALSE)
  }
  .bayes_bridge_assert_positive_schema(saved_stanfit, context)
  candidates <- intersect(c("c_s0", "c_s0_raw"), context@sim$pars_oi)
  if (!length(candidates)) {
    return(list(stanfit = saved_stanfit, added = character(0L)))
  }
  missing <- setdiff(candidates, saved_stanfit@sim$pars_oi)
  if (!length(missing)) {
    return(list(stanfit = saved_stanfit, added = character(0L)))
  }
  missing_dims <- context@par_dims[missing]
  if (any(!vapply(missing_dims, .bayes_bridge_is_zero_width, logical(1)))) {
    stop("Current bridge context adds a nonzero c_s0 field absent from the ",
         "saved posterior.", call. = FALSE)
  }
  if (any(missing %in% names(saved_stanfit@par_dims))) {
    stop("Saved bridge schema contains an unexpected c_s0 field.",
         call. = FALSE)
  }

  ## S4 assignment creates a private copy under R's copy-on-write rules.  The
  ## empty fields add no draw columns or flat names, but let bridgesampling's
  ## skeleton include the current model's declared zero-width parameters.
  private <- saved_stanfit
  private@par_dims <- c(
    saved_stanfit@par_dims,
    setNames(lapply(missing, function(name) context@par_dims[[name]]), missing)
  )
  private@sim$pars_oi <- c(saved_stanfit@sim$pars_oi, missing)
  private@sim$dims_oi <- c(
    saved_stanfit@sim$dims_oi,
    setNames(lapply(missing, function(name) {
      value <- context@sim$dims_oi[[name]]
      if (is.null(value)) value <- context@par_dims[[name]]
      value
    }), missing)
  )
  for (chain in seq_along(private@sim$samples)) {
    for (name in missing) private@sim$samples[[chain]][[name]] <- numeric(0L)
  }

  before <- rstan::extract(saved_stanfit, permuted = FALSE, inc_warmup = FALSE)
  after <- rstan::extract(private, permuted = FALSE, inc_warmup = FALSE)
  if (!identical(before, after)) {
    stop("Private zero-width bridge adapter changed posterior sample values.",
         call. = FALSE)
  }
  list(stanfit = private, added = missing)
}

.bayes_bridge_probe_live <- function(stanfit) {
  tryCatch({
    n_upars <- as.integer(rstan::get_num_upars(stanfit))
    if (length(n_upars) != 1L || is.na(n_upars) || n_upars < 0L) {
      stop("invalid unconstrained parameter count")
    }
    if (!length(stanfit@inits)) stop("saved stanfit has no initialization")
    upars <- rstan::unconstrain_pars(stanfit, stanfit@inits[[1L]])
    if (!is.numeric(upars) || length(upars) != n_upars) {
      stop("saved initialization has an incompatible unconstrained schema")
    }
    stanfit
  }, error = function(e) NULL)
}

.bayes_bridge_fresh_context <- function(fit, contract) {
  rstan::sampling(
    object = contract$current_model,
    data = .bayes_bridge_private_data(fit, contract),
    chains = 1L,
    iter = 1L,
    warmup = 0L,
    algorithm = "Fixed_param",
    init = 0,
    seed = 1L,
    refresh = 0L,
    cores = 1L
  )
}

.bayes_bridge_theta_at <- function(stanfit, chain, iteration) {
  pars <- stanfit@sim$pars_oi
  flat_names <- unlist(lapply(pars, function(par) {
    stanfit@sim$fnames_oi[
      stanfit@sim$fnames_oi == par |
        startsWith(stanfit@sim$fnames_oi, paste0(par, "["))
    ]
  }), use.names = FALSE)
  samples <- stanfit@sim$samples[[chain]]
  if (!length(flat_names) || !all(flat_names %in% names(samples))) {
    stop("Saved stanfit samples do not match their declared parameter schema.",
         call. = FALSE)
  }
  values <- vapply(flat_names, function(name) {
    as.numeric(samples[[name]][[iteration]])
  }, numeric(1))
  skeleton <- bridgesampling:::.create_skeleton(
    pars, stanfit@par_dims[pars]
  )
  bridgesampling:::.rstan_relist(values, skeleton)
}

.bayes_bridge_saved_rows <- function(stanfit, n_rows = 5L) {
  chains <- stanfit@sim$samples
  if (!length(chains)) stop("Saved stanfit has no posterior chains.", call. = FALSE)
  rows <- lapply(seq_along(chains), function(chain) {
    lp <- chains[[chain]][["lp__"]]
    if (is.null(lp) || !length(lp)) {
      stop("Saved stanfit chain ", chain, " has no lp__ draws.", call. = FALSE)
    }
    n <- length(lp)
    warmup <- if (length(stanfit@sim$warmup2) >= chain) {
      as.integer(stanfit@sim$warmup2[[chain]])
    } else 0L
    ## Depending on the rstan serialization path, samples can contain either
    ## post-warmup draws or all iterations.  Trim an all-iteration chain using
    ## its per-chain warmup count; when samples are already post-warmup (the
    ## usual path, n <= warmup2), retain all saved posterior rows.
    post_start <- if (is.finite(warmup) && warmup >= 0L && n > warmup) {
      warmup + 1L
    } else 1L
    post_idx <- seq.int(post_start, n)
    idx <- unique(as.integer(round(seq(post_idx[1L], post_idx[length(post_idx)],
      length.out = min(as.integer(n_rows), length(post_idx))))))
    data.frame(chain = chain, iteration = idx, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.bayes_bridge_validate_target <- function(saved_stanfit, context) {
  n_upars <- as.integer(rstan::get_num_upars(context))
  if (length(n_upars) != 1L || is.na(n_upars) || n_upars < 0L) {
    stop("Restored context has an invalid unconstrained parameter count.",
         call. = FALSE)
  }
  rows <- .bayes_bridge_saved_rows(saved_stanfit, n_rows = 5L)
  checks <- lapply(seq_len(nrow(rows)), function(i) {
    chain <- rows$chain[[i]]
    iteration <- rows$iteration[[i]]
    theta <- tryCatch(
      .bayes_bridge_theta_at(saved_stanfit, chain, iteration),
      error = function(e) stop("Cannot reconstruct saved draw at chain ", chain,
                               ", iteration ", iteration, ": ",
                               conditionMessage(e), call. = FALSE)
    )
    upars <- tryCatch(
      rstan::unconstrain_pars(context, theta),
      error = function(e) stop("Restored context rejects saved schema at chain ",
                               chain, ", iteration ", iteration, ": ",
                               conditionMessage(e), call. = FALSE)
    )
    if (!is.numeric(upars) || length(upars) != n_upars) {
      stop("Restored context has a nonzero parameter schema mismatch at chain ",
           chain, ", iteration ", iteration, call. = FALSE)
    }
    restored <- tryCatch(
      rstan::log_prob(context, upars, adjust_transform = TRUE, gradient = FALSE),
      error = function(e) stop("Restored log target failed at chain ", chain,
                               ", iteration ", iteration, ": ",
                               conditionMessage(e), call. = FALSE)
    )
    saved <- as.numeric(saved_stanfit@sim$samples[[chain]][["lp__"]][[iteration]])
    tolerance <- 1e-7 + 128 * .Machine$double.eps * max(1, abs(saved))
    if (!is.finite(saved) || !is.finite(restored) ||
        abs(restored - saved) > tolerance) {
      stop("Restored normalized log target differs from saved lp__ at chain ",
           chain, ", iteration ", iteration, " (difference = ",
           format(abs(restored - saved), digits = 16), ", tolerance = ",
           format(tolerance, digits = 16), ").", call. = FALSE)
    }
    data.frame(chain = chain, iteration = iteration, saved_lp = saved,
               restored_lp = restored, abs_difference = abs(restored - saved),
               tolerance = tolerance, stringsAsFactors = FALSE)
  })
  do.call(rbind, checks)
}

.bayes_bridge_restore_context <- function(fit) {
  .bayes_bridge_preserve_rng({
    st <- fit@stanfit
    contract <- .bayes_bridge_code_contract(st)
    context <- if (!isTRUE(contract$migration)) .bayes_bridge_probe_live(st) else NULL
    context_source <- if (!is.null(context)) "saved_live_stanfit" else
      "fresh_fixed_param_context"
    if (is.null(context)) context <- .bayes_bridge_fresh_context(fit, contract)
    sample_view <- .bayes_bridge_clone_sample_view(st, context)
    target_check <- .bayes_bridge_validate_target(sample_view$stanfit, context)
    list(stanfit_model = context, samples = sample_view$stanfit,
         zero_width_fields_added = sample_view$added,
         kernel = contract$kernel,
         saved_hash = contract$saved_hash, current_hash = contract$current_hash,
         migration = contract$migration, source = context_source,
         target_check = target_check)
  })
}

#' Compare Bayesian Models or Exact Ordered-k Restrictions
#'
#' Computes LOO-IC, WAIC, and optionally log marginal likelihoods via bridge
#' sampling for one or more fitted \code{\linkS4class{BayesLogit}} objects.
#' For one exact ordered-k fit, \code{restrictions} instead compares restrictions
#' on k with the exact encompassing-prior identity. Returns a comparison table
#' appropriate to the selected route.
#'
#' @param ... One or more \code{BayesLogit} objects, optionally named. With
#'   \code{universe}, a single named list of source fits is also accepted.
#' @param restrictions Optional named restrictions on ordered-k states. When
#'   supplied, \code{...} must contain exactly one exact ordered-k
#'   \code{\linkS4class{BayesFit}}. See
#'   \code{\link{bayes_restriction_compare}} for the one- and multimarket forms.
#' @param reference Optional restriction name or row index used as the Bayes
#'   factor reference. Only available when \code{restrictions} is supplied.
#' @param restriction_draws Optional generated-quantity draws passed to
#'   \code{\link{bayes_restriction_compare}}. This is mainly useful for testing
#'   or recovering summaries from serialized draws.
#' @param bridge Logical. If \code{TRUE} (default), compute log marginal
#'   likelihoods via \code{\link[bridgesampling]{bridge_sampler}} for Bayes
#'   factor model comparison. If \code{FALSE}, only LOO/WAIC are computed.
#' @param bridge_maxiter Maximum iterations for bridge sampling.
#'   Default 5000.
#' @param bridge_repetitions Number of independently initialized bridge runs
#'   per model. The reported log marginal likelihood is their median.
#' @param loo_unit Natural predictive unit: observation, period, market, or
#'   cluster. Cluster uses the identifiers supplied to \code{bayes_est}.
#' @param prior_weights Optional numeric vector of prior model probabilities.
#'   Must be the same length as the number of models and sum to 1. If
#'   \code{NULL} (default), recorded \code{model_prior} values are used when
#'   all are available; otherwise equal probabilities are used.
#' @param universe Optional \code{bayes_model_universe} descriptor. When
#'   supplied, the named fits in \code{...} are its source fits and the
#'   descriptor's explicit scientific prior is used. This route requires
#'   bridge-sampling evidence.
#'
#' @return With multiple fitted models, a data frame containing model name,
#'   supply model, LOOIC, WAIC, stacking weights, and (if
#'   \code{bridge = TRUE}) log marginal likelihoods, Bayes factors, and
#'   posterior model probabilities. With \code{restrictions}, a
#'   \code{bayes_restriction_comparison} table containing prior and posterior
#'   restriction masses and encompassing-prior Bayes factors.
#'
#' @details
#' Four model comparison approaches are provided:
#' \describe{
#'   \item{LOO-IC}{Leave-one-out cross-validation via \pkg{loo}. Stacking
#'     weights reflect predictive performance.}
#'   \item{WAIC}{Widely Applicable Information Criterion.}
#'   \item{Bridge Sampling}{Computes the log marginal likelihood
#'     \eqn{p(y | M_k)} via \code{\link[bridgesampling]{bridge_sampler}}.
#'     Posterior model probabilities use \code{prior_weights}, or the
#'     model-specific probabilities recorded by \code{bayes_est}.
#'     Bayes factors are computed relative to the best model.}
#'   \item{Exact ordered-k restrictions}{Computes
#'     \eqn{BF_{A,F}=P(A\mid y)/P(A)} from the generated conditional
#'     probabilities in one marginalized-k fit. This route does not use bridge
#'     sampling or a Gelfand-Dey estimator.}
#' }
#'
#' Requires that the Stan model includes a \code{log_lik} generated quantity
#' (for LOO/WAIC). All fits must use the same ordered observations, market
#' side, outside prices, margin mask, and current quantity-space likelihood.
#' This permits Logit/CES comparison without comparing likelihoods defined on
#' different observed outcomes.
#' The current Stan specification retains prior density constants and support
#' normalizers directly. Earlier specification fits must be refitted. Evidence
#' refers to the transformed-outcome working likelihood. Bridge comparison
#' with an active, non-normalized HMT penalty is unsupported. Exact ordered-k
#' LOO clusters must contain whole partition units; observation LOO and period
#' LOO with a market-shared count are unsupported for fringe and Stackelberg.
#'
#' @examples
#' \donttest{
#' # fit1 <- bayes_est(..., supply_model = "bertrand")
#' # fit2 <- bayes_est(..., supply_model = "cournot")
#' # bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # bayes_compare(Bertrand = fit1, Cournot = fit2, bridge = FALSE)
#' # bayes_compare(fit_stack, restrictions = list(k1 = 1L, k2 = 2L,
#' #                                                k3_to_5 = 3:5),
#' #               reference = "k2")
#' }
#'
#' @export
bayes_compare <- function(..., restrictions = NULL, reference = NULL,
                          restriction_draws = NULL,
                          bridge = TRUE, bridge_maxiter = 5000L,
                          bridge_repetitions = 3L,
                          loo_unit = c("cluster", "period", "market", "observation"),
                          prior_weights = NULL, universe = NULL) {
  args <- list(...)
  if (!is.null(universe)) {
    # Accept a named collection as well as named fits in ... for this route.
    if (length(args) == 1L && is.list(args[[1L]]) &&
        !methods::is(args[[1L]], "BayesFit")) args <- args[[1L]]
    if (!is.null(restrictions) || !is.null(reference) ||
        !is.null(restriction_draws)) {
      stop("'universe' cannot be combined with ordered-k 'restrictions', ",
           "'reference', or 'restriction_draws'.", call. = FALSE)
    }
    if (!is.null(prior_weights)) {
      stop("'prior_weights' belongs in the explicit model-universe descriptor ",
           "when 'universe' is supplied.", call. = FALSE)
    }
    loo_unit <- match.arg(loo_unit)
    return(.bayes_compare_universe(
      args, universe = universe, bridge = bridge,
      bridge_maxiter = bridge_maxiter,
      bridge_repetitions = bridge_repetitions, loo_unit = loo_unit
    ))
  }
  if (!is.null(restrictions)) {
    is_fit <- vapply(args, methods::is, logical(1), "BayesFit")
    if (length(args) != 1L || sum(is_fit) != 1L) {
      stop("Restriction comparison requires exactly one BayesFit object in '...'.",
           call. = FALSE)
    }
    if (!is.null(prior_weights)) {
      stop("'prior_weights' compares separate fitted models and is unavailable ",
           "for ordered-k restrictions; the fit's recorded k prior is used.",
           call. = FALSE)
    }
    return(bayes_restriction_compare(
      args[[which(is_fit)]], restrictions = restrictions,
      reference = reference, draws = restriction_draws
    ))
  }
  if (!is.null(reference) || !is.null(restriction_draws)) {
    stop("'reference' and 'restriction_draws' require 'restrictions'.",
         call. = FALSE)
  }
  loo_unit <- match.arg(loo_unit)
  if (length(bridge_repetitions) != 1L || !is.finite(bridge_repetitions) ||
      bridge_repetitions < 1 || bridge_repetitions != as.integer(bridge_repetitions)) {
    stop("'bridge_repetitions' must be a positive integer.")
  }
  bridge_repetitions <- as.integer(bridge_repetitions)

  # Separate BayesLogit objects from named arguments
  fits <- args[vapply(args, is, logical(1), "BayesLogit")]
  if (length(fits) == 0L) {
    stop("At least one BayesLogit fit is required.")
  }
  invisible(lapply(fits, validate_fit_spec))
  validate_comparison_contract(fits)
  collapsed_quality <- any(vapply(fits, .bayes_is_collapsed_quality, logical(1)))
  if (collapsed_quality && !isTRUE(bridge)) {
    stop("LOO/WAIC comparison is unavailable for collapsed quality fits; ",
         "use bridge = TRUE for normalized evidence.", call. = FALSE)
  }
  if (!collapsed_quality) .bayes_comparison_units(fits, loo_unit)
  if (bridge) .bayes_validate_bridge(fits)

  # Handle named and unnamed arguments
  nms <- names(fits)
  if (is.null(nms)) {
    nms <- paste0("Model_", seq_along(fits))
  } else {
    empty <- nms == ""
    nms[empty] <- paste0("Model_", which(empty))
  }

  # Prior model probabilities
  if (is.null(prior_weights)) {
    recorded <- vapply(fits, function(x) {
      value <- x@metadata$model_prior
      if (length(value) == 1L && is.numeric(value)) as.numeric(value) else NA_real_
    }, numeric(1))
    if (all(is.finite(recorded)) && sum(recorded) > 0) {
      prior_weights <- recorded / sum(recorded)
    } else {
      prior_weights <- rep(1 / length(fits), length(fits))
    }
  } else {
    if (length(prior_weights) != length(fits)) {
      stop("'prior_weights' must have the same length as ",
           "the number of models.")
    }
    if (!is.numeric(prior_weights) || any(!is.finite(prior_weights)) ||
        any(prior_weights < 0) || sum(prior_weights) <= 0) {
      stop("'prior_weights' must be finite, non-negative probabilities.")
    }
    if (abs(sum(prior_weights) - 1) > 1e-6) {
      stop("'prior_weights' must sum to 1.")
    }
  }

  supply_names <- c(
    "Bertrand", "Monopolistic", "Cournot", "Auction", "Bargaining"
  )

  loo_list  <- vector("list", length(fits))
  waic_list <- vector("list", length(fits))
  bridge_list <- list()
  bridge_context_info <- vector("list", length(fits))
  parameter_name <- switch(loo_unit,
    observation = "log_lik",
    period = "log_lik_period",
    market = "log_lik_market",
    cluster = "log_lik_cluster"
  )
  loo_available <- !collapsed_quality

  for (i in seq_along(fits)) {
    fit <- fits[[i]]

    if (loo_available) {
      log_lik <- loo::extract_log_lik(fit@stanfit, parameter_name = parameter_name,
                                       merge_chains = FALSE)
      n_units <- dim(log_lik)[3L]
      if (is.null(n_units)) n_units <- 1L
      if (n_units < 2L) {
        loo_available <- FALSE
      } else {
        r_eff <- loo::relative_eff(exp(log_lik))
        loo_list[[i]]  <- loo::loo(log_lik, r_eff = r_eff)
        waic_list[[i]] <- loo::waic(log_lik)
      }
    }

    # Bridge sampling for marginal likelihood
    if (bridge) {
      bridge_context <- tryCatch(
        .bayes_bridge_restore_context(fit),
        error = function(e) {
          warning("Bridge context restoration failed for model ", nms[i],
                  "; evidence marked incomplete: ", conditionMessage(e),
                  call. = FALSE)
          NULL
        }
      )
      if (is.null(bridge_context)) {
        bridge_list[[i]] <- rep(list(NULL), bridge_repetitions)
        bridge_context_info[[i]] <- list(status = "incomplete")
      } else {
        bridge_list[[i]] <- lapply(seq_len(bridge_repetitions), function(rep_id) {
          tryCatch(
            bridgesampling::bridge_sampler(
              bridge_context$samples, stanfit_model = bridge_context$stanfit_model,
              silent = TRUE, maxiter = bridge_maxiter
            ),
            error = function(e) {
              warning("Bridge sampling failed for model ", nms[i],
                      ", repetition ", rep_id, ": ", conditionMessage(e))
              NULL
            }
          )
        })
        bridge_context_info[[i]] <- list(
          status = "validated", source = bridge_context$source,
          kernel = bridge_context$kernel,
          saved_hash = bridge_context$saved_hash,
          current_hash = bridge_context$current_hash,
          migration = bridge_context$migration,
          zero_width_fields_added = bridge_context$zero_width_fields_added,
          target_check = bridge_context$target_check
        )
      }
    }
  }

  names(loo_list) <- names(waic_list) <- nms
  if (bridge) names(bridge_list) <- nms
  if (loo_available) {
    comp <- loo::loo_compare(loo_list)
    weights_stacking <- loo::loo_model_weights(loo_list, method = "stacking")
  } else if (!collapsed_quality) {
    warning("LOO/WAIC require at least two ", loo_unit,
            " units; returning NA predictive summaries.", call. = FALSE)
    comp <- NULL
    weights_stacking <- rep(NA_real_, length(fits))
  } else {
    comp <- NULL
    weights_stacking <- rep(NA_real_, length(fits))
  }

  result <- data.frame(
    model         = nms,
    supply_model  = vapply(fits, function(f)
      supply_names[f@supply_model], character(1)),
    demand_model  = vapply(fits, slot, character(1), "demand_model"),
    market_side   = vapply(fits, slot, character(1), "market_side"),
    share_likelihood = vapply(fits, function(x)
      x@metadata$share_likelihood %||% "student_t", character(1)),
    margin_likelihood = vapply(fits, function(x)
      x@metadata$margin_likelihood %||% "student_t", character(1)),
    quality_parameterization = vapply(fits, function(x)
      x@metadata$quality_parameterization %||% "sampled", character(1)),
    quality_scope = vapply(fits, .quality_scope_from_fit, character(1), strict = FALSE),
    prior_manifest = vapply(fits, function(x) x@metadata$prior_manifest,
                            character(1)),
    model_prior = prior_weights,
    loo_unit = loo_unit,
    looic = if (loo_available) vapply(
      loo_list, function(x) x$estimates["looic", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    looic_se = if (loo_available) vapply(
      loo_list, function(x) x$estimates["looic", "SE"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    p_loo = if (loo_available) vapply(
      loo_list, function(x) x$estimates["p_loo", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    waic = if (loo_available) vapply(
      waic_list, function(x) x$estimates["waic", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    weight_stacking = as.numeric(weights_stacking),
    stringsAsFactors = FALSE
  )

  # Add bridge sampling results
  if (bridge && length(bridge_list) > 0L) {
    bridge_logml <- lapply(bridge_list, function(reps) {
      vapply(reps, function(b) if (!is.null(b)) b$logml else NA_real_, numeric(1))
    })
    logml <- vapply(bridge_logml, function(x) {
      if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
    }, numeric(1))
    logml_sd <- vapply(bridge_logml, function(x) {
      if (sum(is.finite(x)) > 1L) stats::sd(x[is.finite(x)]) else NA_real_
    }, numeric(1))

    finite <- is.finite(logml)
    diff_logml <- bf_vs_best <- post_prob <- rep(NA_real_, length(logml))
    if (any(finite)) {
      max_logml <- max(logml[finite])
      diff_logml[finite] <- logml[finite] - max_logml
      bf_vs_best[finite] <- exp(diff_logml[finite])
      unnorm <- exp(diff_logml[finite]) * prior_weights[finite]
      if (sum(unnorm) > 0) {
        post_prob[finite] <- unnorm / sum(unnorm)
      }
    }

    result$logml        <- logml
    result$logml_sd     <- logml_sd
    result$bf_vs_best   <- bf_vs_best
    result$post_prob    <- post_prob
  }

  rownames(result) <- NULL

  attr(result, "loo_compare") <- comp
  attr(result, "loo_list") <- loo_list
  if (bridge) attr(result, "bridge_list") <- bridge_list
  if (bridge) attr(result, "bridge_context_info") <- bridge_context_info
  attr(result, "prior_weights") <- prior_weights
  attr(result, "loo_unit") <- loo_unit
  result
}

#' Decompose Joint and Conditional-Margin Model Evidence
#'
#' Since the conditional log quantity-share likelihood and shared-parameter
#' priors are common across conduct models, subtracting a separately fitted
#' share-only log marginal likelihood isolates
#' \eqn{\log p(M\mid S,M_k)}.
#'
#' @param comparison A bridge-sampling table returned by \code{bayes_compare}.
#' @param share_logml Scalar or model-aligned vector of log marginal
#'   likelihoods from share-only fits. A scalar is recycled when the share
#'   evidence is common across conduct models.
#' @return The comparison table with common-share and conditional-margin
#'   evidence columns.
#' @export
bayes_evidence_decomposition <- function(comparison, share_logml) {
  if (!"logml" %in% names(comparison)) {
    stop("'comparison' must contain bridge-sampling log marginal likelihoods.")
  }
  if (!is.numeric(share_logml) || !length(share_logml) ||
      any(!is.finite(share_logml)) ||
      !(length(share_logml) %in% c(1L, nrow(comparison)))) {
    stop("'share_logml' must be one finite value or one value per model.")
  }
  if (length(share_logml) == 1L) share_logml <- rep(share_logml, nrow(comparison))
  out <- comparison
  out$logml_share <- share_logml
  out$logml_margin_given_share <- out$logml - share_logml
  out
}

#' Extract Model Weights from a Comparison Table
#'
#' Convenience function to extract weights from a \code{\link{bayes_compare}}
#' result using either bridge sampling posterior probabilities or LOO stacking
#' weights.
#'
#' @param comp A data frame returned by \code{\link{bayes_compare}}.
#' @param method Character: \code{"bridge"} (default) for posterior model
#'   probabilities from bridge sampling, or \code{"stacking"} for LOO
#'   stacking weights. Falls back to stacking if bridge results are
#'   unavailable.
#'
#' @return A named numeric vector of weights summing to 1.
#'
#' @examples
#' \donttest{
#' # comp <- bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # bayes_weights(comp, method = "bridge")
#' # bayes_weights(comp, method = "stacking")
#' }
#'
#' @export
bayes_weights <- function(comp, method = c("bridge", "stacking")) {
  method <- match.arg(method)

  bridge_ok <- "post_prob" %in% names(comp) &&
    any(is.finite(comp$post_prob) & comp$post_prob > 0)
  if (method == "bridge" && bridge_ok) {
    w <- comp$post_prob
    w[!is.finite(w)] <- 0
  } else {
    if (method == "bridge" && !bridge_ok) {
      message("Bridge sampling results not available; using stacking weights.")
    }
    w <- comp$weight_stacking
  }

  names(w) <- comp$model
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop("The requested model weights are not finite non-negative values.")
  }
  w / sum(w)
}


#' Bayesian Model Averaging of Merger Simulations
#'
#' Combines merger simulation results from multiple models using
#' model weights (e.g., from \code{\link{bayes_compare}}). Returns a
#' \code{\linkS4class{BayesSim}} object containing a weighted pool of
#' simulation draws, so that all \code{BayesSim} methods
#' (\code{summary}, \code{elast}, \code{diversion}, \code{CV},
#' \code{HypoMonTest}) work on the model-averaged result.
#'
#' @param sims A list of \code{\linkS4class{BayesSim}} objects.
#' @param weights Numeric vector of model weights (must sum to 1).
#'   Typically \code{post_prob} (from bridge sampling) or
#'   \code{weight_stacking} (from LOO) in the output of
#'   \code{\link{bayes_compare}}.
#' @param n_draws Total number of draws in the pooled result.
#'   Default 1000.
#'
#' @return A \code{\linkS4class{BayesSim}} object with draws resampled
#'   across models proportional to their weights. All \code{BayesSim}
#'   methods (e.g., \code{summary}, \code{elast}, \code{CV}) can be
#'   called on the result.
#'
#' @examples
#' \donttest{
#' # comp <- bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # sim1 <- bayes_sim(fit1, ownerPost = ...)
#' # sim2 <- bayes_sim(fit2, ownerPost = ...)
#' #
#' # # BMA using bridge sampling posterior probabilities:
#' # avg_bridge <- bayes_avg(list(sim1, sim2),
#' #                         weights = bayes_weights(comp, "bridge"))
#' #
#' # # BMA using LOO stacking weights:
#' # avg_loo <- bayes_avg(list(sim1, sim2),
#' #                      weights = bayes_weights(comp, "stacking"))
#' #
#' # summary(avg_bridge)
#' }
#'
#' @export
bayes_avg <- function(sims, weights, n_draws = 1000L) {
  if (!is.list(sims) || length(sims) < 2L) {
    stop("'sims' must be a list of at least 2 BayesSim objects.")
  }
  if (!all(vapply(sims, is, logical(1), "BayesSim"))) {
    stop("Every element of 'sims' must be a BayesSim object.")
  }
  if (length(sims) != length(weights)) {
    stop("'sims' and 'weights' must have the same length.")
  }
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) {
    stop("'weights' must be finite and non-negative, with a positive sum.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("'weights' must sum to 1.")
  }
  weights <- weights / sum(weights)
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
      n_draws != as.integer(n_draws)) {
    stop("'n_draws' must be a positive integer.")
  }
  n_draws <- as.integer(n_draws)

  # Get common markets
  all_markets <- Reduce(intersect, lapply(sims, function(s) s@market_ids))
  if (length(all_markets) == 0L) {
    stop("No common markets across simulation objects.")
  }

  # Largest-remainder allocation preserves zero weights and returns exactly
  # n_draws, unlike independent rounding with a minimum of one per model.
  raw_counts <- weights * n_draws
  n_per_model <- floor(raw_counts)
  remainder <- n_draws - sum(n_per_model)
  if (remainder > 0L) {
    take_extra <- order(raw_counts - n_per_model, decreasing = TRUE)[
      seq_len(remainder)
    ]
    n_per_model[take_extra] <- n_per_model[take_extra] + 1L
  }
  n_per_model <- as.integer(n_per_model)

  # Draw indices are selected once per source model and reused in every
  # market, preserving the joint posterior draw across markets.
  selected <- vector("list", length(sims))
  for (i in seq_along(sims)) {
    if (n_per_model[i] == 0L) next
    common_valid <- Reduce(intersect, lapply(all_markets, function(mk) {
      which(!vapply(sims[[i]]@sims[[mk]], is.null, logical(1)))
    }))
    if (length(common_valid) == 0L) {
      stop("Simulation model ", i,
           " has no draw that is valid in every common market.")
    }
    selected[[i]] <- sample(common_valid, n_per_model[i], replace = TRUE)
  }
  source_order <- rep(seq_along(sims), n_per_model)
  source_position <- stats::ave(source_order, source_order, FUN = seq_along)
  permutation <- sample.int(n_draws)
  source_order <- source_order[permutation]
  source_position <- source_position[permutation]

  pooled_sims     <- list()
  pooled_ppre     <- list()
  pooled_ppost    <- list()
  pooled_pchg     <- list()
  pooled_spre     <- list()
  pooled_spost    <- list()
  pooled_labels   <- list()
  pooled_n_valid  <- integer(length(all_markets))
  names(pooled_n_valid) <- all_markets

  for (mk in all_markets) {
    mk_sims  <- list()
    mk_ppre  <- NULL
    mk_ppost <- NULL
    mk_spre  <- NULL
    mk_spost <- NULL

    reference_labels <- sims[[1L]]@product_labels[[mk]]
    for (i in seq_along(sims)) {
      if (!identical(sims[[i]]@product_labels[[mk]], reference_labels)) {
        stop("Product labels or ordering differ across models in market '",
             mk, "'.")
      }
    }

    for (row in seq_len(n_draws)) {
      i <- source_order[row]
      j <- selected[[i]][source_position[row]]
      s <- sims[[i]]
      mk_sims[[row]] <- s@sims[[mk]][[j]]
      mk_ppre <- rbind(mk_ppre, s@price_pre[[mk]][j, , drop = FALSE])
      mk_ppost <- rbind(mk_ppost, s@price_post[[mk]][j, , drop = FALSE])
      mk_spre <- rbind(mk_spre, s@shares_pre[[mk]][j, , drop = FALSE])
      mk_spost <- rbind(mk_spost, s@shares_post[[mk]][j, , drop = FALSE])
    }

    pooled_sims[[mk]]    <- mk_sims
    pooled_ppre[[mk]]    <- mk_ppre
    pooled_ppost[[mk]]   <- mk_ppost
    pooled_pchg[[mk]]    <- mk_ppost / mk_ppre - 1
    pooled_spre[[mk]]    <- mk_spre
    pooled_spost[[mk]]   <- mk_spost
    pooled_labels[[mk]]  <- reference_labels
    pooled_n_valid[mk]   <- n_draws
  }

  pooled_metadata <- list(
    source_models = names(sims),
    weights = weights,
    component_metadata = lapply(sims, slot, "metadata")
  )

  new("BayesSim",
    sims           = pooled_sims,
    price_pre      = pooled_ppre,
    price_post     = pooled_ppost,
    price_change   = pooled_pchg,
    shares_pre     = pooled_spre,
    shares_post    = pooled_spost,
    market_ids     = all_markets,
    product_labels = pooled_labels,
    n_draws        = n_draws,
    n_valid        = pooled_n_valid,
    failures       = stats::setNames(
      replicate(length(all_markets), rep(NA_character_, n_draws),
                simplify = FALSE),
      all_markets
    ),
    call           = match.call(),
    metadata       = pooled_metadata
  )
}
