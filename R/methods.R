#' @include AllClasses.R
NULL

# ---- Internal helper: apply antitrust method across posterior draws ----

#' Apply an antitrust method across valid posterior draws
#'
#' @param object A \code{BayesSim} object.
#' @param fn A function taking a single antitrust sim object and returning
#'   a numeric result (vector, matrix, or scalar).
#' @param market Which market(s) to compute for. Default \code{NULL} (all).
#' @return A list (one per market) of lists of per-draw results.
#' @keywords internal
.apply_to_sims <- function(object, fn, market = NULL) {
  mks <- object@market_ids
  if (!is.null(market)) {
    mks <- intersect(as.character(market), mks)
  }

  out <- list()
  for (mk in mks) {
    sims_mk <- object@sims[[mk]]
    valid <- which(!vapply(sims_mk, is.null, logical(1)))
    draw_results <- vector("list", length(valid))
    for (j in seq_along(valid)) {
      draw_results[[j]] <- tryCatch(fn(.terminal_model(sims_mk[[valid[j]]])),
        error = function(e) NULL
      )
    }
    draw_results <- draw_results[!vapply(draw_results, is.null, logical(1))]
    out[[mk]] <- draw_results
  }
  out
}

#' Summarise per-draw results into a data frame
#'
#' Takes the output of \code{.apply_to_sims} where each draw result is a
#' named numeric vector, and produces a summary data frame.
#' @param draw_list List of numeric vectors (one per draw).
#' @param posterior Logical. If \code{TRUE}, return the full matrix of draws.
#' @return A data frame with columns mean, sd, q025, median, q975 (or a
#'   matrix if \code{posterior = TRUE}).
#' @keywords internal
.summarise_vector_draws <- function(draw_list, posterior = FALSE) {
  mat <- do.call(rbind, draw_list)
  if (posterior) {
    return(mat)
  }
  data.frame(
    mean = colMeans(mat, na.rm = TRUE),
    sd = apply(mat, 2, sd, na.rm = TRUE),
    q025 = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    median = apply(mat, 2, quantile, probs = 0.5, na.rm = TRUE),
    q975 = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Summarise per-draw matrix results (e.g. elasticity matrices)
#'
#' @param draw_list List of matrices (one per draw).
#' @param posterior Logical. If \code{TRUE}, return a 3D array.
#' @return A list with \code{mean} and \code{sd} matrices, or a 3D array
#'   if \code{posterior = TRUE}.
#' @keywords internal
.summarise_matrix_draws <- function(draw_list, posterior = FALSE) {
  n <- length(draw_list)
  nr <- nrow(draw_list[[1]])
  nc <- ncol(draw_list[[1]])
  rn <- rownames(draw_list[[1]])
  cn <- colnames(draw_list[[1]])

  arr <- array(NA_real_,
    dim = c(nr, nc, n),
    dimnames = list(rn, cn, NULL)
  )
  for (i in seq_len(n)) arr[, , i] <- draw_list[[i]]

  if (posterior) {
    return(arr)
  }
  list(
    mean = apply(arr, c(1, 2), mean, na.rm = TRUE),
    sd   = apply(arr, c(1, 2), sd, na.rm = TRUE)
  )
}


# ---- Willingness-to-Pay Parameter Conversion ----

#' Convert Between Logit Alpha and CES Gamma
#'
#' Transforms posterior draws of the price-sensitivity parameter between
#' logit (alpha) and CES (gamma) scales using the Frechet--Gumbel
#' distributional equivalence.
#'
#' CES demand arises from Frechet-distributed taste shocks with shape
#' parameter \eqn{\sigma = \gamma - 1}, while logit demand arises from
#' Gumbel-distributed shocks. Since CES operates on log-prices and logit
#' on price levels, the conversion at reference price \eqn{\bar{p}} is:
#' \deqn{\alpha = (\gamma - 1) / \bar{p}, \quad \gamma = 1 + \alpha \bar{p}}
#'
#' @param object A \code{\linkS4class{BayesLogit}} object.
#' @param to_gamma Logical. If \code{TRUE} (default), convert logit alpha
#'   draws to CES gamma equivalents. If \code{FALSE}, convert CES gamma
#'   draws to logit alpha equivalents.
#' @return A numeric vector of converted parameter draws (one per posterior
#'   sample).
#'
#' @examples
#' \donttest{
#' # Convert CES gamma draws to logit alpha scale
#' alpha_equiv <- wtp_convert(ces_fit, to_gamma = FALSE)
#'
#' # Convert logit alpha draws to CES gamma scale
#' gamma_equiv <- wtp_convert(logit_fit, to_gamma = TRUE)
#' }
#'
#' @export
setGeneric("wtp_convert", function(object, to_gamma = TRUE)
  standardGeneric("wtp_convert"))

#' @rdname wtp_convert
#' @exportMethod wtp_convert
setMethod("wtp_convert", "BayesLogit", function(object, to_gamma = TRUE) {
  p_bar <- mean(object@data$prices)
  mu_log_a_draws <- as.vector(
    .bayes_backend_extract(object@stanfit, "mu_log_a")[[1]]
  )

  if (to_gamma) {
    if (object@demand_model != "logit") {
      stop("to_gamma = TRUE requires a logit fit.")
    }
    # alpha = exp(mu_log_a) / price_sd  =>  gamma = 1 + alpha * p_bar
    alpha_draws <- exp(mu_log_a_draws) / object@price_sd
    gamma_draws <- 1 + alpha_draws * p_bar
    gamma_draws
  } else {
    if (object@demand_model != "ces") {
      stop("to_gamma = FALSE requires a CES fit.")
    }
    # gamma = 1 + exp(mu_log_a) / log_price_sd  =>  alpha = (gamma - 1) / p_bar
    gamma_draws <- 1 + exp(mu_log_a_draws) / object@price_sd
    alpha_draws <- (gamma_draws - 1) / p_bar
    alpha_draws
  }
})


# ---- BayesLogit methods ----

#' @rdname BayesFit-class
#' @exportMethod show
setMethod("show", "BayesLogit", function(object) {
  supply_names <- c(
    "Bertrand", "Monopolistic", "Cournot", "Auction", "Bargaining"
  )
  demand_label <- if (length(object@demand_model) > 0 &&
    object@demand_model == "ces") {
    "CES"
  } else {
    "logit"
  }
  cat("Bayesian", demand_label, "Demand Estimation\n")
  if (!is.null(object@metadata$backend))
    cat("  Sampler backend:", object@metadata$backend, "\n")
  cat("  Supply model:", supply_names[object@supply_model], "\n")
  cat("  Observations:", nrow(object@data), "\n")
  cat("  Markets:", length(unique(object@data$market)), "\n")
  cat("  Groups:", length(unique(object@data$group)), "\n")
  cat("  Quality scope:", .quality_scope_from_fit(object, strict = FALSE), "\n")
  intercept_mode <- .market_intercept_parameterization_from_fit(object)
  if (identical(intercept_mode, "collapsed"))
    cat("  Market intercept: collapsed Gaussian integration\n")
  if (!is.null(object@stan_data$N_quality))
    cat("  Persistent quality units:", object@stan_data$N_quality, "\n")
  price_label <- if (demand_label == "CES") "Log-price SD" else "Price SD"
  cat(" ", price_label, "(scaling):", round(object@price_sd, 4), "\n")

  # Basic MCMC diagnostics
  fit_summary <- .bayes_backend_summary(object@stanfit,
    variables = c("mu_log_a", "sigma_share", "sigma_margin")
  )
  cat("\nKey parameter estimates:\n")
  print(fit_summary[, c("variable", "mean", "sd", "q5", "q95",
                         "rhat", "ess_bulk")])
  if (identical(object@metadata$inference, "pathfinder")) {
    cat("  Inference: exploratory Pathfinder (Pareto k = ",
        format(object@metadata$pathfinder_pareto_k), ")\n", sep = "")
  }

  # For CES fits, show equivalent logit alpha
  if (demand_label == "CES") {
    alpha_equiv <- wtp_convert(object, to_gamma = FALSE)
    cat("\nEquivalent logit alpha (Frechet-Gumbel):\n")
    cat("  mean =", round(mean(alpha_equiv), 4),
        " sd =", round(sd(alpha_equiv), 4),
        " 95% CI = [", round(quantile(alpha_equiv, 0.025), 4), ",",
        round(quantile(alpha_equiv, 0.975), 4), "]\n")
  }
  invisible(object)
})

#' @rdname BayesFit-class
#' @param object A \code{BayesLogit} object.
#' @exportMethod summary
setMethod("summary", "BayesLogit", function(object) {
  s <- .bayes_backend_summary(object@stanfit)
  print(s)
  invisible(s)
})


# ---- BayesSim methods: elast, diversion, CV, HypoMonTest ----

#' Posterior Distribution of Elasticities
#'
#' Computes own- and cross-price elasticities for each posterior draw stored
#' in a \code{\linkS4class{BayesSim}} object.
#'
#' @param object A \code{BayesSim} object.
#' @param preMerger Logical. If \code{TRUE} (default), compute pre-merger
#'   elasticities; if \code{FALSE}, post-merger.
#' @param market Which market(s) to compute for. Default \code{NULL} (first).
#' @param posterior Logical. If \code{TRUE}, return the full 3D array
#'   (products x products x draws). If \code{FALSE} (default), return a
#'   list with \code{mean} and \code{sd} matrices.
#'
#' @return If a single market: a list with \code{mean} and \code{sd}
#'   matrices (or a 3D array if \code{posterior = TRUE}).
#'   If multiple markets: a named list of such results.
#'
#' @seealso \code{\link[antitrust]{elast}}
#' @rdname BayesSim-methods
#' @exportMethod elast
setMethod(
  "elast", "BayesSim",
  function(object, preMerger = TRUE, market = NULL, posterior = FALSE) {
    if (is.null(market)) market <- object@market_ids[1]
    fn <- function(sim) elast(sim, preMerger = preMerger)
    draws <- .apply_to_sims(object, fn, market = market)

    out <- lapply(draws, function(dl) {
      if (length(dl) == 0L) {
        return(NULL)
      }
      .summarise_matrix_draws(dl, posterior = posterior)
    })
    if (length(out) == 1L) out[[1]] else out
  }
)

#' Posterior Distribution of Diversion Ratios
#'
#' Computes diversion ratios for each posterior draw stored in a
#' \code{\linkS4class{BayesSim}} object.
#'
#' @param object A \code{BayesSim} object.
#' @param preMerger Logical. If \code{TRUE} (default), compute pre-merger
#'   diversions.
#' @param revenue Logical. If \code{TRUE}, compute revenue diversion.
#'   Default \code{FALSE}.
#' @param market Which market(s). Default \code{NULL} (first).
#' @param posterior Logical. If \code{TRUE}, return the full 3D array.
#'
#' @return Same structure as \code{elast}.
#'
#' @seealso \code{\link[antitrust]{diversion}}
#' @rdname BayesSim-methods
#' @exportMethod diversion
setMethod(
  "diversion", "BayesSim",
  function(object, preMerger = TRUE, revenue = FALSE,
           market = NULL, posterior = FALSE) {
    if (is.null(market)) market <- object@market_ids[1]
    fn <- function(sim) {
      diversion(sim,
        preMerger = preMerger,
        revenue = revenue
      )
    }
    draws <- .apply_to_sims(object, fn, market = market)

    out <- lapply(draws, function(dl) {
      if (length(dl) == 0L) {
        return(NULL)
      }
      .summarise_matrix_draws(dl, posterior = posterior)
    })
    if (length(out) == 1L) out[[1]] else out
  }
)

#' Posterior Distribution of Compensating Variation
#'
#' Computes consumer welfare (compensating variation) for each posterior draw
#' stored in a \code{\linkS4class{BayesSim}} object.
#'
#' @param object A \code{BayesSim} object.
#' @param market Which market(s). Default \code{NULL} (all).
#' @param posterior Logical. If \code{TRUE}, return the full vector of draws.
#'
#' @return If a single market: a named numeric with \code{mean}, \code{sd},
#'   \code{q025}, \code{median}, \code{q975} (or a numeric vector of draws
#'   if \code{posterior = TRUE}).
#'   If multiple markets: a named list or data frame of such results.
#'
#' @seealso \code{\link[antitrust]{CV}}
#' @rdname BayesSim-methods
#' @exportMethod CV
setMethod(
  "CV", "BayesSim",
  function(object, market = NULL, posterior = FALSE) {
    fn <- function(sim) CV(sim)
    draws <- .apply_to_sims(object, fn, market = market)

    out <- lapply(draws, function(dl) {
      if (length(dl) == 0L) {
        return(NULL)
      }
      vals <- unlist(dl)
      if (posterior) {
        return(vals)
      }
      c(
        mean = mean(vals, na.rm = TRUE),
        sd = sd(vals, na.rm = TRUE),
        q025 = unname(quantile(vals, 0.025, na.rm = TRUE)),
        median = unname(quantile(vals, 0.5, na.rm = TRUE)),
        q975 = unname(quantile(vals, 0.975, na.rm = TRUE))
      )
    })
    if (length(out) == 1L) out[[1]] else out
  }
)

#' Posterior Distribution of Hypothetical Monopolist Test
#'
#' Evaluates the Hypothetical Monopolist Test for each posterior draw stored
#' in a \code{\linkS4class{BayesSim}} object, returning the proportion of
#' draws for which the test passes.
#'
#' @param object A \code{BayesSim} object.
#' @param prodIndex Numeric vector of product indices forming the candidate
#'   market.
#' @param ssnip SSNIP threshold. Default 0.05.
#' @param market Which market(s). Default \code{NULL} (first).
#' @param posterior Logical. If \code{TRUE}, return the full logical vector
#'   of per-draw results.
#' @param ... Additional arguments passed to
#'   \code{\link[antitrust]{HypoMonTest}}.
#'
#' @return If \code{posterior = FALSE} (default): proportion of draws where
#'   the HMT passes (i.e., the probability that a SSNIP is profitable).
#'   If \code{posterior = TRUE}: a logical vector of per-draw results.
#'
#' @seealso \code{\link[antitrust]{HypoMonTest}}
#' @rdname BayesSim-methods
#' @exportMethod HypoMonTest
setMethod(
  "HypoMonTest", "BayesSim",
  function(object, prodIndex, ssnip = 0.05,
           market = NULL, posterior = FALSE, ...) {
    if (is.null(market)) market <- object@market_ids[1]
    fn <- function(sim) {
      HypoMonTest(sim,
        prodIndex = prodIndex,
        ssnip = ssnip, ...
      )
    }
    draws <- .apply_to_sims(object, fn, market = market)

    out <- lapply(draws, function(dl) {
      if (length(dl) == 0L) {
        return(NULL)
      }
      vals <- unlist(dl)
      if (posterior) {
        return(vals)
      }
      mean(vals, na.rm = TRUE)
    })
    if (length(out) == 1L) out[[1]] else out
  }
)


# ---- BayesSim: show and summary ----

#' @rdname BayesSim-class
#' @exportMethod show
setMethod("show", "BayesSim", function(object) {
  cat("Bayesian Counterfactual Simulation Results\n")
  cat("  Markets:", length(object@market_ids), "\n")
  cat("  Posterior draws:", object@n_draws, "\n")
  cat("  Valid simulations per market:\n")

  for (mk in object@market_ids) {
    n_ok <- object@n_valid[mk]
    cat("    Market", mk, ":", n_ok, "/", object@n_draws, "\n")
  }
  cat("\n")

  # Print market-level summary (means only)
  s <- summary(object, market = TRUE)
  if (!is.null(s)) {
    cat("Market-level summary (posterior means):\n")
    print(s$mean, digits = 3, row.names = FALSE)
  }
  invisible(object)
})


.bayes_sim_draw_summary <- function(model, market, levels, revenue, ...) {
  invisible(capture.output(res <- summary(.terminal_model(model),
    market = market, levels = levels, revenue = revenue, digits = 10, ...)))
  res
}

#' Summarise Bayesian Counterfactual Simulation Results
#'
#' Mirrors the \code{summary} method for \pkg{antitrust} \code{Bertrand}
#' objects, but computes posterior means and standard deviations across
#' all valid draws.
#'
#' @param object A \code{BayesSim} object.
#' @param market Logical. If \code{TRUE}, return market-level statistics
#'   (HHI change, industry price change, compensating variation, etc.).
#'   If \code{FALSE} (default), return product-level statistics
#'   (pre/post prices, shares, price deltas).
#' @param levels Logical. If \code{TRUE}, report level changes instead of
#'   percent changes. Default \code{FALSE}.
#' @param revenue Logical. If \code{TRUE}, report revenue shares.
#'   Default \code{TRUE}.
#' @param digits Integer; number of significant digits. Default 2.
#' @param ... Additional arguments passed to the underlying
#'   \pkg{antitrust} \code{summary} method.
#'
#' @return A list with components \code{mean} and \code{sd}, each a
#'   data frame with the same structure as the \pkg{antitrust}
#'   \code{summary} output. If \code{market = TRUE}, each is a
#'   single-row data frame of market-level statistics; if
#'   \code{market = FALSE}, each has one row per product.
#'   The \code{diagnostics} attribute records requested draw IDs, summary
#'   failures, and valid/missing counts by estimand. Attrition emits a warning.
#'
#' @rdname BayesSim-class
#' @exportMethod summary
setMethod(
  "summary", "BayesSim",
  function(object, market = FALSE, levels = FALSE, revenue = TRUE,
           digits = 2, ...) {
    # Collect summary data frames from each valid draw
    results_by_market <- list()
    diagnostics <- list()

    for (mk in object@market_ids) {
      sims_mk <- object@sims[[mk]]
      valid <- which(!vapply(sims_mk, is.null, logical(1)))
      ids <- object@metadata$draw_indices
      if (is.null(ids)) ids <- seq_along(sims_mk)
      errors <- rep(NA_character_, length(sims_mk))
      errors[setdiff(seq_along(sims_mk), valid)] <- "simulation_unavailable"
      diagnostics[[mk]] <- list(requested = length(sims_mk), draw_ids = ids,
                                errors = errors, summarized = 0L, estimands = NULL)
      if (length(valid) == 0L) next

      draw_dfs <- vector("list", length(valid))
      for (j in seq_along(valid)) {
        value <- tryCatch(.bayes_sim_draw_summary(sims_mk[[valid[j]]],
          market, levels, revenue, ...), error = identity)
        if (inherits(value, "error") || !is.data.frame(value)) {
          diagnostics[[mk]]$errors[valid[j]] <- if (inherits(value, "error"))
            conditionMessage(value) else "summary_unavailable"
          value <- NULL
        }
        draw_dfs[j] <- list(value)
      }
      keep <- !vapply(draw_dfs, is.null, logical(1))
      diagnostics[[mk]]$summary_draw_ids <- ids[valid[keep]]
      diagnostics[[mk]]$summarized <- sum(keep)
      draw_dfs <- draw_dfs[!vapply(draw_dfs, is.null, logical(1))]
      if (length(draw_dfs) == 0L) next

      results_by_market[[mk]] <- draw_dfs
    }

    if (length(results_by_market) == 0L) {
      reasons <- unlist(lapply(diagnostics, `[[`, "errors"), use.names = FALSE)
      counts <- table(reasons, useNA = "no")
      condition <- structure(list(message = paste0(
        "No simulation draws could be summarized (", length(reasons), " requested): ",
        paste(names(counts), counts, sep = "=", collapse = "; ")),
        call = NULL, diagnostics = diagnostics),
        class = c("bayes_summary_attrition", "warning", "condition"))
      warning(condition)
      return(NULL)
    }

    # Aggregate: compute mean and sd of each numeric column across draws
    .aggregate_draw_dfs <- function(draw_dfs, mk) {
      ref <- draw_dfs[[1]]
      num_cols <- vapply(ref, is.numeric, logical(1))

      # For non-market summaries, drop the isParty column if present
      if ("isParty" %in% names(ref)) {
        num_cols["isParty"] <- FALSE
      }

      num_names <- names(which(num_cols))
      n <- length(draw_dfs)

      # Stack numeric columns into a 3D array: products x columns x draws
      nr <- nrow(ref)
      nc <- length(num_names)
      arr <- array(NA_real_, dim = c(nr, nc, n))
      for (i in seq_len(n)) {
        for (k in seq_len(nc)) {
          arr[, k, i] <- as.numeric(draw_dfs[[i]][[num_names[k]]])
        }
      }

      finite <- is.finite(arr)
      counts <- expand.grid(row = seq_len(nr), estimand = num_names,
                            stringsAsFactors = FALSE)
      counts$valid <- as.vector(apply(finite, c(1, 2), sum))
      counts$missing <- diagnostics[[mk]]$requested - counts$valid
      diagnostics[[mk]]$estimands <<- counts
      bad <- which(!finite, arr.ind = TRUE)
      diagnostics[[mk]]$nonfinite_values <<- data.frame(
        row = bad[, 1L], estimand = num_names[bad[, 2L]],
        draw_id = diagnostics[[mk]]$summary_draw_ids[bad[, 3L]])
      arr[!finite] <- NA_real_
      mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
      sd_mat <- apply(arr, c(1, 2), sd, na.rm = TRUE)

      # Build output data frames matching antitrust summary format
      mean_df <- ref
      sd_df <- ref
      for (k in seq_len(nc)) {
        mean_df[[num_names[k]]] <- round(mean_mat[, k], digits)
        sd_df[[num_names[k]]] <- round(sd_mat[, k], digits)
      }

      list(mean = mean_df, sd = sd_df)
    }

    # Process each market
    out_list <- lapply(names(results_by_market), function(mk)
      .aggregate_draw_dfs(results_by_market[[mk]], mk))
    names(out_list) <- names(results_by_market)

    result <- if (length(out_list) == 1L) out_list[[1]] else out_list
    attr(result, "diagnostics") <- diagnostics
    if (any(vapply(diagnostics, function(d) any(!is.na(d$errors)) ||
                   any(d$estimands$missing > 0L), logical(1))))
      warning("Simulation summary omitted draws or values; see attr(result, 'diagnostics').",
              call. = FALSE)

    # Print output (needed for Rscript where auto-printing doesn't occur)
    if (length(out_list) == 1L) {
      cat("Posterior Mean:\n")
      print(result$mean, digits = digits, row.names = !market)
      cat("\nPosterior SD:\n")
      print(result$sd, digits = digits, row.names = !market)
    } else {
      for (mk in names(result)) {
        cat("Market:", mk, "\n")
        cat("  Posterior Mean:\n")
        print(result[[mk]]$mean, digits = digits, row.names = !market)
        cat("  Posterior SD:\n")
        print(result[[mk]]$sd, digits = digits, row.names = !market)
        cat("\n")
      }
    }

    invisible(result)
  }
)
