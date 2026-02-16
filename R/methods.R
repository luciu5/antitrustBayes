#' @include AllClasses.R
NULL

# ---- BayesLogit methods ----

#' @rdname BayesLogit-class
#' @exportMethod show
setMethod("show", "BayesLogit", function(object) {
  supply_names <- c("Bertrand", "Monopolistic", "Cournot", "Competitive")
  cat("Bayesian Logit Demand Estimation\n")
  cat("  Supply model:", supply_names[object@supply_model], "\n")
  cat("  Observations:", nrow(object@data), "\n")
  cat("  Markets:", length(unique(object@data$market)), "\n")
  cat("  Groups:", length(unique(object@data$group)), "\n")
  cat("  Price SD (scaling):", round(object@price_sd, 4), "\n")

  # Basic MCMC diagnostics
  fit_summary <- rstan::summary(object@stanfit,
                                pars = c("mu_log_a", "sigma_share",
                                         "sigma_margin"))$summary
  cat("\nKey parameter estimates:\n")
  print(round(fit_summary[, c("mean", "sd", "2.5%", "97.5%", "Rhat",
                               "n_eff")], 3))
  invisible(object)
})

#' @rdname BayesLogit-class
#' @param object A \code{BayesLogit} object.
#' @exportMethod summary
setMethod("summary", "BayesLogit", function(object) {
  rstan::summary(object@stanfit)$summary
})

# ---- BayesSim methods ----

#' @rdname BayesSim-class
#' @exportMethod show
setMethod("show", "BayesSim", function(object) {
  cat("Bayesian Merger Simulation Results\n")
  cat("  Markets:", length(object@market_ids), "\n")
  cat("  Posterior draws:", object@n_draws, "\n")
  cat("  Valid simulations per market:\n")

  for (mk in object@market_ids) {
    n_ok <- object@n_valid[mk]
    cat("    Market", mk, ":", n_ok, "/", object@n_draws, "\n")

    if (n_ok > 0L) {
      pchg <- object@price_change[[mk]]
      valid_rows <- which(rowSums(is.na(pchg)) == 0)
      if (length(valid_rows) > 0L) {
        pchg_valid <- pchg[valid_rows, , drop = FALSE]
        means <- colMeans(pchg_valid)
        sds   <- apply(pchg_valid, 2, sd)
        tbl <- data.frame(
          Product = object@product_labels[[mk]],
          Mean    = sprintf("%.4f", means),
          SD      = sprintf("%.4f", sds),
          stringsAsFactors = FALSE
        )
        cat("    Price changes (proportional):\n")
        print(tbl, row.names = FALSE, right = FALSE)
      }
    }
    cat("\n")
  }
  invisible(object)
})

#' @rdname BayesSim-class
#' @param object A \code{BayesSim} object.
#' @exportMethod summary
setMethod("summary", "BayesSim", function(object) {
  results <- list()

  for (mk in object@market_ids) {
    n_ok <- object@n_valid[mk]
    if (n_ok == 0L) next

    pchg <- object@price_change[[mk]]
    valid_rows <- which(rowSums(is.na(pchg)) == 0)
    if (length(valid_rows) == 0L) next

    pchg_valid <- pchg[valid_rows, , drop = FALSE]
    labs <- object@product_labels[[mk]]

    tbl <- data.frame(
      market  = mk,
      product = labs,
      mean    = colMeans(pchg_valid),
      sd      = apply(pchg_valid, 2, sd),
      q025    = apply(pchg_valid, 2, quantile, probs = 0.025),
      q250    = apply(pchg_valid, 2, quantile, probs = 0.25),
      median  = apply(pchg_valid, 2, quantile, probs = 0.5),
      q750    = apply(pchg_valid, 2, quantile, probs = 0.75),
      q975    = apply(pchg_valid, 2, quantile, probs = 0.975),
      prob_positive = colMeans(pchg_valid > 0),
      n_valid = n_ok,
      stringsAsFactors = FALSE
    )
    rownames(tbl) <- NULL
    results[[mk]] <- tbl
  }

  do.call(rbind, results)
})
