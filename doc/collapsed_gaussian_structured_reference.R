## Independent structured reference for a centered, sparse-loading Gaussian.
##
## The production package is deliberately not sourced here.  For observations
## l = 1,...,L, the model is
##
##   z_l = x_l b_{g_l} + eps_l,
##   eps_l ~ N(0, v_l),       b ~ N(0, tau^2 P_G),
##   P_G = I_G - 1_G 1_G' / G.
##
## The structured evaluator uses only group sufficient statistics and is
## O(L + G).  The dense evaluator is an oracle for small matrices only.  This
## is an algorithmic reference and benchmark, not a Stan performance claim.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

stop_if_not_scalar <- function(x, name) {
  if (length(x) != 1L || !is.finite(x)) {
    stop(name, " must be one finite scalar")
  }
}

group_sum <- function(values, group, G) {
  ## tabulate() has no weighted form in base R.  This loop has the same
  ## linear complexity and retains zeroes for explicitly empty groups.
  out <- numeric(G)
  for (i in seq_along(values)) out[group[i]] <- out[group[i]] + values[i]
  out
}

validate_inputs <- function(z, v, group, x, tau, G) {
  if (length(z) == 0L) stop("at least one observation is required")
  if (length(v) != length(z) || length(group) != length(z) ||
      length(x) != length(z)) {
    stop("z, v, group, and x must have equal lengths")
  }
  if (any(!is.finite(z)) || any(!is.finite(v)) || any(v <= 0) ||
      any(!is.finite(x))) {
    stop("z, x, and v must be finite, with v > 0")
  }
  if (any(group != as.integer(group)) || any(group < 1L) ||
      any(group > G)) {
    stop("group must contain integer labels in 1,...,G")
  }
  stop_if_not_scalar(tau, "tau")
  if (tau < 0) stop("tau must be nonnegative")
  stop_if_not_scalar(G, "G")
  if (G < 1 || G != as.integer(G)) stop("G must be a positive integer")
  invisible(NULL)
}

## Exact O(L + G) normalized Gaussian log density.
structured_gaussian_log_density <- function(z, v, group, x, tau,
                                            G = max(group)) {
  G <- as.integer(G)
  validate_inputs(z, v, group, x, tau, G)
  L <- length(z)
  tau2 <- tau * tau
  inv_v <- 1 / v
  a <- group_sum(x * x * inv_v, group, G)
  c_stat <- group_sum(x * z * inv_v, group, G)
  e <- sum(z * z * inv_v)
  group_mean <- numeric(G)
  positive <- a > 0
  group_mean[positive] <- c_stat[positive] / a[positive]
  ## Independent two-pass residual calculation. Production uses one-pass
  ## weighted regression updates; neither subtracts nearly equal large sums.
  within <- sum((z - x * group_mean[group])^2 * inv_v)

  w <- if (tau == 0) rep(1, G) else 1 / (1 + tau2 * a)
  d <- mean(w)
  bsum <- sum(c_stat * w)

  ## For G=1, P_G is exactly zero.  The explicit branch avoids cancellation
  ## between the two centered-rank-one terms and is also the tau=0 limit.
  if (G == 1L || tau == 0) {
    logdet <- sum(log(v))
    quadratic <- e
    return(list(
      log_density = -0.5 * (L * log(2 * pi) + logdet + quadratic),
      logdet = logdet,
      quadratic = quadratic,
      a = a,
      c = c_stat,
      e = e,
      w = w,
      d = d,
      bsum = bsum
    ))
  }

  sum_w <- sum(w)
  logdet_V <- sum(log(v)) + sum(log1p(tau2 * a))
  q0 <- within + sum(group_mean^2 * a * w)
  quadratic <- q0 + tau2 * bsum * bsum / sum_w
  logdet <- logdet_V + log(d)
  list(
    log_density = -0.5 * (L * log(2 * pi) + logdet + quadratic),
    logdet = logdet,
    quadratic = quadratic,
    a = a,
    c = c_stat,
    e = e,
    w = w,
    d = d,
    bsum = bsum
  )
}

make_loading_matrix <- function(group, x, G) {
  X <- matrix(0, nrow = length(group), ncol = G)
  X[cbind(seq_along(group), as.integer(group))] <- x
  X
}

log_density_from_covariance <- function(z, covariance) {
  R <- chol(covariance)
  solved <- forwardsolve(t(R), z)
  logdet <- 2 * sum(log(diag(R)))
  quadratic <- sum(solved * solved)
  list(
    log_density = -0.5 * (length(z) * log(2 * pi) + logdet + quadratic),
    logdet = logdet,
    quadratic = quadratic
  )
}

## Dense K = diag(v) + tau^2 X P X' oracle.  Use only for small N.
dense_gaussian_log_density <- function(z, v, group, x, tau,
                                       G = max(group)) {
  G <- as.integer(G)
  validate_inputs(z, v, group, x, tau, G)
  X <- make_loading_matrix(group, x, G)
  P <- diag(G) - matrix(1 / G, nrow = G, ncol = G)
  covariance <- diag(v, nrow = length(z), ncol = length(z)) +
    tau * tau * X %*% P %*% t(X)
  log_density_from_covariance(z, covariance)
}

rho_augmented_rows <- function(z_share, group_share, margin_residual,
                               matched_rows, sigma_share, sigma_margin,
                               rho) {
  L <- length(z_share)
  if (length(group_share) != L) stop("share groups must match z_share")
  if (length(margin_residual) != length(matched_rows)) {
    stop("margin residuals must match matched_rows")
  }
  if (any(matched_rows < 1L) || any(matched_rows > L) ||
      any(matched_rows != as.integer(matched_rows))) {
    stop("matched_rows must be valid integer share-row indices")
  }
  stop_if_not_scalar(sigma_share, "sigma_share")
  stop_if_not_scalar(sigma_margin, "sigma_margin")
  stop_if_not_scalar(rho, "rho")
  if (sigma_share <= 0 || sigma_margin <= 0 || abs(rho) >= 1) {
    stop("sigma values must be positive and |rho| < 1")
  }
  slope <- rho * sigma_margin / sigma_share
  ## The transformed margin is m - nu - slope * (y - mu).  The triangular
  ## transformation has unit Jacobian, so its joint density is unchanged.
  list(
    z = c(z_share, margin_residual - slope * z_share[matched_rows]),
    group = c(group_share, group_share[matched_rows]),
    x = c(rep(1, L), rep(-slope, length(matched_rows))),
    v = c(rep(sigma_share * sigma_share, L),
          rep(sigma_margin * sigma_margin * (1 - rho * rho),
              length(matched_rows))),
    slope = slope,
    matched_rows = matched_rows
  )
}

rho_structured_log_density <- function(z_share, group_share, margin_residual,
                                       matched_rows, tau, sigma_share,
                                       sigma_margin, rho, G = max(group_share)) {
  augmented <- rho_augmented_rows(
    z_share, group_share, margin_residual, matched_rows,
    sigma_share, sigma_margin, rho
  )
  out <- structured_gaussian_log_density(
    augmented$z, augmented$v, augmented$group, augmented$x, tau, G
  )
  out$augmented <- augmented
  out
}

rho_dense_original_log_density <- function(z_share, group_share,
                                           margin_residual, matched_rows,
                                           tau, sigma_share, sigma_margin,
                                           rho, G = max(group_share)) {
  L <- length(z_share)
  q <- length(margin_residual)
  if (length(group_share) != L || length(matched_rows) != q) {
    stop("rho inputs have inconsistent lengths")
  }
  if (any(matched_rows < 1L) || any(matched_rows > L) ||
      any(matched_rows != as.integer(matched_rows))) {
    stop("matched_rows must be valid integer share-row indices")
  }
  stop_if_not_scalar(sigma_share, "sigma_share")
  stop_if_not_scalar(sigma_margin, "sigma_margin")
  stop_if_not_scalar(rho, "rho")
  if (sigma_share <= 0 || sigma_margin <= 0 || abs(rho) >= 1) {
    stop("sigma values must be positive and |rho| < 1")
  }
  X_share <- make_loading_matrix(group_share, rep(1, L), G)
  P <- diag(G) - matrix(1 / G, nrow = G, ncol = G)
  share_covariance <- diag(sigma_share * sigma_share, nrow = L, ncol = L) +
    tau * tau * X_share %*% P %*% t(X_share)
  R_match <- matrix(0, nrow = L, ncol = q)
  if (q) R_match[cbind(matched_rows, seq_len(q))] <- 1
  cross_covariance <- rho * sigma_share * sigma_margin * R_match
  margin_covariance <- diag(sigma_margin * sigma_margin,
                            nrow = q, ncol = q)
  covariance <- rbind(
    cbind(share_covariance, cross_covariance),
    cbind(t(cross_covariance), margin_covariance)
  )
  log_density_from_covariance(c(z_share, margin_residual), covariance)
}

max_component_error <- function(a, b) {
  max(abs(c(a$log_density, a$logdet, a$quadratic) -
            c(b$log_density, b$logdet, b$quadratic)))
}

run_structured_validation <- function(verbose = TRUE) {
  ## The three fixtures intentionally cover single-product firms, repeated
  ## products within groups, and exact repeated panel rows.
  fixtures <- list(
    one_product = list(
      group = seq_len(9L),
      x = rep(1, 9L),
      description = "9 one-product groups"
    ),
    multiproduct = list(
      group = rep(seq_len(4L), each = 3L),
      x = c(0.55, 1.10, -0.80, 1.30, 0.70, 0.95,
             1.25, -0.60, 0.85, 0.90, 1.40, 0.65),
      description = "4 groups, 3 products per group"
    ),
    repeated_panel = list(
      group = rep(c(1L, 1L, 2L, 3L), times = 3L),
      x = rep(c(0.80, 1.20, -0.65, 1.45), times = 3L),
      description = "4-row panel pattern repeated over 3 periods"
    )
  )
  records <- list()
  append_record <- function(check, cases, error, tolerance, notes) {
    records[[length(records) + 1L]] <<- data.frame(
      check = check, cases = cases, max_abs_error = error,
      tolerance = tolerance, status = if (error <= tolerance) "passed" else "failed",
      notes = notes
    )
  }

  errors <- c(log_density = 0, logdet = 0, quadratic = 0)
  case_count <- 0L
  seeds <- c(101L, 211L, 307L, 401L, 503L, 601L, 709L, 809L,
             907L, 1009L, 1103L, 1201L, 1301L, 1409L, 1511L, 1601L)
  tau_values <- c(0, 0.13, 0.71, 1.90)
  for (seed in seeds) {
    set.seed(seed)
    for (fixture_name in names(fixtures)) {
      fixture <- fixtures[[fixture_name]]
      L <- length(fixture$group)
      z <- rnorm(L, sd = 1.1)
      v <- exp(runif(L, log(0.07), log(2.4)))
      ## Every case has unequal variances and non-unit sparse loadings except
      ## the one-product fixture's first base loading; all are overwritten by
      ## this deterministic weighted perturbation.
      x <- fixture$x * (0.60 + runif(L, min = 0, max = 0.90))
      for (tau in tau_values) {
        structured <- structured_gaussian_log_density(
          z, v, fixture$group, x, tau, G = max(fixture$group)
        )
        dense <- dense_gaussian_log_density(
          z, v, fixture$group, x, tau, G = max(fixture$group)
        )
        component_error <- abs(c(
          structured$log_density - dense$log_density,
          structured$logdet - dense$logdet,
          structured$quadratic - dense$quadratic
        ))
        errors <- pmax(errors, component_error)
        case_count <- case_count + 1L
      }
    }
  }
  append_record(
    "structured_vs_dense_generic", case_count, max(errors), 2e-10,
    "deterministic seeds; one-product, multiproduct, repeated-panel, unequal-v, weighted-x fixtures"
  )

  ## Explicit centered-boundary checks.  These exercise both branches that
  ## remove the P_G rank-one cancellation.
  set.seed(1717L)
  z1 <- rnorm(8L)
  v1 <- exp(runif(8L, log(0.1), log(2)))
  x1 <- runif(8L, 0.3, 1.8)
  g1 <- rep(1L, 8L)
  g1_errors <- numeric()
  for (tau in c(0, 0.2, 3.1)) {
    s <- structured_gaussian_log_density(z1, v1, g1, x1, tau, G = 1L)
    d <- dense_gaussian_log_density(z1, v1, g1, x1, tau, G = 1L)
    g1_errors <- c(g1_errors, max_component_error(s, d))
  }
  append_record(
    "G1_centered_projection_boundary", length(g1_errors), max(g1_errors), 2e-12,
    "P_1 = 0 for tau=0, positive tau, and unequal variances"
  )
  tau0_errors <- numeric()
  set.seed(1818L)
  for (fixture_name in names(fixtures)) {
    fixture <- fixtures[[fixture_name]]
    L <- length(fixture$group)
    z <- rnorm(L)
    v <- exp(runif(L, log(0.05), log(3)))
    x <- runif(L, -1.5, 1.5)
    s <- structured_gaussian_log_density(z, v, fixture$group, x, 0)
    expected <- -0.5 * (L * log(2 * pi) + sum(log(v)) + sum(z * z / v))
    tau0_errors <- c(tau0_errors, abs(s$log_density - expected))
  }
  append_record(
    "tau_zero_diagonal_limit", length(tau0_errors), max(tau0_errors), 2e-14,
    "tau=0 removes the centered quality contribution exactly"
  )

  ## The augmented rho construction is checked against the original dense
  ## share/margin block, including missing margins and the three components.
  rho_errors <- c(log_density = 0, logdet = 0, quadratic = 0)
  rho_cases <- 0L
  rho_seeds <- c(1913L, 2017L, 2111L, 2203L, 2309L, 2411L, 2503L, 2609L,
                 2707L, 2801L, 2903L, 3001L)
  for (seed in rho_seeds) {
    set.seed(seed)
    fixture <- fixtures[[if (seed %% 2L) "multiproduct" else "repeated_panel"]]
    L <- length(fixture$group)
    z_share <- rnorm(L)
    z_margin_all <- rnorm(L)
    matched_rows <- seq(2L, L, by = 2L)
    z_margin <- z_margin_all[seq_along(matched_rows)]
    sigma_share <- exp(runif(1, log(0.15), log(0.7)))
    sigma_margin <- exp(runif(1, log(0.12), log(0.8)))
    rho <- sample(c(-0.68, -0.27, 0.41, 0.79), 1L)
    for (tau in c(0, 0.37, 1.25)) {
      structured <- rho_structured_log_density(
        z_share, fixture$group, z_margin, matched_rows, tau,
        sigma_share, sigma_margin, rho, G = max(fixture$group)
      )
      dense <- rho_dense_original_log_density(
        z_share, fixture$group, z_margin, matched_rows, tau,
        sigma_share, sigma_margin, rho, G = max(fixture$group)
      )
      rho_errors <- pmax(rho_errors, abs(c(
        structured$log_density - dense$log_density,
        structured$logdet - dense$logdet,
        structured$quadratic - dense$quadratic
      )))
      rho_cases <- rho_cases + 1L
    }
  }
  append_record(
    "rho_augmented_vs_dense_original", rho_cases, max(rho_errors), 3e-10,
    "share rows plus matched observed margins; transformed margin has unit Jacobian"
  )

  validation <- do.call(rbind, records)
  if (verbose) {
    print(validation, row.names = FALSE)
    cat(sprintf("\nMaximum generic component error: %.3e\n", max(errors)))
    cat(sprintf("Maximum rho augmented component error: %.3e\n", max(rho_errors)))
    if (any(validation$status != "passed")) stop("structured validation failed")
    cat("All structured and dense reference checks passed.\n")
  }
  if (any(validation$status != "passed")) stop("structured validation failed")
  validation
}

median_elapsed <- function(fun, repetitions, batch = FALSE,
                           target_seconds = 0.05) {
  if (repetitions < 1L) {
    return(list(median_seconds = NA_real_, inner_iterations = 0L))
  }
  invisible(fun())
  inner_iterations <- 1L
  if (batch) {
    ## A single R timing tick is too coarse for the O(L+G) call at small N.
    ## Calibrate a fixed batch once, then report seconds per call.
    max_inner <- 1000000L
    repeat {
      elapsed_probe <- unname(system.time({
        for (j in seq_len(inner_iterations)) invisible(fun())
      })[["elapsed"]])
      if (elapsed_probe >= target_seconds || inner_iterations >= max_inner) break
      inner_iterations <- min(max_inner, inner_iterations * 2L)
    }
  }
  elapsed <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    elapsed[i] <- unname(system.time({
      for (j in seq_len(inner_iterations)) invisible(fun())
    })[["elapsed"]]) / inner_iterations
  }
  list(median_seconds = median(elapsed), inner_iterations = inner_iterations)
}

make_benchmark_case <- function(N, layout, seed) {
  set.seed(seed)
  if (layout == "one_product") {
    group <- seq_len(N)
    x <- rep(1, N)
    layout_note <- "G=N; one product per group"
  } else if (layout == "multiproduct_panel") {
    ## Four rows per group, corresponding to repeated products/panel rows;
    ## the final group may have fewer than four rows.
    G <- ceiling(N / 4)
    group <- rep(seq_len(G), each = 4L, length.out = N)
    x <- 0.70 + 0.10 * ((seq_len(N) - 1L) %% 11L) +
      0.03 * ((seq_len(N) - 1L) %% 3L)
    layout_note <- "G=ceiling(N/4); four repeated rows per group"
  } else {
    stop("unknown benchmark layout: ", layout)
  }
  z <- rnorm(N, sd = 0.9)
  v <- exp(runif(N, log(0.08), log(2.0)))
  tau <- 0.63
  matched_rows <- seq(2L, N, by = 2L)
  margin_residual <- rnorm(length(matched_rows), sd = 0.75)
  rho <- 0.43
  sigma_share <- 0.31
  sigma_margin <- 0.47
  augmented <- rho_augmented_rows(
    z, group, margin_residual, matched_rows,
    sigma_share, sigma_margin, rho
  )
  list(
    N = N, layout = layout, layout_note = layout_note, G = max(group),
    z = z, v = v, group = group, x = x, tau = tau,
    rho = list(
      z_share = z, group_share = group, margin_residual = margin_residual,
      matched_rows = matched_rows, sigma_share = sigma_share,
      sigma_margin = sigma_margin, rho = rho, augmented = augmented
    )
  )
}

dense_bytes <- function(n) as.double(n) * as.double(n) * 8
format_bytes <- function(bytes) {
  if (bytes < 1024) return(sprintf("%.0f B", bytes))
  if (bytes < 1024^2) return(sprintf("%.1f KiB", bytes / 1024))
  if (bytes < 1024^3) return(sprintf("%.1f MiB", bytes / 1024^2))
  sprintf("%.2f GiB", bytes / 1024^3)
}

benchmark_one_case <- function(case, rho_flag, structured_repetitions,
                               dense_repetitions, dense_limit = 1000L) {
  if (rho_flag) {
    input <- case$rho$augmented
    structured_fun <- function() structured_gaussian_log_density(
      input$z, input$v, input$group, input$x, case$tau, G = case$G
    )
    n_dense <- length(input$z)
    dense_fun <- function() rho_dense_original_log_density(
      case$rho$z_share, case$rho$group_share, case$rho$margin_residual,
      case$rho$matched_rows, case$tau, case$rho$sigma_share,
      case$rho$sigma_margin, case$rho$rho, G = case$G
    )
    rho_note <- "augmented share plus every second margin; dense oracle is original block"
  } else {
    input <- case
    structured_fun <- function() structured_gaussian_log_density(
      input$z, input$v, input$group, input$x, case$tau, G = case$G
    )
    n_dense <- case$N
    dense_fun <- function() dense_gaussian_log_density(
      input$z, input$v, input$group, input$x, case$tau, G = case$G
    )
    rho_note <- ""
  }
  equality_error <- NA_real_
  dense_status <- "ok"
  if (case$N <= dense_limit) {
    ## One equality calculation is kept separate from timing so both method
    ## rows carry the same directly comparable error.
    structured_once <- structured_fun()
    dense_once <- dense_fun()
    equality_error <- max_component_error(structured_once, dense_once)
  } else {
    dense_status <- paste0(
      "skipped_dense_covariance_estimate_", format_bytes(dense_bytes(n_dense)),
      "_exceeds_safe_limit"
    )
  }
  structured_time <- median_elapsed(
    structured_fun, structured_repetitions, batch = TRUE
  )
  ## Small dense matrices can also fall below the elapsed-time clock tick.
  ## Batch only N=100, where the dense allocation is at most a few hundred
  ## KiB; retain one expensive dense call per replicate at N=1000.
  dense_batch <- case$N == 100L
  dense_time <- if (case$N <= dense_limit) {
    median_elapsed(
      dense_fun, dense_repetitions, batch = dense_batch,
      target_seconds = 0.01
    )
  } else {
    list(median_seconds = NA_real_, inner_iterations = 0L)
  }
  bytes <- dense_bytes(n_dense)
  common <- data.frame(
    N = case$N,
    layout = case$layout,
    rho = if (rho_flag) "TRUE" else "FALSE",
    repetitions = structured_repetitions,
    median_seconds = structured_time$median_seconds,
    inner_iterations = structured_time$inner_iterations,
    dense_byte_estimate = bytes,
    status = "ok",
    equality_error = equality_error,
    notes = paste(case$layout_note, rho_note),
    method = "structured"
  )
  dense_row <- data.frame(
    N = case$N,
    layout = case$layout,
    rho = if (rho_flag) "TRUE" else "FALSE",
    repetitions = if (case$N <= dense_limit) dense_repetitions else 0L,
    median_seconds = dense_time$median_seconds,
    inner_iterations = dense_time$inner_iterations,
    dense_byte_estimate = bytes,
    status = dense_status,
    equality_error = equality_error,
    notes = paste(case$layout_note, rho_note),
    method = "dense"
  )
  rbind(common, dense_row)
}

run_benchmark <- function(verbose = TRUE) {
  sizes <- c(100L, 1000L, 10000L)
  layouts <- c("one_product", "multiproduct_panel")
  rows <- list()
  row_id <- 0L
  for (layout_id in seq_along(layouts)) {
    layout <- layouts[[layout_id]]
    for (size_id in seq_along(sizes)) {
      N <- sizes[[size_id]]
      case <- make_benchmark_case(N, layout, seed = 9000L +
                                    100L * layout_id + size_id)
      structured_repetitions <- if (N >= 10000L) 5L else 7L
      dense_repetitions <- 3L
      for (rho_flag in c(FALSE, TRUE)) {
        row_id <- row_id + 1L
        rows[[row_id]] <- benchmark_one_case(
          case, rho_flag, structured_repetitions, dense_repetitions
        )
      }
    }
  }
  benchmark <- do.call(rbind, rows)
  benchmark <- benchmark[, c(
    "N", "layout", "method", "rho", "repetitions", "inner_iterations",
    "median_seconds", "dense_byte_estimate", "status", "equality_error", "notes"
  )]
  if (verbose) {
    cat("\nStructured/dense timing benchmark (seconds are median elapsed time):\n")
    print(benchmark, row.names = FALSE)
    cat("\nDense covariance byte estimates are covariance storage only; actual\n")
    cat("dense evaluations allocate additional working matrices.\n")
  }
  benchmark
}

write_outputs <- function(validation, benchmark,
                          validation_path = file.path("doc", "collapsed_gaussian_structured_validation.csv"),
                          benchmark_path = file.path("doc", "collapsed_gaussian_structured_benchmark.csv")) {
  write.csv(validation, validation_path, row.names = FALSE, quote = TRUE)
  write.csv(benchmark, benchmark_path, row.names = FALSE, quote = TRUE)
  invisible(list(validation = validation_path, benchmark = benchmark_path))
}

main <- function() {
  validation <- run_structured_validation(verbose = TRUE)
  benchmark <- run_benchmark(verbose = TRUE)
  output_args <- commandArgs(trailingOnly = TRUE)
  output_args <- output_args[grepl("^--(validation|benchmark)=", output_args)]
  validation_path <- file.path("doc", "collapsed_gaussian_structured_validation.csv")
  benchmark_path <- file.path("doc", "collapsed_gaussian_structured_benchmark.csv")
  for (arg in output_args) {
    if (startsWith(arg, "--validation=")) validation_path <- sub("^--validation=", "", arg)
    if (startsWith(arg, "--benchmark=")) benchmark_path <- sub("^--benchmark=", "", arg)
  }
  write_outputs(validation, benchmark, validation_path, benchmark_path)
  cat("\nWrote ", validation_path, "\n", sep = "")
  cat("Wrote ", benchmark_path, "\n", sep = "")
  invisible(list(validation = validation, benchmark = benchmark))
}

if (identical(environment(), globalenv())) main()
