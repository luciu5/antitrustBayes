## Normalized Gaussian densities on within-cell orthonormal share contrasts.
## Cell levels are computational nuisance coordinates, never model parameters.

.relative_share_labels <- function(x, n, name) {
  if (!is.atomic(x) || !is.null(dim(x)) || length(x) != n || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x)))) {
    stop("'", name, "' must contain ", n, " finite, non-missing labels.",
         call. = FALSE)
  }
  x <- as.character(x)
  if (any(!nzchar(x))) stop("'", name, "' contains empty labels.", call. = FALSE)
  match(x, unique(x))
}

.relative_share_basis <- function(cell) {
  n <- length(cell)
  if (!n) stop("'cell' must contain at least one observation.", call. = FALSE)
  ids <- .relative_share_labels(cell, n, "cell")
  indices <- split(seq_len(n), factor(ids, levels = seq_len(max(ids))))
  out <- matrix(0, n, n - length(indices))
  column <- 0L
  for (ii in indices) if (length(ii) > 1L) {
    for (k in seq_len(length(ii) - 1L)) {
      out[ii[seq_len(k)], column + k] <- 1 / sqrt(k * (k + 1))
      out[ii[k + 1L], column + k] <- -k / sqrt(k * (k + 1))
    }
    column <- column + length(ii) - 1L
  }
  out
}

.relative_share_harmonic <- function(x, y) {
  smaller <- pmin(x, y)
  smaller / (1 + smaller / pmax(x, y))
}

.relative_share_center_quality <- function(x, blocks) {
  for (ii in blocks) x[ii, ] <- sweep(x[ii, , drop = FALSE], 2L,
                                     colMeans(x[ii, , drop = FALSE]), "-")
  x
}

.relative_share_factor <- function(variance, quality_index, quality_block,
                                   cell, tau) {
  n <- length(quality_index)
  g <- length(quality_block)
  if (!n || !g || !is.numeric(quality_index) ||
      any(!is.finite(quality_index)) || any(quality_index != floor(quality_index)) ||
      any(quality_index < 1L | quality_index > g)) {
    stop("'quality_index' must index the nonempty 'quality_block' vector.",
         call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0 ||
      tau > sqrt(.Machine$double.xmax)) {
    stop("'tau' must be a finite nonnegative, numerically usable scalar.",
         call. = FALSE)
  }
  if (!is.numeric(variance) || !length(variance) %in% c(1L, n) ||
      any(!is.finite(variance)) || any(variance <= 0)) {
    stop("'variance' must contain one or N finite positive variances.", call. = FALSE)
  }
  variance <- rep(variance, length.out = n)
  block_id <- .relative_share_labels(quality_block, g, "quality_block")
  cell_id <- .relative_share_labels(cell, n, "cell")
  original_cells <- split(seq_len(n), factor(cell_id, levels = seq_len(max(cell_id))))
  if (any(vapply(original_cells, function(ii)
    length(unique(block_id[quality_index[ii]])) != 1L, logical(1)))) {
    stop("Each contrast cell must lie within one quality-centering block.",
         call. = FALSE)
  }
  blocks <- split(seq_len(g), factor(block_id, levels = seq_len(max(block_id))))
  ## A singleton has no observed share contrast. Keep its quality prior and
  ## original output row, but exclude it from the numerical factorization.
  active_cells <- original_cells[lengths(original_cells) > 1L]
  active <- sort(as.integer(unlist(active_cells, use.names = FALSE)))
  out <- list(n = n, g = g, tau = tau, tau2 = tau^2, blocks = blocks,
              active = active, dimension = n - length(original_cells))
  if (!length(active)) return(out)
  ids <- match(cell_id[active], unique(cell_id[active]))
  out$cells <- split(seq_along(active), factor(ids, levels = seq_len(max(ids))))
  out$cell_id <- ids
  out$group <- as.integer(quality_index[active])
  out$groups <- split(seq_along(active), factor(out$group, levels = seq_len(g)))
  out$variance <- variance[active]
  out$weight <- 1 / out$variance
  if (any(!is.finite(out$weight))) {
    stop("Share variances are too small for finite Gaussian precision.", call. = FALSE)
  }
  out$total <- unname(vapply(out$groups, function(ii) sum(out$weight[ii]), numeric(1)))
  if (any(!is.finite(out$total))) {
    stop("Accumulated share precision is not finite.", call. = FALSE)
  }
  out$h <- if (out$tau2 == 0) rep(0, g) else
    .relative_share_harmonic(rep(out$tau2, g), 1 / out$total)
  out$a <- if (out$tau2 == 0) out$total else
    .relative_share_harmonic(out$total, rep(1 / out$tau2, g))
  nc <- length(out$cells)
  ridge <- numeric(nc)
  edge <- matrix(0, nc, nc)
  out$group_cell <- vector("list", g)
  for (gg in seq_len(g)) {
    ii <- out$groups[[gg]]
    if (!length(ii)) next
    parts <- split(out$weight[ii], ids[ii])
    jj <- as.integer(names(parts))
    b <- vapply(parts, sum, numeric(1))
    p <- b / out$total[gg]
    out$group_cell[[gg]] <- list(cell = jj, fraction = unname(p))
    ridge[jj] <- ridge[jj] + out$a[gg] * p
    if (length(jj) > 1L && out$tau2 > 0) {
      ## h*b_c*b_d is evaluated as (t*h)*b_c*(b_d/t), avoiding overflow.
      contribution <- (out$total[gg] * out$h[gg]) * tcrossprod(b, p)
      edge[jj, jj] <- edge[jj, jj] + contribution
    }
  }
  diag(edge) <- 0
  edge <- (edge + t(edge)) / 2
  s <- -edge
  diag(s) <- ridge + rowSums(edge)
  if (any(!is.finite(s)) || any(!is.finite(ridge)) || any(ridge <= 0)) {
    stop("Cell-level Gaussian precision is not numerically usable.", call. = FALSE)
  }
  r <- tryCatch(chol(s), error = function(e) NULL)
  use_qr <- is.null(r) || rcond(r, triangular = TRUE) < .Machine$double.eps^.25
  out$pivot <- seq_len(nc)
  if (use_qr) {
    ## Form the positive square root instead of rounding a tiny ridge into a
    ## much larger Laplacian diagonal. This changes no variance or target.
    pairs <- which(upper.tri(edge) & edge > 0, arr.ind = TRUE)
    f <- matrix(0, nc + nrow(pairs), nc)
    f[cbind(seq_len(nc), seq_len(nc))] <- sqrt(ridge)
    if (nrow(pairs)) for (j in seq_len(nrow(pairs))) {
      value <- sqrt(edge[pairs[j, 1L], pairs[j, 2L]])
      f[nc + j, pairs[j, ]] <- c(value, -value)
    }
    qr_f <- qr(f, LAPACK = TRUE)
    r <- qr.R(qr_f)[seq_len(nc), , drop = FALSE]
    out$pivot <- qr_f$pivot
  }
  if (any(!is.finite(r)) || any(diag(r) == 0)) {
    stop("Cell-level Gaussian factorization is singular.", call. = FALSE)
  }
  out$factor <- r
  out$factorization <- if (use_qr) "positive_qr" else "cholesky"
  out$logdet_s <- 2 * sum(log(abs(diag(r))))
  extra <- if (out$tau2 == 0) rep(0, g) else {
    x <- log(out$tau2) + log(out$total)
    pmax(x, 0) + log1p(exp(-abs(x)))
  }
  out$logdet_a <- sum(log(out$variance)) + sum(extra)
  out$logdet <- out$logdet_a + out$logdet_s - sum(log(lengths(out$cells)))
  out
}

.relative_share_solve <- function(system, x) {
  x <- as.matrix(x)
  value <- backsolve(system$factor,
    forwardsolve(t(system$factor), x[system$pivot, , drop = FALSE]))
  out <- value
  out[system$pivot, ] <- value
  out
}

.relative_share_a_inverse <- function(system, x) {
  out <- matrix(0, nrow(x), ncol(x))
  for (g in seq_len(system$g)) {
    ii <- system$groups[[g]]
    if (!length(ii)) next
    p <- system$weight[ii] / system$total[g]
    average <- as.numeric(crossprod(p, x[ii, , drop = FALSE]))
    out[ii, ] <- system$weight[ii] * sweep(x[ii, , drop = FALSE], 2L, average, "-") +
      (system$a[g] * p) * rep(average, each = length(ii))
  }
  out
}

.relative_share_profile <- function(system, x) {
  x <- as.matrix(x)
  for (ii in system$cells) x[ii, ] <- sweep(x[ii, , drop = FALSE], 2L,
                                           colMeans(x[ii, , drop = FALSE]), "-")
  ax <- .relative_share_a_inverse(system, x)
  rhs <- do.call(rbind, lapply(system$cells, function(ii)
    colSums(ax[ii, , drop = FALSE])))
  coefficient <- .relative_share_solve(system, rhs)
  residual <- x - coefficient[system$cell_id, , drop = FALSE]
  list(residual = residual, precision = .relative_share_a_inverse(system, residual))
}

.relative_share_quality_update <- function(system, residual) {
  out <- matrix(0, system$g, ncol(residual))
  for (g in seq_len(system$g)) {
    ii <- system$groups[[g]]
    if (!length(ii)) next
    average <- as.numeric(crossprod(system$weight[ii] / system$total[g],
                                    residual[ii, , drop = FALSE]))
    out[g, ] <- (system$total[g] * system$h[g]) * average
  }
  out
}

.relative_share_input <- function(x, n, name = "residual") {
  if (!is.numeric(x) || any(!is.finite(x)) ||
      (is.null(dim(x)) && length(x) != n) ||
      (!is.null(dim(x)) && (length(dim(x)) != 2L || nrow(x) != n))) {
    stop("'", name, "' must be a finite length-N vector or N-row matrix.",
         call. = FALSE)
  }
  as.matrix(x)
}

.relative_share_precision <- function(x, variance, quality_index, quality_block,
                                      cell, tau) {
  system <- .relative_share_factor(variance, quality_index, quality_block, cell, tau)
  matrix_input <- !is.null(dim(x))
  x <- .relative_share_input(x, system$n, "x")
  out <- matrix(0, system$n, ncol(x))
  if (system$dimension) out[system$active, ] <-
    .relative_share_profile(system, x[system$active, , drop = FALSE])$precision
  if (matrix_input) out else as.numeric(out)
}

.relative_share_gaussian <- function(residual, variance, quality_index,
                                     quality_block, cell, tau,
                                     return_moments = FALSE) {
  system <- .relative_share_factor(variance, quality_index, quality_block, cell, tau)
  x <- .relative_share_input(residual, system$n)
  if (ncol(x) != 1L) stop("'residual' must contain one observation vector.", call. = FALSE)
  if (!is.logical(return_moments) || length(return_moments) != 1L || is.na(return_moments))
    stop("'return_moments' must be TRUE or FALSE.", call. = FALSE)
  out <- list(log_density = 0, dimension = as.integer(system$dimension))
  if (system$dimension) {
    profile <- .relative_share_profile(system, x[system$active, , drop = FALSE])
    r <- as.numeric(profile$residual)
    quadratic <- 0
    for (g in seq_len(system$g)) {
      ii <- system$groups[[g]]
      if (!length(ii)) next
      average <- sum((system$weight[ii] / system$total[g]) * r[ii])
      quadratic <- quadratic + sum(system$weight[ii] * (r[ii] - average)^2) +
        system$a[g] * average^2
    }
    out$log_density <- as.numeric(-.5 *
      (system$dimension * log(2 * pi) + system$logdet + quadratic))
  }
  if (return_moments) {
    mean <- numeric(system$g)
    variance_q <- numeric(system$g)
    if (system$dimension) {
      mean <- as.numeric(.relative_share_quality_update(system, profile$residual))
      f <- matrix(0, system$g, length(system$cells))
      for (g in seq_len(system$g)) if (length(system$groups[[g]])) {
        b <- system$group_cell[[g]]
        f[g, b$cell] <- (system$total[g] * system$h[g]) * b$fraction
      }
      f <- .relative_share_center_quality(f, system$blocks)
      white <- forwardsolve(t(system$factor), t(f)[system$pivot, , drop = FALSE])
      variance_q <- colSums(white^2)
      for (ii in system$blocks) if (length(ii) > 1L) {
        variance_q[ii] <- variance_q[ii] + system$h[ii] * (1 - 2 / length(ii)) +
          sum(system$h[ii]) / length(ii)^2
      } else variance_q[ii] <- 0
    } else {
      for (ii in system$blocks) variance_q[ii] <- system$tau2 * (1 - 1 / length(ii))
    }
    out$conditional_mean <- as.numeric(.relative_share_center_quality(matrix(mean, ncol = 1L), system$blocks))
    out$conditional_variance <- variance_q
  }
  out
}

.relative_share_quality_draw <- function(residual, variance, quality_index,
                                        quality_block, cell, tau,
                                        n_draws = 1L, seed = NULL) {
  system <- .relative_share_factor(variance, quality_index, quality_block, cell, tau)
  x <- .relative_share_input(residual, system$n)
  if (ncol(x) != 1L) stop("'residual' must contain one observation vector.", call. = FALSE)
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1L || n_draws != floor(n_draws) || n_draws > .Machine$integer.max)
    stop("'n_draws' must be a positive integer.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max))
    stop("'seed' must be NULL or a nonnegative integer.", call. = FALSE)
  calculate <- function() {
    prior <- .relative_share_center_quality(
      matrix(stats::rnorm(system$g * n_draws), system$g, n_draws) * tau,
      system$blocks)
    if (system$dimension && system$tau2 > 0) {
      simulated <- prior[system$group, , drop = FALSE] +
        sqrt(system$variance) * matrix(stats::rnorm(length(system$active) * n_draws),
                                      length(system$active), n_draws)
      difference <- as.numeric(x[system$active, 1L]) - simulated
      profile <- .relative_share_profile(system, difference)
      update <- .relative_share_quality_update(system, profile$residual)
      prior <- prior + .relative_share_center_quality(update, system$blocks)
    }
    t(.relative_share_center_quality(prior, system$blocks))
  }
  if (is.null(seed)) return(calculate())
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  calculate()
}
