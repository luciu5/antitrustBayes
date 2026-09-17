## Independent numerical and timing check for the sampled-quality Helmert
## contrast in doc/performance_audit/sampled_quality_contrast.md.
## This is a reference calculation only; it does not modify Stan or R source.

out_dir <- "/tmp/performance-stan-profile"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dense_Q <- function(G) {
  G <- as.integer(G)
  if (G == 1L) return(matrix(numeric(0), nrow = 1L, ncol = 0L))
  Q <- matrix(0, nrow = G, ncol = G - 1L)
  for (k in seq_len(G - 1L)) {
    scale <- 1 / sqrt(as.numeric(k) * (as.numeric(k) + 1))
    Q[seq_len(k), k] <- scale
    Q[k + 1L, k] <- -k * scale
  }
  Q
}

recurrence <- function(beta) {
  G <- length(beta) + 1L
  if (G == 1L) return(0)
  k <- seq_len(G - 1L)
  v <- beta / sqrt(as.numeric(k) * (as.numeric(k) + 1))
  out <- numeric(G)
  a <- 0
  for (j in G:2L) {
    out[j] <- a - (j - 1L) * v[j - 1L]
    a <- a + v[j - 1L]
  }
  out[1L] <- a
  out
}

max_abs <- function(x) if (length(x)) max(abs(x)) else 0
max_rel <- function(x, ref) {
  if (!length(x)) return(0)
  max(abs(x) / (1 + abs(ref)))
}

Gs <- c(1L, 3L, 10L, 100L, 1000L)
checks <- list()
timings <- list()
set.seed(48201L)
for (G in Gs) {
  beta <- if (G == 1L) numeric(0) else sin(seq_len(G - 1L) / 7)
  direction <- if (G == 1L) numeric(0) else cos(seq_len(G - 1L) / 11)
  Q <- dense_Q(G)
  dense <- if (G == 1L) 0 else as.numeric(Q %*% beta)
  fast <- recurrence(beta)
  dense_direction <- if (G == 1L) 0 else as.numeric(Q %*% direction)
  fast_direction <- recurrence(direction)
  h <- 1e-6
  dense_fd <- if (G == 1L) 0 else as.numeric(
    (Q %*% (beta + h * direction) - Q %*% (beta - h * direction)) / (2 * h)
  )
  fast_fd <- (recurrence(beta + h * direction) -
                recurrence(beta - h * direction)) / (2 * h)
  checks[[length(checks) + 1L]] <- data.frame(
    G = G,
    output_max_abs = max_abs(dense - fast),
    output_max_rel = max_rel(dense - fast, dense),
    dense_centered_sum = sum(dense), fast_centered_sum = sum(fast),
    centered_sum_abs = abs(sum(fast)),
    directional_max_abs = max_abs(dense_direction - fast_direction),
    directional_max_rel = max_rel(dense_direction - fast_direction, dense_direction),
    directional_centered_sum = sum(fast_direction),
    finite_difference_dense_error = max_abs(dense_fd - dense_direction),
    finite_difference_recurrence_error = max_abs(fast_fd - fast_direction),
    stringsAsFactors = FALSE
  )

  reps <- if (G >= 1000L) 10L else 100L
  dense_time <- system.time(replicate(reps, dense_Q(G)))[["elapsed"]] / reps
  recurrence_time <- system.time(replicate(reps, recurrence(beta)))[["elapsed"]] / reps
  timings[[length(timings) + 1L]] <- data.frame(
    G = G, repetitions = reps, dense_seconds = dense_time,
    recurrence_seconds = recurrence_time,
    dense_over_recurrence = dense_time / recurrence_time,
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, checks), file.path(out_dir, "helmert_check_v2.csv"), row.names = FALSE)
write.csv(do.call(rbind, timings), file.path(out_dir, "helmert_timing_v2.csv"), row.names = FALSE)

## The G=10,000 dense matrix is deliberately not allocated here.
storage <- data.frame(
  G = 10000L, status = "skipped_storage_cap",
  dense_Q_bytes = 8 * 10000 * 9999,
  dense_Q_mib = (8 * 10000 * 9999) / 2^20,
  note = "Reference storage estimate; no dense native allocation.",
  stringsAsFactors = FALSE
)
write.csv(storage, file.path(out_dir, "helmert_storage_estimate_v2.csv"), row.names = FALSE)
cat("helmert check complete\n")
