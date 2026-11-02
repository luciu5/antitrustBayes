# A linear-algebra reference benchmark; no package implementation is changed.
# Run with single-threaded BLAS for a controlled comparison:
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript ownership_rhs_reference.R out.csv
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
rows <- list()
for (n in c(100L, 300L, 1000L)) {
  for (products_per_firm in c(1L, 3L)) {
    firm <- ceiling(seq_len(n) / products_per_firm)
    owner <- outer(firm, firm, "==") * 1
    share <- rep(0.7 / n, n)
    price <- seq(2, 3, length.out = n)
    elasticity <- matrix(share * price, n, n, byrow = TRUE)
    diag(elasticity) <- diag(elasticity) - price
    a <- t(elasticity) * owner
    revenue <- share * price
    rhs <- revenue * diag(owner)
    inverse <- function() as.vector(solve(a) %*% rhs)
    direct <- function() as.vector(solve(a, rhs))
    blocks <- split(seq_len(n), firm)
    block_reference <- function() {
      ans <- numeric(n)
      for (index in blocks) ans[index] <- solve(a[index, index, drop = FALSE], rhs[index])
      ans
    }
    # Warm both paths, then alternate their order in paired repetitions.
    invisible(inverse()); invisible(direct()); invisible(block_reference())
    for (repetition in seq_len(3L)) {
      order <- if (repetition %% 2L) c("inverse", "direct", "block") else c("block", "direct", "inverse")
      times <- numeric(3L); names(times) <- c("inverse", "direct", "block")
      values <- list()
      for (method in order) {
        fun <- switch(method, inverse = inverse, direct = direct, block = block_reference)
        times[method] <- system.time(values[[method]] <- fun())[["elapsed"]]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        n, products_per_firm, repetition,
        inverse_seconds = unname(times["inverse"]),
        direct_seconds = unname(times["direct"]),
        block_seconds = unname(times["block"]),
        max_solution_difference = max(abs(values$inverse - values$direct)),
        max_margin_difference = max(abs((values$inverse - values$direct) / revenue)),
        direct_residual = max(abs(a %*% values$direct - rhs)),
        max_block_margin_difference = max(abs((values$inverse - values$block) / revenue)),
        stage = "reference_only_no_equilibrium_or_fallback_test"
      )
    }
  }
}
write.csv(do.call(rbind, rows), args[[1L]], row.names = FALSE)
