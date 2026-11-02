## Lightweight independent FOC check for the collapsed-Gaussian benchmark.
## Run from the repository root with: Rscript doc/collapsed_gaussian_oracle_validation.R
## It performs no compilation, fitting, sampling, or package installation.

benchmark_path <- file.path("doc", "collapsed_gaussian_benchmark.R")
lines <- readLines(benchmark_path, warn = FALSE)
entry <- which(grepl("if (sys.nframe() == 0L)", lines, fixed = TRUE))[1L]
if (!is.finite(entry) || entry <= 1L) stop("benchmark CLI guard not found")
env <- new.env(parent = globalenv())
eval(parse(text = lines[seq_len(entry - 1L)]), envir = env)

out <- list()
push <- function(family, structure, supply_model, n_draws, max_abs_error) {
  out[[length(out) + 1L]] <<- data.frame(
    family = family, structure = structure, supply_model = supply_model,
    n_draws = as.integer(n_draws), max_abs_error = max_abs_error,
    stringsAsFactors = FALSE
  )
}

set.seed(11)
for (structure in c("simultaneous", "fringe", "stackelberg")) {
  for (supply_model in c("bertrand", "cournot")) {
    max_error <- 0
    for (i in seq_len(1000L)) {
      n <- sample(3:10, 1L)
      q <- runif(n, .01, .9); q <- q / sum(q) * runif(1, .3, .9)
      s0 <- runif(1, .01, .5)
      rank <- order(order(-q))
      k <- if (structure == "simultaneous") NULL else sample(1:(n - 1L), 1L)
      core <- if (structure == "simultaneous") rep(TRUE, n) else rank <= k
      S <- q * (1 - s0)
      response <- if (structure == "stackelberg" && supply_model == "bertrand") {
        1 / (1 - sum(S[!core]^2 / (1 - S[!core] + S[!core]^2)))
      } else if (structure == "stackelberg") sum(S[!core]) else 0
      h <- if (structure == "fringe") {
        if (supply_model == "bertrand") ifelse(core, 1 / (1 - S), 1) else
          ifelse(core, 1 + S / s0, 1)
      } else if (structure == "stackelberg") {
        if (supply_model == "bertrand") ifelse(core, 1 / (1 - response * S), 1 / (1 - S)) else
          ifelse(core, 1 + S / (s0 + response), 1 + S / s0)
      } else if (supply_model == "bertrand") 1 / (1 - S) else 1 + S / s0
      production <- 1 / h
      oracle <- env$cga_role_conduct_logit(q, s0, supply_model, structure, rank, k)
      max_error <- max(max_error, max(abs(production - oracle)))
    }
    push("logit", structure, supply_model, 1000L, max_error)
  }
}

for (structure in c("simultaneous", "fringe", "stackelberg")) {
  for (supply_model in c("bertrand", "cournot")) {
    max_error <- 0
    for (i in seq_len(1000L)) {
      n <- sample(3:10, 1L)
      p <- runif(n, .5, 2); q <- runif(n, .01, .9)
      q <- q / sum(q) * runif(1, .3, .9)
      rev <- p * q / sum(p * q)
      s0 <- runif(1, .01, .5); p0 <- runif(1, .5, 2)
      avg <- sum(p * q)
      R0 <- s0 * p0 / (s0 * p0 + (1 - s0) * avg)
      gamma <- runif(1, 1.2, 8); rank <- order(order(-rev))
      k <- if (structure == "simultaneous") NULL else sample(1:(n - 1L), 1L)
      core <- if (structure == "simultaneous") rep(TRUE, n) else rank <= k
      R <- rev * (1 - R0); h <- gamma - 1; H <- 1 + h * R0
      b <- h / (gamma * H)
      if (structure == "fringe") {
        production <- ifelse(core,
          if (supply_model == "bertrand") gamma - h * R else 1 / (1 / gamma + b * R),
          gamma)
      } else if (structure == "stackelberg") {
        if (supply_model == "bertrand") {
          a <- h^2 * R / ((gamma - h * R) * (gamma - h * R - 1) + h^2 * R)
          response <- 1 / (1 - sum(R[!core] * a[!core]))
          production <- ifelse(core, gamma - h * response * R, gamma - h * R)
        } else {
          a <- h * (H - 2 * R + gamma * h * R0 * R / H) /
            (H * (H + (gamma - 2) * R))
          response <- sum(R[!core] * a[!core])
          production <- 1 / ifelse(core, 1 / gamma + b * R / (1 + response),
                                   1 / gamma + b * R)
        }
      } else if (supply_model == "bertrand") production <- gamma - h * R else
        production <- 1 / (1 / gamma + h * R / (gamma * (1 + h * R0)))
      oracle <- env$cga_role_inverse_lerner_ces(
        rev, R0, gamma, supply_model, structure, rank, k
      )
      max_error <- max(max_error, max(abs(production - oracle)))
    }
    push("ces", structure, supply_model, 1000L, max_error)
  }
}

results <- do.call(rbind, out)
csv_path <- file.path("doc", "collapsed_gaussian_oracle_validation.csv")
utils::write.csv(results, csv_path, row.names = FALSE)
print(results, row.names = FALSE)
