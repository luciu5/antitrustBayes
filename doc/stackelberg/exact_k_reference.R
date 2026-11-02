# Independent R oracle for exact ordered-k deterministic density checks.
# This file is intentionally Stan-free and does not perform posterior sampling.

.exact_k_match <- function(x, choices, name) {
  ans <- match.arg(x, choices)
  if (length(ans) != 1L) stop(name, " must be scalar")
  ans
}

.exact_k_firm_table <- function(firm, size, rank) {
  firm <- as.character(firm)
  stopifnot(length(firm) == length(size), length(size) == length(rank))
  idx <- !duplicated(firm)
  tab <- data.frame(firm = firm[idx], size = size[idx], rank = rank[idx])
  for (ff in tab$firm) {
    ii <- firm == ff
    stopifnot(length(unique(size[ii])) == 1L, length(unique(rank[ii])) == 1L)
  }
  tab
}

exact_k_logit_inverse_lerner <- function(price, firm, firm_share_in, rank,
                                          s0, alpha, k,
                                          structure = c("fringe", "stackelberg"),
                                          conduct = c("bertrand", "cournot")) {
  structure <- .exact_k_match(structure, c("fringe", "stackelberg"), "structure")
  conduct <- .exact_k_match(conduct, c("bertrand", "cournot"), "conduct")
  stopifnot(length(price) == length(firm), length(firm) == length(firm_share_in),
            length(rank) == length(firm), is.finite(s0), s0 > 0, s0 < 1,
            is.finite(alpha), alpha > 0, length(k) == 1L, k >= 0)
  ft <- .exact_k_firm_table(firm, firm_share_in, rank)
  S <- ft$size * (1 - s0)
  is_core <- ft$rank <= k
  h <- numeric(nrow(ft))
  if (structure == "fringe") {
    if (conduct == "bertrand") h <- ifelse(is_core, 1 / (1 - S), 1)
    else h <- ifelse(is_core, 1 + S / s0, 1)
  } else {
    stopifnot(k >= 1)
    if (conduct == "bertrand") {
      B <- sum(S[!is_core]^2 / (1 - S[!is_core] + S[!is_core]^2))
      phi <- 1 / (1 - B)
      h <- ifelse(is_core, 1 / (1 - phi * S), 1 / (1 - S))
    } else {
      SF <- sum(S[!is_core])
      h <- ifelse(is_core, 1 + S / (s0 + SF), 1 + S / s0)
    }
  }
  price * alpha / h[match(as.character(firm), ft$firm)]
}

exact_k_ces_inverse_lerner <- function(firm, firm_rev_share_in, rank,
                                        outside_quantity_share,
                                        outside_price, inside_price_index,
                                        gamma, k,
                                        structure = c("fringe", "stackelberg"),
                                        conduct = c("bertrand", "cournot")) {
  structure <- .exact_k_match(structure, c("fringe", "stackelberg"), "structure")
  conduct <- .exact_k_match(conduct, c("bertrand", "cournot"), "conduct")
  stopifnot(length(firm) == length(firm_rev_share_in), length(rank) == length(firm),
            outside_quantity_share > 0, outside_quantity_share < 1,
            outside_price > 0, inside_price_index > 0, gamma > 1,
            length(k) == 1L, k >= 0)
  ft <- .exact_k_firm_table(firm, firm_rev_share_in, rank)
  R0 <- plogis(qlogis(outside_quantity_share) + log(outside_price) -
                 log(inside_price_index))
  R <- ft$size * (1 - R0)
  core <- ft$rank <= k
  h <- gamma - 1
  H <- 1 + h * R0
  b <- h / (gamma * H)
  m <- numeric(nrow(ft))
  if (structure == "fringe") {
    if (conduct == "bertrand") m <- ifelse(core, 1 / (gamma - h * R), 1 / gamma)
    else m <- ifelse(core, 1 / gamma + b * R, 1 / gamma)
  } else {
    stopifnot(k >= 1)
    if (conduct == "bertrand") {
      t <- gamma - h * R
      a <- h^2 * R / (t * (t - 1) + h^2 * R)
      phi <- 1 / (1 - sum(R[!core] * a[!core]))
      m <- ifelse(core, 1 / (gamma - h * phi * R), 1 / (gamma - h * R))
    } else {
      a <- h * (H - 2 * R + gamma * h * R0 * R / H) /
        (H * (H + (gamma - 2) * R))
      C <- sum(R[!core] * a[!core])
      m <- ifelse(core, 1 / gamma + b * R / (1 + C), 1 / gamma + b * R)
    }
  }
  1 / m[match(as.character(firm), ft$firm)]
}

exact_k_market_mixture <- function(log_lik_by_k, market, k_support, k_prior) {
  log_lik_by_k <- as.matrix(log_lik_by_k)
  market <- as.character(market)
  stopifnot(nrow(log_lik_by_k) == length(market),
            ncol(log_lik_by_k) == length(k_support),
            length(k_prior) == length(k_support),
            all(is.finite(k_prior)), all(k_prior >= 0),
            abs(sum(k_prior) - 1) <= 1e-8)
  log_prior <- ifelse(k_prior == 0, -Inf, log(k_prior))
  markets <- unique(market)
  scores <- posterior <- matrix(NA_real_, length(markets), length(k_support),
                                 dimnames = list(markets, as.character(k_support)))
  log_mix <- setNames(numeric(length(markets)), markets)
  for (mm in seq_along(markets)) {
    ii <- market == markets[mm]
    scores[mm, ] <- log_prior + colSums(log_lik_by_k[ii, , drop = FALSE])
    mx <- max(scores[mm, ])
    log_mix[mm] <- mx + log(sum(exp(scores[mm, ] - mx)))
    posterior[mm, ] <- exp(scores[mm, ] - log_mix[mm])
  }
  list(score = scores, posterior = posterior, log_mix = log_mix)
}
