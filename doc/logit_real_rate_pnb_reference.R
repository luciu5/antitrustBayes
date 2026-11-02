#!/usr/bin/env Rscript

## Independent latent-quality PNB real-price oracle.
## Reads saved fits/plans/scenarios/arrays only.  It does not call a production
## price solver, native target, bridge, Stan sampling, HMC, or nonlinear solver.

suppressPackageStartupMessages(library(antitrustBayes))
args <- commandArgs(trailingOnly = TRUE)
mq <- normalizePath(if (length(args)) args[[1L]] else
  file.path("doc", "statistical_toolkit_runs", "20260916", "performance-audit",
            "recovery", "market-quality"),
  mustWork = TRUE)
out <- if (length(args) >= 2L) args[[2L]] else
  file.path(mq, paste0("rate-domain-oracle-", format(Sys.time(), "%Y%m%dT%H%M%S")))
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
archive <- file.path(mq, "pnb-archived-retry-01")
policy <- file.path(archive, "postprocess-run-01", "analysis", "policy")
stopf <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(ok, ...) if (!isTRUE(ok)) stopf(...)

logsumexp <- function(x) {
  assert(length(x) > 0L && all(is.finite(x)), "invalid logsumexp input")
  z <- max(x); z + log(sum(exp(x - z)))
}

## Solve log(x) + x = ell for x > 0, i.e. W(exp(ell)), on the real branch.
w_from_log <- function(ell) {
  assert(length(ell) == 1L && is.finite(ell), "invalid W log argument")
  lo <- min(-50, ell - 1); hi <- max(1, log(max(1, ell)) + 1)
  exp(uniroot(function(z) z + exp(z) - ell, c(lo, hi), tol = 1e-13)$root)
}

logit_shares <- function(delta, alpha, p_ref, p) {
  u <- delta + alpha * (p - p_ref)
  z <- max(c(0, u)); ei <- exp(u - z); e0 <- exp(-z)
  den <- e0 + sum(ei)
  list(inside = ei / den, outside = e0 / den, utility = u)
}

observed_h <- function(structure, owner, roles, odds) {
  if (identical(structure, "fringe"))
    return(ifelse(owner %in% roles, 1 + odds, 1))
  assert(identical(structure, "stackelberg"), "unexpected structure")
  reserve <- 1 + sum(odds[!owner %in% roles])
  ifelse(owner %in% roles, 1 + odds / reserve, 1 + odds)
}

## Exact CF/CS ownership-state oracle. `logA` is fixed across the merger.
solve_oracle <- function(structure, owner, strategic, logA, value, alpha,
                          p_ref, delta) {
  firms <- unique(owner)
  logAf <- setNames(vapply(firms, function(f)
    logsumexp(logA[owner == f]), numeric(1)), firms)
  x <- setNames(numeric(length(firms)), firms)
  h <- setNames(numeric(length(firms)), firms)
  if (identical(structure, "fringe")) {
    for (f in firms) {
      if (f %in% strategic) {
        x[[f]] <- w_from_log(logAf[[f]] - 1); h[[f]] <- 1 + x[[f]]
      } else {
        x[[f]] <- exp(logAf[[f]] - 1); h[[f]] <- 1
      }
    }
  } else {
    followers <- setdiff(firms, strategic)
    x[followers] <- vapply(followers, function(f)
      w_from_log(logAf[[f]] - 1), numeric(1))
    reserve <- 1 + sum(x[followers])
    leaders <- intersect(firms, strategic)
    x[leaders] <- vapply(leaders, function(f)
      reserve * w_from_log(logAf[[f]] - 1 - log(reserve)), numeric(1))
    h[followers] <- 1 + x[followers]
    h[leaders] <- 1 + x[leaders] / reserve
  }
  prices <- value - h[owner] / alpha
  sh <- logit_shares(delta, alpha, p_ref, prices)
  firm_odds <- setNames(vapply(firms, function(f)
    sum(sh$inside[owner == f]) / sh$outside, numeric(1)), firms)
  h_check <- setNames(numeric(length(firms)), firms)
  if (identical(structure, "fringe")) {
    h_check[strategic] <- 1 + firm_odds[strategic]
    h_check[setdiff(firms, strategic)] <- 1
  } else {
    followers <- setdiff(firms, strategic)
    reserve_check <- 1 + sum(firm_odds[followers])
    h_check[followers] <- 1 + firm_odds[followers]
    h_check[intersect(firms, strategic)] <-
      1 + firm_odds[intersect(firms, strategic)] / reserve_check
  }
  foc <- prices - (value - h_check[owner] / alpha)
  list(prices = as.numeric(prices), shares = as.numeric(sh$inside),
       outside = sh$outside, utility = sh$utility, foc = as.numeric(foc),
       h = h[owner], h_check = h_check[owner])
}

stage_class <- function(pre, post) {
  if (any(!is.finite(c(pre$prices, post$prices)))) return("nonfinite_oracle")
  a <- any(pre$prices <= 0); b <- any(post$prices <= 0)
  if (!a && !b) return("pre_positive_post_positive")
  if (a && b) return("pre_nonpositive_post_nonpositive")
  if (a) return("pre_nonpositive_post_positive")
  "pre_positive_post_nonpositive"
}

failure_stage <- function(native_failed, pre, post) {
  if (!native_failed) return("native_success")
  if (any(pre$prices <= 0)) return("pre")
  if (any(post$prices <= 0)) return("post")
  "native_failure_with_positive_oracle"
}

one_case <- function(case, anchor = c("latent_quality", "observed")) {
  anchor <- match.arg(anchor)
  fit <- readRDS(file.path(archive, paste0(case, "-fit.rds")))
  plan <- readRDS(file.path(policy, paste0(case, "-draw-plan.rds")))
  saved <- readRDS(file.path(policy, paste0(case, "-", anchor, "-draws.rds")))
  sim <- readRDS(file.path(policy, paste0(case, "-", anchor, "-sim.rds")))
  scenario <- readRDS(file.path(policy, paste0(case, "-policy-status.rds")))
  assert(identical(as.character(fit@market_side), "input"), case, ": not input")
  assert(identical(as.character(fit@demand_model), "logit"), case, ": not Logit")
  structure <- as.character(fit@metadata$structure)
  assert(structure %in% c("fringe", "stackelberg"), case, ": bad structure")
  assert(identical(as.character(scenario$acquirer), "298816") &&
           identical(as.character(scenario$target), "889018"),
         case, ": bad saved PNB scenario")
  assert(identical(as.integer(plan$draw_indices), as.integer(saved$draw_indices)),
         case, ": plan/array draw mismatch")
  product <- as.character(fit@data$group)
  assert(identical(as.character(saved$product_labels), product), case, ": product mismatch")

  ## This is the canonical rstan permuted draw axis.  The saved plan indices
  ## refer to this axis, not to raw chain-row order.
  draws <- rstan::extract(fit@stanfit,
    pars = c("a_market", "s0_prob", "share_mean"),
    permuted = TRUE, inc_warmup = FALSE)
  share_mean <- as.matrix(draws$share_mean)
  a_market <- as.matrix(draws$a_market)
  s0_prob <- as.matrix(draws$s0_prob)
  assert(nrow(share_mean) >= max(plan$draw_indices), case, ": short draw axis")
  assert(ncol(share_mean) == nrow(fit@data), case, ": share axis mismatch")

  p <- as.numeric(fit@data$prices)
  r <- as.numeric(fit@data$quantities) / sum(fit@data$quantities)
  ranked <- product[order(-r, product)]
  post_owner <- product
  post_owner[post_owner == as.character(scenario$target)] <- as.character(scenario$acquirer)
  q_col <- match(product, as.character(plan$quality_units))
  assert(!anyNA(q_col), case, ": quality/product map mismatch")
  native_failure <- as.character(saved$failures)
  rows <- vector("list", length(plan$draw_indices))
  oracle_rows <- vector("list", length(plan$draw_indices))

  for (pos in seq_along(plan$draw_indices)) {
    di <- as.integer(plan$draw_indices[[pos]])
    s0 <- as.numeric(s0_prob[di, 1L])
    alpha <- as.numeric(a_market[di, 1L]) / as.numeric(fit@price_sd)
    assert(is.finite(s0) && s0 > 0 && s0 < 1, case, " position ", pos, ": s0")
    assert(is.finite(alpha) && alpha > 0, case, " position ", pos, ": alpha")
    if (identical(anchor, "latent_quality")) {
      q <- as.numeric(plan$quality_draws[pos, q_col])
      score <- as.numeric(share_mean[di, ]) + q
      z <- score - max(score)
      latent_ref_share <- (1 - s0) * exp(z) / sum(exp(z))
    } else {
      ## Observed anchor is the same conditional share normalization used by
      ## build_conduct_object: observed conditional shares and retained s0.
      q <- rep(0, length(product))
      score <- log((1 - s0) * r / s0)
      latent_ref_share <- (1 - s0) * r
    }
    delta <- log(latent_ref_share / s0)
    k <- as.integer(plan$k_cache[[as.character(di)]][["PNB"]])
    assert(length(k) == 1L && k >= 1L && k <= length(ranked),
           case, " position ", pos, ": k")
    roles <- ranked[seq_len(k)]
    h_obs <- observed_h(structure, product, roles, (1 - s0) * r / s0)
    value <- p + h_obs / alpha
    ## At p_observed, delta is actual latent utility; observed FOC inversion
    ## contributes h_obs, so log A = delta + h_obs.
    logA <- delta + h_obs
    pre <- solve_oracle(structure, product, roles, logA, value, alpha, p, delta)
    post_roles <- unique(post_owner[product %in% roles])
    post <- solve_oracle(structure, post_owner, post_roles, logA, value, alpha, p, delta)
    oracle_rows[[pos]] <- data.frame(
      case = case, anchor = anchor, draw_position = pos,
      draw_index = di, product_position = seq_along(product), product_label = product,
      owner_pre = product, owner_post = post_owner,
      strategic_pre = product %in% roles, strategic_post = post_owner %in% post_roles,
      latent_reference_share = latent_ref_share, delta_reference = delta,
      value_primitive = value, logA = logA,
      pre_price = pre$prices, post_price = post$prices,
      pre_share = pre$shares, post_share = post$shares,
      pre_utility = pre$utility, post_utility = post$utility,
      pre_foc_residual = pre$foc, post_foc_residual = post$foc,
      stringsAsFactors = FALSE
    )

    success <- is.na(native_failure[[pos]])
    np <- as.numeric(sim@price_pre$PNB[pos, ])
    nq <- as.numeric(sim@price_post$PNB[pos, ])
    ns <- as.numeric(sim@shares_pre$PNB[pos, ])
    nt <- as.numeric(sim@shares_post$PNB[pos, ])
    npc <- as.numeric(sim@price_change$PNB[pos, ])
    if (success) assert(all(is.finite(c(np, nq, ns, nt, npc))),
                        case, " position ", pos, ": malformed success array")
    pc <- post$prices / pre$prices - 1
    i_pre <- which.min(pre$prices); i_post <- which.min(post$prices)
    pre_resid <- abs(pre$foc); post_resid <- abs(post$foc)
    rows[[pos]] <- data.frame(
      case = case, anchor = anchor, draw_position = pos, draw_index = di, k = k,
      native_success = success,
      native_failure = if (success) NA_character_ else native_failure[[pos]],
      oracle_stage = stage_class(pre, post),
      native_failure_stage = failure_stage(!success, pre, post),
      alpha = alpha, s0 = s0,
      min_pre_price = min(pre$prices), min_pre_product = product[[i_pre]],
      min_post_price = min(post$prices), min_post_product = post_owner[[i_post]],
      pre_nonpositive = any(pre$prices <= 0), post_nonpositive = any(post$prices <= 0),
      max_abs_pre_foc_residual = max(pre_resid), max_abs_post_foc_residual = max(post_resid),
      max_rel_pre_foc_residual = max(pre_resid / pmax(abs(pre$prices), 1)),
      max_rel_post_foc_residual = max(post_resid / pmax(abs(post$prices), 1)),
      max_abs_native_pre_price_diff = if (success) max(abs(pre$prices - np)) else NA_real_,
      max_abs_native_post_price_diff = if (success) max(abs(post$prices - nq)) else NA_real_,
      max_abs_native_pre_share_diff = if (success) max(abs(pre$shares - ns)) else NA_real_,
      max_abs_native_post_share_diff = if (success) max(abs(post$shares - nt)) else NA_real_,
      max_abs_native_price_change_diff = if (success) max(abs(pc - npc)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  list(ledger = do.call(rbind, rows), oracle_long = do.call(rbind, oracle_rows))
}

case_runs <- unlist(lapply(c("latent_quality", "observed"), function(anchor)
  lapply(c("CF", "CS"), one_case, anchor = anchor)), recursive = FALSE)
result <- do.call(rbind, lapply(case_runs, `[[`, "ledger"))
oracle_long <- do.call(rbind, lapply(case_runs, `[[`, "oracle_long"))
rownames(result) <- NULL
rownames(oracle_long) <- NULL
latent_result <- result[result$anchor == "latent_quality", , drop = FALSE]
observed_result <- result[result$anchor == "observed", , drop = FALSE]
latent_long <- oracle_long[oracle_long$anchor == "latent_quality", , drop = FALSE]
observed_long <- oracle_long[oracle_long$anchor == "observed", , drop = FALSE]
write.csv(result, file.path(out, "pnb-rate-classification.csv"), row.names = FALSE)
write.csv(oracle_long, file.path(out, "pnb-oracle-long.csv"), row.names = FALSE)
write.csv(latent_result, file.path(out, "latent-rate-classification.csv"), row.names = FALSE)
write.csv(observed_result, file.path(out, "observed-rate-classification.csv"), row.names = FALSE)
write.csv(latent_long, file.path(out, "latent-oracle-long.csv"), row.names = FALSE)
write.csv(observed_long, file.path(out, "observed-oracle-long.csv"), row.names = FALSE)
saveRDS(list(ledger = result, oracle_long = oracle_long,
             latent = list(ledger = latent_result, oracle_long = latent_long),
             observed = list(ledger = observed_result, oracle_long = observed_long)),
        file.path(out, "pnb-oracle-matrices.rds"))
saveRDS(list(ledger = latent_result, oracle_long = latent_long),
        file.path(out, "latent-oracle-matrices.rds"))
write.csv(latent_result[!latent_result$native_success, , drop = FALSE],
          file.path(out, "latent-failures.csv"), row.names = FALSE)
write.csv(observed_result[!observed_result$native_success, , drop = FALSE],
          file.path(out, "observed-failures.csv"), row.names = FALSE)

stage_summary <- as.data.frame(table(result$anchor, result$case, result$oracle_stage))
names(stage_summary) <- c("anchor", "case", "oracle_stage", "n")
native_summary <- as.data.frame(table(result$anchor, result$case, result$native_success))
names(native_summary) <- c("anchor", "case", "native_success", "n")
fail_summary <- as.data.frame(table(
  result$anchor[!result$native_success], result$case[!result$native_success],
  result$native_failure_stage[!result$native_success]
))
names(fail_summary) <- c("anchor", "case", "native_failure_stage", "n")
write.csv(stage_summary, file.path(out, "pnb-oracle-stage-summary.csv"), row.names = FALSE)
write.csv(native_summary, file.path(out, "pnb-native-summary.csv"), row.names = FALSE)
write.csv(fail_summary, file.path(out, "pnb-native-failure-stage-summary.csv"), row.names = FALSE)
write.csv(stage_summary[stage_summary$anchor == "latent_quality", , drop = FALSE],
          file.path(out, "latent-oracle-stage-summary.csv"), row.names = FALSE)
write.csv(native_summary[native_summary$anchor == "latent_quality", , drop = FALSE],
          file.path(out, "latent-native-summary.csv"), row.names = FALSE)
write.csv(fail_summary[fail_summary$anchor == "latent_quality", , drop = FALSE],
          file.path(out, "latent-native-failure-stage-summary.csv"), row.names = FALSE)

max_or_na <- function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_
case_summary <- do.call(rbind, lapply(split(result, interaction(result$anchor, result$case)), function(z) {
  data.frame(anchor = z$anchor[[1L]], case = z$case[[1L]], n_draws = nrow(z),
    n_native_failures = sum(!z$native_success),
    n_pre_nonpositive = sum(z$pre_nonpositive), n_post_nonpositive = sum(z$post_nonpositive),
    max_abs_pre_foc_residual = max(z$max_abs_pre_foc_residual),
    max_abs_post_foc_residual = max(z$max_abs_post_foc_residual),
    max_success_pre_price_diff = max_or_na(z$max_abs_native_pre_price_diff),
    max_success_post_price_diff = max_or_na(z$max_abs_native_post_price_diff),
    max_success_pre_share_diff = max_or_na(z$max_abs_native_pre_share_diff),
    max_success_post_share_diff = max_or_na(z$max_abs_native_post_share_diff),
    max_success_price_change_diff = max_or_na(z$max_abs_native_price_change_diff),
    stringsAsFactors = FALSE)
}))
write.csv(case_summary, file.path(out, "pnb-oracle-case-summary.csv"), row.names = FALSE)
write.csv(case_summary[case_summary$anchor == "latent_quality", , drop = FALSE],
          file.path(out, "latent-oracle-case-summary.csv"), row.names = FALSE)

report <- file(file.path(out, "pnb-rate-classification-report.md"), open = "wt")
on.exit(close(report), add = TRUE)
writeLines(c(
  "# PNB real-price rate classification", "",
  "Generated by `reference_latent_rate_classification.R` from the archived CF/CS",
  "fits, saved policy plans, merger scenarios, and observed/latent simulation arrays.",
  "The script uses the permuted `rstan::extract` draw axis indexed by the saved",
  "plan indices. It does not call a production equilibrium solver, native target,",
  "bridge, Stan sampling, HMC, or nonlinear equilibrium evaluation.", "",
  "For each draw, latent reference shares are `(1-s0) * softmax(share_mean +",
  "saved_quality)`, with `delta = log(latent_reference_share/s0)`. The observed",
  "input FOC inversion fixes `value = p_observed + h_observed/alpha`; hence",
  "`log A = delta + h_observed`. The saved ordered-k role realization is used",
  "for pre and post ownership. The CF/CS real-rate equilibrium is solved through",
  "the scalar positive-real W branch, then shares are substituted into the FOCs.", "",
  "## Case summary", "", capture.output(print(case_summary, row.names = FALSE)), "",
  "## Oracle stage counts", "", capture.output(print(stage_summary, row.names = FALSE)), "",
  "## Native status counts", "", capture.output(print(native_summary, row.names = FALSE)), "",
  "## Native failure stage counts", "", capture.output(print(fail_summary, row.names = FALSE)), "",
  "## Every saved native failure", ""
), report)
fail_rows <- result[!result$native_success, , drop = FALSE]
if (nrow(fail_rows)) {
  writeLines(capture.output(print(fail_rows[, c(
    "anchor", "case", "draw_position", "draw_index", "k", "native_failure",
    "oracle_stage", "native_failure_stage", "min_pre_price", "min_pre_product",
    "min_post_price", "min_post_product", "max_abs_pre_foc_residual",
    "max_abs_post_foc_residual"), drop = FALSE], row.names = FALSE)), report)
} else {
  writeLines("No saved native failures.", report)
}
writeLines(c("", "The latent ledger has 400 rows in `latent-rate-classification.csv`;",
  "the combined observed/latent ledger has 800 rows in `pnb-rate-classification.csv`.",
  "`latent-failures.csv` and `observed-failures.csv` retain every saved native failure.",
  "Successful saved",
  "price and share differences are in the ledger and case summary. Residual",
  "bounds include all oracle states, including nonpositive real-rate states."), report)

print(case_summary)
print(stage_summary)
print(fail_summary)
print(fail_rows[, c("anchor", "case", "draw_position", "draw_index", "oracle_stage",
                    "native_failure_stage", "min_pre_price", "min_post_price"),
               drop = FALSE])
