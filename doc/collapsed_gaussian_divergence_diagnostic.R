# Post hoc diagnostic for the saved collapsed-Gaussian two-market ladder fit.
# This script reads stored draws only; it does not sample, refit, or modify a
# model.  It is intentionally exploratory and does not turn ten divergences
# into an architectural conclusion.

fit_path <- "doc/collapsed_gaussian_runs/20260916/cga-ladder-benchmark-20260916/stanfit/logit_bertrand_noisy_271828__B__ladder_multi_market.rds"
csv_path <- "doc/collapsed_gaussian_divergence_diagnostic.csv"
txt_path <- "doc/collapsed_gaussian_divergence_diagnostic.txt"

fit <- readRDS(fit_path)
sf <- fit@stanfit
draw_array <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE)
param_names <- dimnames(draw_array)[[3]]

wanted_patterns <- c(
  "^sigma_b_strat$",
  "^sigma_share$",
  "^s0_prob\\[[0-9]+\\]$",
  "^alpha_structural_market\\[[0-9]+\\]$",
  "^b_market_raw\\[[0-9]+\\]$",
  "^b_market\\[[0-9]+\\]$",
  "^mu_b_market_vec\\[[0-9]+\\]$",
  "^sigma_b_market_vec\\[[0-9]+\\]$"
)
keep <- unique(unlist(lapply(wanted_patterns, function(pattern) {
  grep(pattern, param_names, value = TRUE)
})))
if (length(keep) == 0L) {
  stop("No requested parameters were found in the saved fit")
}

post <- do.call(rbind, lapply(seq_len(dim(draw_array)[2]), function(chain) {
  out <- as.data.frame(draw_array[, chain, match(keep, param_names), drop = FALSE])
  names(out) <- keep
  out$chain <- chain
  out$iter <- seq_len(nrow(out))
  out$divergent <- as.integer(
    rstan::get_sampler_params(sf, inc_warmup = FALSE)[[chain]][, "divergent__"] > 0
  )
  out
}))

qfun <- function(x) {
  if (!length(x)) return(c(q05 = NA_real_, median = NA_real_, q95 = NA_real_))
  unname(stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE, type = 8))
}

records <- list()
add_record <- function(...) records[[length(records) + 1L]] <<- data.frame(...)

div_count <- aggregate(divergent ~ chain, post, sum)
names(div_count)[2] <- "divergences"
for (i in seq_len(nrow(div_count))) {
  add_record(
    record_type = "chain_divergence", parameter = NA_character_,
    scope = "postwarmup", chain = div_count$chain[i], iter = NA_integer_,
    n = nrow(post[post$chain == div_count$chain[i], , drop = FALSE]),
    q05 = NA_real_, median = NA_real_, q95 = NA_real_, mean = NA_real_,
    sd = NA_real_, value = NA_real_, rank_percentile = NA_real_,
    divergence_count = div_count$divergences[i]
  )
}

for (parameter in keep) {
  x <- post[[parameter]]
  for (scope in c("all", "divergent")) {
    z <- if (scope == "all") x else x[post$divergent == 1L]
    qq <- qfun(z)
    add_record(
      record_type = "parameter_summary", parameter = parameter,
      scope = scope, chain = NA_integer_, iter = NA_integer_, n = length(z),
      q05 = qq[1], median = qq[2], q95 = qq[3],
      mean = if (length(z)) mean(z) else NA_real_,
      sd = if (length(z) > 1L) stats::sd(z) else NA_real_,
      value = NA_real_, rank_percentile = NA_real_, divergence_count = NA_real_
    )
  }

  div_rows <- which(post$divergent == 1L)
  for (row in div_rows) {
    chain_values <- x[post$chain == post$chain[row]]
    rank_pct <- (sum(chain_values <= x[row]) - 0.5 * sum(chain_values == x[row])) /
      length(chain_values)
    add_record(
      record_type = "divergent_draw", parameter = parameter,
      scope = "divergent", chain = post$chain[row], iter = post$iter[row], n = 1L,
      q05 = NA_real_, median = NA_real_, q95 = NA_real_, mean = NA_real_, sd = NA_real_,
      value = x[row], rank_percentile = rank_pct, divergence_count = NA_real_
    )
  }
}

out <- do.call(rbind, records)
utils::write.csv(out, csv_path, row.names = FALSE, na = "")

div_rows <- post[post$divergent == 1L, c("chain", "iter", keep), drop = FALSE]
point_biserial <- vapply(keep, function(parameter) {
  suppressWarnings(stats::cor(post$divergent, post[[parameter]]))
}, numeric(1))

div_corr <- if (nrow(div_rows) >= 3L && length(keep) >= 2L) {
  suppressWarnings(stats::cor(div_rows[, keep, drop = FALSE]))
} else {
  matrix(NA_real_, nrow = length(keep), ncol = length(keep),
         dimnames = list(keep, keep))
}

# The saved ladder has two market-years per market: mky 1--2 belong to market
# 1 and mky 3--4 belong to market 2.  These all-draw correlations are separate
# from the small, divergent-only correlation matrix above.
market_s0_pairs <- list(
  c("b_market[1]", "s0_prob[1]"), c("b_market[1]", "s0_prob[2]"),
  c("b_market[2]", "s0_prob[3]"), c("b_market[2]", "s0_prob[4]")
)
all_corr <- data.frame(
  relation = vapply(market_s0_pairs, function(pair) paste(pair, collapse = " ~ "), character(1)),
  correlation = vapply(market_s0_pairs, function(pair) {
    stats::cor(post[[pair[1]]], stats::qlogis(post[[pair[2]]]))
  }, numeric(1)),
  stringsAsFactors = FALSE
)
all_corr <- rbind(
  all_corr,
  data.frame(
    relation = c(
      "mu_b_market_vec[1] ~ sigma_b_strat",
      "mu_b_market_vec[1] ~ sigma_share",
      "mu_b_market_vec[1] ~ sigma_b_market_vec[1]"
    ),
    correlation = c(
      stats::cor(post[["mu_b_market_vec[1]"]], post[["sigma_b_strat"]]),
      stats::cor(post[["mu_b_market_vec[1]"]], post[["sigma_share"]]),
      stats::cor(post[["mu_b_market_vec[1]"]], post[["sigma_b_market_vec[1]"]])
    ),
    stringsAsFactors = FALSE
  )
)

capture.output({
  cat("Collapsed-Gaussian saved-draw divergence diagnostic\n")
  cat("Fit:", fit_path, "\n")
  cat("Stored post-warmup draws:", nrow(post), "across", length(unique(post$chain)), "chains\n")
  cat("Total divergences:", sum(post$divergent), "\n\n")
  cat("Chain divergence counts\n")
  print(div_count, row.names = FALSE)
  cat("\nRequested parameters found:\n")
  cat(paste(keep, collapse = ", "), "\n\n")
  cat("All-draw and divergent-draw summaries are in the CSV. Divergent rank is\n")
  cat("the within-chain percentile of each divergent value (ties mid-ranked).\n\n")
  cat("Point-biserial correlations with the divergence indicator\n")
  print(round(point_biserial, 4))
  cat("\nCorrelation matrix among the ten divergent draws (exploratory; n=", nrow(div_rows), ")\n", sep = "")
  print(round(div_corr, 3))
  cat("\nAll-draw correlations for market intercepts and matched outside shares\n")
  print(data.frame(relation = all_corr$relation,
                   correlation = round(all_corr$correlation, 4)), row.names = FALSE)
  cat("\nThese use the stored market-year mapping (mky 1--2 in market 1;\n")
  cat("mky 3--4 in market 2). They describe posterior dependence in this\n")
  cat("failed fit; they are not evidence of a required reparameterization.\n")
  cat("\nInterpretation guard:\n")
  cat("This is a descriptive analysis of stored draws from a fit that failed the\n")
  cat("diagnostic gate. It does not establish a causal source of divergences and\n")
  cat("does not evaluate the new s0 auto preset: this fit used explicit unpooled\n")
  cat("controls. No refit or sampler tuning was performed.\n")
}, file = txt_path)

message("Wrote ", csv_path, " and ", txt_path)
