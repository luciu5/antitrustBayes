#!/usr/bin/env Rscript

## Prespecified collapsed-Gaussian benchmark runner.
##
## The runner is intentionally self-contained.  The DGP is an independent
## base-R oracle for the current global centered quality and standardized
## price equations.  Sampling is preparation-only until bayes_est() exposes
## the share_likelihood and quality_parameterization arguments.  The default
## invocation is a dry run; --run is required before any fit is attempted.

cga_or_null <- function(x, y) if (is.null(x)) y else x

cga_runner_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[startsWith(a, "--file=")]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]),
                                      mustWork = FALSE)) else getwd()
})

cga_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript collapsed_gaussian_benchmark.R [--run [--resume] --out=/tmp/cga-run]",
    "",
    "The default prints and validates the frozen 21-fit primary plan.",
    "--run attempts each fit exactly once and requires --out under /tmp.",
    "Fit, log, and RDS artifacts are written only below that caller-supplied path.",
    "The four-period ladder and compact multi-market fixture are validated in",
    "the dry run and are outside the 21 primary cells.",
    sep = "\n"
  ), "\n")
}

cga_option <- function(args, name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (!length(hit)) return(default)
  sub(key, "", hit[[1L]], fixed = TRUE)
}

cga_parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("--help", "-h"))) return(list(help = TRUE))
  unknown <- args[startsWith(args, "--") &
    !args %in% c("--run", "--resume", "--help", "-h") &
    !startsWith(args, "--out=")]
  if (length(unknown)) stop("Unknown option: ", unknown[[1L]])
  if ("--resume" %in% args && !"--run" %in% args)
    stop("--resume requires --run.")
  list(
    help = FALSE,
    run = any(args == "--run"),
    resume = any(args == "--resume"),
    out = cga_option(args, "out")
  )
}

cga_require_tmp <- function(path, resume = FALSE) {
  if (is.null(path) || !nzchar(path)) stop("--run requires --out=/tmp/...")
  path <- normalizePath(path, mustWork = FALSE)
  if (!identical(path, "/tmp") && !startsWith(path, "/tmp/")) {
    stop("--out must be a caller-supplied directory under /tmp")
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Unable to create output directory: ", path)
  exists_run <- file.exists(file.path(path, "benchmark_run_manifest.csv"))
  if (exists_run && !resume)
    stop("An existing benchmark run must be preserved; choose a new --out directory.")
  if (resume && !exists_run) stop("--resume requires an existing benchmark manifest.")
  path
}

cga_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(as.integer(seed))
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(expr)
}

cga_derived_seed <- function(seed, ...) {
  bits <- paste(..., collapse = "|")
  bytes <- utf8ToInt(bits)
  z <- as.double(seed) * 104729 +
    sum(as.double(seq_along(bytes)) * as.double(bytes))
  as.integer(1 + z %% 2147483646)
}

cga_lse <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cga_softmax <- function(x) {
  z <- x - cga_lse(x)
  out <- exp(z)
  out / sum(out)
}

cga_helmert <- function(n_group) {
  out <- matrix(0, nrow = n_group, ncol = max(0L, n_group - 1L))
  if (n_group > 1L) {
    for (k in seq_len(n_group - 1L)) {
      out[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
      out[k + 1L, k] <- -k / sqrt(k * (k + 1))
    }
  }
  out
}

cga_group_index <- function(group) {
  lev <- unique(as.character(group))
  match(as.character(group), lev)
}

cga_safe_sd <- function(x) {
  z <- stats::sd(as.numeric(x))
  if (!is.finite(z) || z <= 0) stop("DGP standardization has zero variance")
  z
}

cga_supply_name <- function(x) {
  x <- match.arg(x, c("bertrand", "cournot"))
  x
}

cga_role_conduct_logit <- function(firm_share, s0, supply_model,
                                    structure = "simultaneous", rank = NULL,
                                    k = NULL) {
  supply_model <- cga_supply_name(supply_model)
  n <- length(firm_share)
  core <- rep(TRUE, n)
  if (identical(structure, "fringe")) core <- rank <= k
  if (identical(structure, "stackelberg")) core <- rank <= k
  S <- firm_share * (1 - s0)
  if (identical(structure, "simultaneous")) {
    return(if (supply_model == "bertrand") 1 - S else s0 / (s0 + S))
  }
  if (identical(structure, "fringe")) {
    if (supply_model == "bertrand") return(ifelse(core, 1 - S, 1))
    return(ifelse(core, s0 / (s0 + S), 1))
  }
  if (supply_model == "bertrand") {
    follower <- !core
    B_F <- sum(S[follower]^2 / (1 - S[follower] + S[follower]^2))
    phi <- 1 / (1 - B_F)
    return(ifelse(core, 1 - phi * S, 1 - S))
  }
  S_F <- sum(S[!core])
  ifelse(core, (s0 + S_F) / (s0 + S_F + S), s0 / (s0 + S))
}

cga_role_inverse_lerner_ces <- function(firm_rev_share, R0, gamma,
                                          supply_model,
                                          structure = "simultaneous",
                                          rank = NULL, k = NULL) {
  supply_model <- cga_supply_name(supply_model)
  core <- rep(TRUE, length(firm_rev_share))
  if (identical(structure, "fringe") || identical(structure, "stackelberg"))
    core <- rank <= k
  R <- firm_rev_share * (1 - R0)
  h <- gamma - 1
  H <- 1 + h * R0
  b <- h / (gamma * H)
  if (identical(structure, "simultaneous")) {
    if (supply_model == "bertrand") return(gamma - h * R)
    L <- 1 / gamma + (h * R) / (gamma * (1 + h * R0))
    return(1 / L)
  }
  if (identical(structure, "fringe")) {
    if (supply_model == "bertrand") return(ifelse(core, gamma - h * R, gamma))
    L <- ifelse(core, 1 / gamma + b * R, 1 / gamma)
    return(1 / L)
  }
  if (supply_model == "bertrand") {
    follower <- !core
    t <- gamma - h * R
    a <- h^2 * R / (t * (t - 1) + h^2 * R)
    B_F <- sum(R[follower] * a[follower])
    phi <- 1 / (1 - B_F)
    return(ifelse(core, gamma - h * phi * R, gamma - h * R))
  }
  follower <- !core
  a <- h * (H - 2 * R + gamma * h * R0 * R / H) /
    (H * (H + (gamma - 2) * R))
  C_F <- sum(R[follower] * a[follower])
  L <- ifelse(core, 1 / gamma + b * R / (1 + C_F), 1 / gamma + b * R)
  1 / L
}

cga_make_plan <- function() {
  ## The order is frozen so the six groups are easy to audit.  A/B/C share
  ## the same rows within each data design; only their likelihood/API fields
  ## differ.
  add_modes <- function(x) {
    do.call(rbind, lapply(c("A", "B", "C"), function(mode) {
      y <- x
      y$mode <- mode
      y$share_likelihood <- if (mode == "C") "student_t" else "gaussian"
      y$margin_likelihood <- if (mode == "C") "student_t" else "gaussian"
      y$quality_parameterization <- if (mode == "B") "collapsed" else "sampled"
      y$bridge_requested <- identical(x$design, "logit_bertrand_noisy_271828")
      y
    }))
  }
  cells <- list(
    data.frame(design = "logit_bertrand_exact_271828", demand_model = "logit",
               supply_model = "bertrand", structure = "simultaneous",
               share_scenario = "exact", seed = 271828L,
               partition_prior = "none", stringsAsFactors = FALSE),
    data.frame(design = "logit_bertrand_noisy_271828", demand_model = "logit",
               supply_model = "bertrand", structure = "simultaneous",
               share_scenario = "noisy", seed = 271828L,
               partition_prior = "none", stringsAsFactors = FALSE),
    data.frame(design = "logit_bertrand_noisy_314159", demand_model = "logit",
               supply_model = "bertrand", structure = "simultaneous",
               share_scenario = "noisy", seed = 314159L,
               partition_prior = "none", stringsAsFactors = FALSE),
    data.frame(design = "logit_cournot_noisy_271828", demand_model = "logit",
               supply_model = "cournot", structure = "simultaneous",
               share_scenario = "noisy", seed = 271828L,
               partition_prior = "none", stringsAsFactors = FALSE),
    data.frame(design = "logit_BF_noisy_271828", demand_model = "logit",
               supply_model = "bertrand", structure = "fringe",
               share_scenario = "noisy", seed = 271828L,
               partition_prior = "ordered_k", stringsAsFactors = FALSE),
    data.frame(design = "logit_BS_noisy_271828", demand_model = "logit",
               supply_model = "bertrand", structure = "stackelberg",
               share_scenario = "noisy", seed = 271828L,
               partition_prior = "ordered_k", stringsAsFactors = FALSE),
    data.frame(design = "ces_bertrand_noisy_271828", demand_model = "ces",
               supply_model = "bertrand", structure = "simultaneous",
               share_scenario = "noisy", seed = 271828L,
               partition_prior = "none", stringsAsFactors = FALSE)
  )
  rows <- do.call(rbind, lapply(cells, add_modes))
  rows$k_support <- ifelse(rows$structure == "fringe", "0:5",
                           ifelse(rows$structure == "stackelberg", "1:5", NA_character_))
  rows$k_prior <- ifelse(rows$structure == "fringe", "uniform_1/6",
                         ifelse(rows$structure == "stackelberg", "uniform_1/5", NA_character_))
  rows$k_true <- ifelse(rows$structure %in% c("fringe", "stackelberg"), 2L, NA_integer_)
  rows$demand <- rows$demand_model
  rows$conduct <- rows$supply_model
  rows$use_rho <- FALSE
  rows$use_market_resid_scale <- FALSE
  rows$use_alpha_market_pool <- FALSE
  rows$use_s0_market_pool <- FALSE
  rows$use_s0_market_year_pool <- FALSE
  rows$use_centered_alpha <- FALSE
  rows$use_centered_s0_pool <- FALSE
  rows$use_hmt <- FALSE
  rows$df <- 5
  rows$chains <- 4L
  rows$iter <- 1000L
  rows$warmup <- 500L
  rows$adapt_delta <- .95
  rows$max_treedepth <- 13L
  rows$metric <- "diag_e"
  rows$run_profile <- "primary"
  rows$case_id <- paste(rows$design, rows$mode, sep = "__")
  rows <- rows[order(match(rows$design, vapply(cells, `[[`, character(1), "design")),
                     match(rows$mode, c("A", "B", "C"))), , drop = FALSE]
  rownames(rows) <- NULL
  if (nrow(rows) != 21L) stop("Frozen primary plan must contain exactly 21 fits")
  rows
}

cga_base_prices <- function(seed, n_products = 24L, n_markets = 1L) {
  cga_with_seed(cga_derived_seed(seed, "prices", n_products, n_markets), {
    base <- matrix(9.8 + stats::runif(n_products * n_markets, 0, 0.4),
                   nrow = n_markets, ncol = n_products)
    dimnames(base) <- list(paste0("market_", seq_len(n_markets)),
                           paste0("firm_", sprintf("%02d", seq_len(n_products))))
    base
  })
}

cga_make_dgp <- function(seed = 271828L,
                         demand_model = c("logit", "ces"),
                         supply_model = c("bertrand", "cournot"),
                         structure = c("simultaneous", "fringe", "stackelberg"),
                         share_scenario = c("exact", "noisy"),
                         n_markets = 1L, n_periods = 1L,
                         outside_price = 1.3,
                         truth_k = 2L,
                         n_products = 24L,
                         persistent_groups = TRUE) {
  demand_model <- match.arg(demand_model)
  supply_model <- cga_supply_name(supply_model)
  structure <- match.arg(structure)
  share_scenario <- match.arg(share_scenario)
  n_markets <- as.integer(n_markets); n_periods <- as.integer(n_periods)
  n_products <- as.integer(n_products)
  if (n_markets < 1L || n_periods < 1L || n_products < 2L)
    stop("n_markets, n_periods, and n_products must be positive")
  if (!isTRUE(persistent_groups)) stop("The benchmark requires persistent groups")
  if (structure != "simultaneous" &&
      (length(truth_k) != 1L || truth_k < 1L || truth_k > n_products))
    stop("truth_k must be feasible for exact role designs")
  outside_price <- rep(as.numeric(outside_price), length.out = n_markets * n_periods)
  if (any(!is.finite(outside_price) | outside_price <= 0))
    stop("outside_price must be positive and finite")

  base <- cga_base_prices(seed, n_products, n_markets)
  n_group <- n_markets * n_products
  group_levels <- if (n_markets == 1L) {
    paste0("firm_", sprintf("%02d", seq_len(n_products)))
  } else {
    unlist(lapply(seq_len(n_markets), function(m) paste0(
      "market_", m, "_firm_", sprintf("%02d", seq_len(n_products))
    )), use.names = FALSE)
  }
  ## The global Helmert coordinate is the quality prior used by Stan.  Each
  ## market receives its own group identifiers in the multi-market fixture,
  ## while a group's product rows repeat across periods.
  quality_draw <- cga_with_seed(cga_derived_seed(seed, "quality", n_group), {
    C <- cga_helmert(n_group)
    beta <- stats::rnorm(max(0L, n_group - 1L), 0, .15)
    list(C = C, beta = beta,
         group = if (n_group > 1L) as.numeric(C %*% beta) else 0)
  })
  quality_group <- quality_draw$group
  names(quality_group) <- group_levels

  rec <- vector("list", n_markets * n_periods)
  truths <- vector("list", n_markets * n_periods)
  at <- 0L
  truth_records <- 0L
  all_prices <- numeric(0L)
  for (m in seq_len(n_markets)) for (tt in seq_len(n_periods)) {
    at <- at + 1L
    ## Fixed, deterministic period/market price movement makes the repeated
    ## rows informative without changing the demand or conduct equations.
    p <- base[m, ] + .18 * (tt - 1L) + .27 * (m - 1L) +
      .08 * sin(seq_len(n_products) + tt + m)
    all_prices <- c(all_prices, p)
    idx <- ((m - 1L) * n_products + 1L):(m * n_products)
    groups <- group_levels[idx]
    q <- unname(quality_group[groups])
    market_name <- paste0("market_", m)
    year_name <- as.character(tt)
    rec[[at]] <- data.frame(
      prices = p, quantities = NA_real_, margins = NA_real_,
      market = market_name, year = year_name, group = groups,
      product = paste0(market_name, "_product_", sprintf("%02d", seq_len(n_products))),
      outside_price = outside_price[at], stringsAsFactors = FALSE
    )
  }
  price_all <- unlist(lapply(rec, `[[`, "prices"), use.names = FALSE)
  price_sd <- cga_safe_sd(price_all)
  price_std <- (price_all - mean(price_all)) / price_sd
  log_price_all <- log(price_all)
  log_price_sd <- cga_safe_sd(log_price_all)
  log_price_std <- (log_price_all - mean(log_price_all)) / log_price_sd
  desired_alpha <- 6
  desired_gamma <- 3.5
  a_std <- if (demand_model == "logit") desired_alpha * price_sd else
    (desired_gamma - 1) * log_price_sd
  gamma <- if (demand_model == "logit") NA_real_ else desired_gamma
  row_pos <- 0L
  for (m in seq_len(n_markets)) for (tt in seq_len(n_periods)) {
    truth_records <- truth_records + 1L
    row_idx <- ((truth_records - 1L) * n_products + 1L):(truth_records * n_products)
    market_name <- paste0("market_", m)
    year_name <- as.character(tt)
    q <- unname(quality_group[group_levels[((m - 1L) * n_products + 1L):(m * n_products)]])
    p <- price_all[row_idx]
    std <- if (demand_model == "logit") price_std[row_idx] else log_price_std[row_idx]
    p0 <- outside_price[truth_records]
    delta <- if (demand_model == "logit") {
      -a_std * std + q
    } else {
      q - (a_std + log_price_sd) * std + gamma * log(p0)
    }
    ## The single-market zero-intercept equation identifies logit(s0) from
    ## the centered mean: mean(log share) = logit(s0), since mean(q)=0 and
    ## mean(standardized price)=0.  This is the exact multi-market/year
    ## analogue with b_market fixed to zero in the oracle.
    s0 <- stats::plogis(-cga_lse(delta))
    latent_share <- cga_softmax(delta)
    observed_share <- if (share_scenario == "exact") latent_share else
      cga_with_seed(cga_derived_seed(seed, "share", m, tt), {
        cga_softmax(delta + stats::rnorm(n_products, 0, .04))
      })
    inside_price_index <- sum(p * observed_share)
    revenue_share <- p * observed_share / sum(p * observed_share)
    firm_metric <- if (demand_model == "logit") observed_share else revenue_share
    rank_order <- order(-firm_metric, group_levels[((m - 1L) * n_products + 1L):(m * n_products)])
    rank <- integer(n_products); rank[rank_order] <- seq_len(n_products)
    if (structure == "simultaneous") role_core <- rep(TRUE, n_products) else
      role_core <- rank <= truth_k
    R0 <- if (demand_model == "ces") stats::plogis(
      stats::qlogis(s0) + log(p0) - log(inside_price_index)
    ) else NA_real_
    inv_lerner <- if (demand_model == "logit") {
      p * cga_role_conduct_logit(observed_share, s0, supply_model, structure,
                                rank, truth_k) * desired_alpha
    } else {
      cga_role_inverse_lerner_ces(revenue_share, R0, gamma, supply_model,
                                  structure, rank, truth_k)
    }
    ## Log inverse-margin noise is Gaussian for every A/B/C row.  C changes
    ## the fitted likelihood family only; its observed rows are identical.
    margin_noise <- cga_with_seed(cga_derived_seed(seed, "margin", m, tt),
                                  stats::rnorm(n_products, 0, .04))
    margin_truth <- 1 / inv_lerner
    observed_margin <- exp(-(log(inv_lerner) + margin_noise))
    if (any(!is.finite(observed_margin) | observed_margin <= 0))
      stop("DGP produced invalid margins")
    if (demand_model %in% c("logit", "ces") && any(observed_margin >= 1))
      stop("DGP produced an output Lerner margin >= 1")
    rec[[truth_records]]$quantities <- 1000 * observed_share
    rec[[truth_records]]$margins <- observed_margin
    truths[[truth_records]] <- list(
      market = market_name, year = year_name, mky = truth_records,
      prices = p, quality = q, quality_centered = q,
      quality_group = quality_group[group_levels[((m - 1L) * n_products + 1L):(m * n_products)]],
      price_std = std, price_sd = if (demand_model == "logit") price_sd else log_price_sd,
      price_mean = if (demand_model == "logit") mean(price_all) else mean(log_price_all),
      stan_mean = stats::qlogis(s0) + delta,
      latent_log_score = delta, latent_log_share = log(latent_share),
      latent_share_in = latent_share, observed_share_in = observed_share,
      s0 = s0, outside_price = p0, inside_price_index = inside_price_index,
      outside_revenue = R0, revenue_share_in = revenue_share,
      role_metric = firm_metric, role_rank = rank,
      role_core = role_core, role_labels = if (structure == "stackelberg")
        ifelse(role_core, "leader", "follower") else if (structure == "fringe")
          ifelse(role_core, "strategic", "fringe") else rep("simultaneous", n_products),
      truth_k = if (structure == "simultaneous") NA_integer_ else as.integer(truth_k),
      inverse_margin_truth = inv_lerner, margin_truth = margin_truth,
      margin_noise = margin_noise, inverse_margin_observed = 1 / observed_margin
    )
  }
  ## Exact BF/BS uses one observable firm ranking for the complete market
  ## panel.  Repeated-period fixtures therefore rank by equal-period mean
  ## conditional share (Logit) or conditional revenue share (CES), then
  ## recompute every hard-role FOC using that fixed rank.
  if (structure != "simultaneous" && n_periods > 1L) {
    for (m in seq_len(n_markets)) {
      jj <- which(vapply(truths, function(x) identical(x$market,
        paste0("market_", m)), logical(1)))
      metric_matrix <- do.call(rbind, lapply(truths[jj], `[[`, "role_metric"))
      labels <- names(truths[[jj[[1L]]]]$quality_group)
      mean_metric <- colMeans(metric_matrix)
      fixed_order <- order(-mean_metric, labels)
      fixed_rank <- integer(n_products)
      fixed_rank[fixed_order] <- seq_len(n_products)
      for (j in jj) {
        tr <- truths[[j]]
        tr$role_rank <- fixed_rank
        tr$role_core <- fixed_rank <= truth_k
        tr$role_labels <- if (structure == "stackelberg")
          ifelse(tr$role_core, "leader", "follower") else
            ifelse(tr$role_core, "strategic", "fringe")
        inv <- if (demand_model == "logit") {
          tr$prices * cga_role_conduct_logit(tr$observed_share_in, tr$s0,
            supply_model, structure, fixed_rank, truth_k) * desired_alpha
        } else {
          cga_role_inverse_lerner_ces(tr$revenue_share_in,
            tr$outside_revenue, gamma, supply_model, structure,
            fixed_rank, truth_k)
        }
        tr$inverse_margin_truth <- inv
        tr$margin_truth <- 1 / inv
        tr$inverse_margin_observed <- exp(log(inv) + tr$margin_noise)
        rec[[j]]$margins <- exp(-(log(inv) + tr$margin_noise))
        truths[[j]] <- tr
      }
    }
  }
  dat <- do.call(rbind, rec)
  rownames(dat) <- NULL
  list(
    data = dat,
    truth = truths,
    quality_group = quality_group,
    demand_model = demand_model,
    supply_model = supply_model,
    structure = structure,
    share_scenario = share_scenario,
    seed = as.integer(seed),
    parameters = list(
      alpha = if (demand_model == "logit") desired_alpha else NA_real_,
      gamma = if (demand_model == "ces") desired_gamma else NA_real_,
      a_std = a_std, price_sd = price_sd, log_price_sd = log_price_sd,
      outside_price = outside_price, share_noise_sd = if (share_scenario == "noisy") .04 else 0,
      margin_noise_sd = .04, margin_noise_family = "gaussian",
      margin_noise_df = NA_real_, quality_prior = "global_helmert_centered",
      quality_basis = quality_draw$C, quality_coordinates = quality_draw$beta,
      persistent_groups = TRUE, truth_k = if (structure == "simultaneous") NA_integer_ else truth_k
    ),
    normalization = list(
      global_quality_center_error = abs(mean(quality_group)),
      repeated_quality = all(vapply(truths, function(x) max(abs(x$quality_group -
        quality_group[names(x$quality_group)])), numeric(1)) < 1e-14)
    ),
    design = list(n_markets = n_markets, n_periods = n_periods,
                  n_products = n_products, group_levels = group_levels)
  )
}

cga_observed_shares <- function(dat) {
  split_idx <- interaction(dat$market, dat$year, drop = TRUE, lex.order = TRUE)
  out <- numeric(nrow(dat))
  for (z in levels(split_idx)) {
    ii <- which(split_idx == z)
    out[ii] <- dat$quantities[ii] / sum(dat$quantities[ii])
  }
  out
}

cga_validate_dgp <- function(panel, tolerance = 2e-10) {
  if (!is.list(panel) || !is.data.frame(panel$data) || !length(panel$truth))
    stop("DGP must contain data and truth records")
  dat <- panel$data
  if (any(!is.finite(dat$prices)) || any(dat$prices <= 0) ||
      any(!is.finite(dat$quantities)) || any(dat$quantities <= 0) ||
      any(!is.finite(dat$margins)) || any(dat$margins <= 0))
    stop("DGP observations are not finite and positive")
  if (any(dat$margins >= 1)) stop("Output DGP margins must be below one")
  if (length(panel$truth) != length(unique(interaction(dat$market, dat$year))))
    stop("DGP truth does not have one record per market-year")
  if (!isTRUE(panel$parameters$persistent_groups))
    stop("DGP quality groups must persist across periods")
  q_by_group <- split(
    unlist(lapply(panel$truth, `[[`, "quality_group"), use.names = FALSE),
    unlist(lapply(panel$truth, function(x) names(x$quality_group)),
           use.names = FALSE)
  )
  if (max(vapply(q_by_group, function(x) max(x) - min(x), numeric(1))) > tolerance)
    stop("Quality is not persistent within group")
  if (abs(mean(unlist(lapply(panel$truth, function(x) x$quality_group)))) > tolerance)
    stop("Global centered quality is not centered")
  if (ncol(panel$parameters$quality_basis) > 0L &&
      max(abs(as.numeric(panel$parameters$quality_basis %*%
        panel$parameters$quality_coordinates) -
        as.numeric(panel$quality_group))) > tolerance)
    stop("Helmert quality coordinates are inconsistent")
  obs_share <- cga_observed_shares(dat)
  p_sd <- cga_safe_sd(dat$prices)
  lp <- log(dat$prices)
  lp_sd <- cga_safe_sd(lp)
  for (j in seq_along(panel$truth)) {
    tr <- panel$truth[[j]]
    ii <- which(dat$market == tr$market & dat$year == tr$year)
    if (length(ii) != length(tr$prices)) stop("Truth/data market-year mismatch")
    p <- dat$prices[ii]
    std <- if (panel$demand_model == "logit") (p - mean(dat$prices)) / p_sd else
      (log(p) - mean(lp)) / lp_sd
    if (max(abs(std - tr$price_std)) > tolerance)
      stop("DGP does not use actual global standardized prices")
    if (panel$demand_model == "logit") {
      expected_mean <- stats::qlogis(tr$s0) + tr$quality -
        panel$parameters$a_std * std
    } else {
      expected_mean <- stats::qlogis(tr$s0) + tr$quality -
        panel$parameters$gamma * (log(p) - mean(lp)) +
        panel$parameters$gamma * log(tr$outside_price)
    }
    if (max(abs(expected_mean - tr$stan_mean)) > tolerance)
      stop("Stan mean reconstruction failed")
    expected_latent <- cga_softmax(expected_mean - stats::qlogis(tr$s0))
    if (max(abs(expected_latent - tr$latent_share_in)) > tolerance)
      stop("Latent share reconstruction failed")
    if (max(abs(obs_share[ii] - tr$observed_share_in)) > tolerance)
      stop("Observed composition does not match quantities")
    if (abs(sum(tr$observed_share_in) - 1) > tolerance)
      stop("Observed shares are not compositional")
    if (panel$share_scenario == "exact" &&
        max(abs(tr$observed_share_in - tr$latent_share_in)) > tolerance)
      stop("Exact shares are not latent structural shares")
    ## Reconstruct the package's observable role metric and fixed rank.
    metric <- if (panel$demand_model == "logit") tr$observed_share_in else
      tr$revenue_share_in
    if (max(abs(metric - tr$role_metric)) > tolerance)
      stop("Role metric reconstruction failed")
    expected_rank <- integer(length(metric))
    if (panel$structure == "simultaneous") {
      expected_rank[order(-metric, dat$group[ii])] <- seq_along(metric)
    } else {
      same_market <- which(vapply(panel$truth, function(x)
        identical(x$market, tr$market), logical(1)))
      mean_metric <- colMeans(do.call(rbind, lapply(panel$truth[same_market],
        `[[`, "role_metric")))
      labels <- names(tr$quality_group)
      expected_rank[order(-mean_metric, labels)] <- seq_along(metric)
    }
    if (!identical(as.integer(expected_rank), as.integer(tr$role_rank)))
      stop("Role ranking reconstruction failed")
    if (panel$demand_model == "ces") {
      expected_r0 <- stats::plogis(stats::qlogis(tr$s0) +
        log(tr$outside_price) - log(sum(p * tr$observed_share_in)))
      if (abs(expected_r0 - tr$outside_revenue) > tolerance)
        stop("CES outside-price normalization reconstruction failed")
      inv <- cga_role_inverse_lerner_ces(tr$revenue_share_in,
        expected_r0, panel$parameters$gamma, panel$supply_model,
        panel$structure, tr$role_rank, tr$truth_k)
    } else {
      inv <- p * cga_role_conduct_logit(tr$observed_share_in, tr$s0,
        panel$supply_model, panel$structure, tr$role_rank, tr$truth_k) *
        panel$parameters$alpha
    }
    if (max(abs(inv - tr$inverse_margin_truth)) > tolerance)
      stop("FOC inverse-margin reconstruction failed")
    if (max(abs(1 / tr$margin_truth - tr$inverse_margin_truth)) > tolerance)
      stop("FOC margin truth is inconsistent")
    if (max(abs(log(1 / dat$margins[ii]) -
               (log(tr$inverse_margin_truth) + tr$margin_noise))) > tolerance)
      stop("Observed margin noise is inconsistent")
    if (any(!is.finite(tr$role_rank)) ||
        any(tr$role_rank < 1L | tr$role_rank > length(tr$role_rank)))
      stop("Invalid role rank")
  }
  if (panel$share_scenario == "noisy" &&
      !any(vapply(panel$truth, function(x) max(abs(x$observed_share_in -
        x$latent_share_in)) > 1e-12, logical(1))))
    stop("Noisy DGP did not perturb any composition")
  list(
    status = "passed",
    n_observations = nrow(dat), n_market_year = length(panel$truth),
    max_stan_mean_error = max(vapply(panel$truth, function(x) {
      ii <- which(dat$market == x$market & dat$year == x$year)
      p <- dat$prices[ii]
      if (panel$demand_model == "logit") exp <- stats::qlogis(x$s0) + x$quality -
        panel$parameters$a_std * x$price_std else exp <- stats::qlogis(x$s0) +
        x$quality - panel$parameters$gamma * (log(p) - mean(lp)) +
        panel$parameters$gamma * log(x$outside_price)
      max(abs(exp - x$stan_mean))
    }, numeric(1))),
    max_quality_center_error = abs(mean(unlist(lapply(panel$truth,
      `[[`, "quality_group")))),
    max_share_sum_error = max(vapply(panel$truth, function(x)
      abs(sum(x$observed_share_in) - 1), numeric(1))),
    max_foc_error = max(vapply(panel$truth, function(x)
      max(abs(1 / x$margin_truth - x$inverse_margin_truth)), numeric(1)))
  )
}

cga_ladder_fixtures <- function(seed = 271828L) {
  list(
    four_period = cga_make_dgp(seed = seed, demand_model = "logit",
      supply_model = "bertrand", structure = "simultaneous",
      share_scenario = "noisy", n_markets = 1L, n_periods = 4L),
    multi_market = cga_make_dgp(seed = cga_derived_seed(seed, "fixture"),
      demand_model = "logit", supply_model = "bertrand",
      structure = "simultaneous", share_scenario = "noisy",
      n_markets = 2L, n_periods = 2L)
  )
}

cga_validate_plan <- function(plan = cga_make_plan()) {
  required <- c("case_id", "mode", "demand_model", "supply_model",
                "structure", "share_scenario", "seed", "k_support", "k_prior",
                "use_rho", "use_market_resid_scale", "use_alpha_market_pool",
                "use_s0_market_pool", "use_s0_market_year_pool",
                "use_centered_alpha", "use_centered_s0_pool", "use_hmt",
                "chains", "iter", "warmup", "adapt_delta", "max_treedepth",
                "metric")
  if (!all(required %in% names(plan))) stop("Primary plan is missing fields")
  if (nrow(plan) != 21L) stop("Primary plan must contain exactly 21 fits")
  if (!identical(as.integer(plan$chains), rep(4L, 21L)) ||
      !identical(as.integer(plan$iter), rep(1000L, 21L)) ||
      !identical(as.integer(plan$warmup), rep(500L, 21L)) ||
      !identical(as.integer(plan$max_treedepth), rep(13L, 21L)) ||
      !all(plan$adapt_delta == .95) || !all(plan$metric == "diag_e"))
    stop("Frozen sampler settings changed")
  if (any(plan$use_rho | plan$use_market_resid_scale |
          plan$use_alpha_market_pool | plan$use_s0_market_pool |
          plan$use_s0_market_year_pool | plan$use_centered_alpha |
          plan$use_centered_s0_pool | plan$use_hmt))
    stop("Primary pooling/rho/HMT flags must be false")
  mode_counts <- table(plan$mode)
  if (!all(c("A", "B", "C") %in% names(mode_counts)) ||
      any(as.integer(mode_counts[c("A", "B", "C")]) != 7L))
    stop("Primary modes are not balanced A/B/C")
  if (sum(plan$structure == "fringe") != 3L ||
      sum(plan$structure == "stackelberg") != 3L)
    stop("Primary exact BF/BS cells are missing")
  invisible(TRUE)
}

cga_fit_args <- function(panel, case, settings = NULL) {
  if (is.null(settings)) settings <- list(
    chains = case$chains, iter = case$iter, warmup = case$warmup,
    cores = 1L, refresh = 0L,
    control = list(adapt_delta = case$adapt_delta,
                   max_treedepth = case$max_treedepth, metric = "diag_e")
  )
  args <- list(
    prices = panel$data$prices, quantities = panel$data$quantities,
    margins = panel$data$margins, market = panel$data$market,
    year = panel$data$year, group = panel$data$group,
    demand_model = case$demand_model, market_side = "output",
    supply_model = case$supply_model, structure = case$structure,
    use_rho = FALSE, df = 5,
    use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_s0_market_pool = FALSE,
    use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
    use_hmt = FALSE, prior_manifest = "baseline", likelihood = "joint",
    ## bayes_est() expects one outside price per market-year, whereas the
    ## observation table repeats that value on every product row.
    outside_price = as.numeric(panel$parameters$outside_price),
    share_likelihood = case$share_likelihood,
    quality_parameterization = case$quality_parameterization,
    margin_likelihood = case$margin_likelihood,
    stan_options = settings
  )
  if (case$structure == "fringe") {
    args$partition_prior <- "ordered_k"
    args$k_support <- 0:5
    args$k_prior <- rep(1 / 6, 6)
    args$k_scope <- "market"
  } else if (case$structure == "stackelberg") {
    args$partition_prior <- "ordered_k"
    args$k_support <- 1:5
    args$k_prior <- rep(1 / 5, 5)
    args$k_scope <- "market"
  }
  seed_bytes <- utf8ToInt(case$case_id)
  args$stan_options$seed <- as.integer(1 + (as.double(case$seed) * 104729 +
    sum(seq_along(seed_bytes) * seed_bytes)) %% 2147483646)
  args
}

cga_api_ready <- function(estimator = NULL) {
  if (!is.null(estimator)) return(TRUE)
  if (!requireNamespace("antitrustBayes", quietly = TRUE)) return(FALSE)
  f <- getExportedValue("antitrustBayes", "bayes_est")
  n <- names(formals(f))
  all(c("share_likelihood", "quality_parameterization") %in% n)
}

cga_unwrap_stanfit <- function(fit) {
  if (is.null(fit)) return(NULL)
  if (methods::is(fit, "BayesFit")) return(fit@stanfit)
  if (methods::is(fit, "stanfit")) return(fit)
  if (is.list(fit) && !is.null(fit$stanfit)) return(cga_unwrap_stanfit(fit$stanfit))
  NULL
}

cga_extract <- function(stanfit, name) {
  if (is.null(stanfit) || !requireNamespace("rstan", quietly = TRUE)) return(NULL)
  tryCatch(rstan::extract(stanfit, pars = name, permuted = TRUE)[[1L]],
           error = function(e) NULL)
}

cga_draw_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  if (is.null(d)) return(matrix(as.numeric(x), nrow = length(x), ncol = 1L))
  if (!length(d) || d[[1L]] < 1L) return(NULL)
  matrix(as.numeric(x), nrow = d[[1L]], ncol = prod(d[-1L]))
}

cga_first_ess <- function(summary_table, patterns, column = c("n_eff", "ess_bulk")) {
  if (is.null(summary_table) || !nrow(summary_table)) return(NA_real_)
  column <- column[[1L]]
  hit <- grepl(paste(patterns, collapse = "|"), rownames(summary_table))
  if (!any(hit)) return(NA_real_)
  cols <- intersect(c(column, if (column == "n_eff") "ess_bulk" else "n_eff"),
                    colnames(summary_table))
  if (!length(cols)) return(NA_real_)
  z <- as.numeric(summary_table[hit, cols[[1L]]])
  if (!length(z) || !any(is.finite(z))) NA_real_ else min(z, na.rm = TRUE)
}

cga_sampler_diagnostics <- function(stanfit, max_treedepth = 13L) {
  if (is.null(stanfit) || !requireNamespace("rstan", quietly = TRUE))
    return(c(divergences = NA_real_, treedepth_hits = NA_real_, min_bfmi = NA_real_))
  sp <- tryCatch(rstan::get_sampler_params(stanfit, inc_warmup = FALSE),
                 error = function(e) list())
  if (!length(sp)) return(c(divergences = NA_real_, treedepth_hits = NA_real_,
                            min_bfmi = NA_real_))
  div <- sum(vapply(sp, function(x)
    if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"] > 0) else 0,
    numeric(1)))
  hits <- sum(vapply(sp, function(x)
    if ("treedepth__" %in% colnames(x)) sum(x[, "treedepth__"] >= max_treedepth) else 0,
    numeric(1)))
  bfmi <- vapply(sp, function(x) {
    if (!"energy__" %in% colnames(x)) return(NA_real_)
    e <- as.numeric(x[, "energy__"])
    if (length(e) < 2L || !is.finite(stats::var(e)) || stats::var(e) <= 0)
      return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  c(divergences = div, treedepth_hits = hits,
    min_bfmi = if (any(is.finite(bfmi))) min(bfmi, na.rm = TRUE) else NA_real_)
}

cga_predictive_errors <- function(stanfit, panel, fit = NULL, seed = 271828L) {
  n <- nrow(panel$data)
  share <- cga_draw_matrix(cga_extract(stanfit, "pred_shareIn"))
  if (!is.null(fit) && identical(fit@metadata$quality_parameterization, "collapsed")) {
    quality <- antitrustBayes::posterior_quality(fit, seed = seed)
    mu <- cga_draw_matrix(cga_extract(stanfit, "share_mean"))
    share <- exp(mu + quality$draws[, fit@stan_data$group, drop = FALSE])
  }
  margin_inv <- cga_draw_matrix(cga_extract(stanfit, "pred_marginNormInv"))
  obs_share <- cga_observed_shares(panel$data)
  out <- c(share_predictive_median_mae = NA_real_,
           share_predictive_median_max_abs = NA_real_,
           margin_predictive_median_mae = NA_real_,
           margin_predictive_median_max_abs = NA_real_)
  if (!is.null(share) && ncol(share) == n) {
    err <- abs(apply(share, 2L, stats::median, na.rm = TRUE) - obs_share)
    out[["share_predictive_median_mae"]] <- mean(err, na.rm = TRUE)
    out[["share_predictive_median_max_abs"]] <- max(err, na.rm = TRUE)
  }
  if (!is.null(margin_inv) && ncol(margin_inv) == n) {
    med <- apply(margin_inv, 2L, stats::median, na.rm = TRUE)
    err <- abs(1 / med - panel$data$margins)
    out[["margin_predictive_median_mae"]] <- mean(err, na.rm = TRUE)
    out[["margin_predictive_median_max_abs"]] <- max(err, na.rm = TRUE)
  }
  out
}

cga_fit_diagnostics <- function(fit, panel, case, elapsed,
                                fit_error = "", warnings = "") {
  base <- data.frame(
    case_id = case$case_id, design = case$design, mode = case$mode,
    demand_model = case$demand_model, supply_model = case$supply_model,
    structure = case$structure, seed = case$seed,
    share_scenario = case$share_scenario,
    status = if (grepl("^pending production API:", fit_error)) "pending_api"
      else if (nzchar(fit_error)) "fit_error" else "fit_complete",
    error = fit_error, warnings = warnings, elapsed_seconds = elapsed,
    gradient_evaluations = NA_real_, stringsAsFactors = FALSE
  )
  base$diagnostic_method <- "rank-normalized posterior package"
  base$gradient_evaluation_measure <- "sum warmup and sampling n_leapfrog__"
  base$stan_seconds <- NA_real_
  base$timing_note <- if (is.finite(elapsed)) "recorded fit wall time" else
    "fit wall time unavailable after postprocessing interruption; native Stan timing retained"
  base$predictive_quality_conditioning <- if (case$mode == "B")
    "conditional Gaussian quality draws given shares" else "sampled quality draws"
  for (name in c("alpha", "gamma", "s0", "tau", "sigma_share", "sigma_margin"))
    for (q in c("_q05", "_q50", "_q95")) base[[paste0(name, q)]] <- NA_real_
  d <- c(divergences = NA_real_, treedepth_hits = NA_real_, min_bfmi = NA_real_,
         max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
         alpha_ess_per_sec = NA_real_, gamma_ess_per_sec = NA_real_,
         alpha_or_gamma_ess_per_sec = NA_real_, s0_ess_per_sec = NA_real_,
         tau_ess_per_sec = NA_real_,
         sigma_b_strat_ess_per_sec = NA_real_, sigma_share_ess_per_sec = NA_real_,
         tau_sigma_correlation = NA_real_)
  p <- c(share_predictive_median_mae = NA_real_,
         share_predictive_median_max_abs = NA_real_,
         margin_predictive_median_mae = NA_real_,
         margin_predictive_median_max_abs = NA_real_)
  if (!is.null(fit) && !nzchar(fit_error)) {
    st <- cga_unwrap_stanfit(fit)
    if (!is.null(st) && requireNamespace("rstan", quietly = TRUE)) {
      d[names(cga_sampler_diagnostics(st, case$max_treedepth))] <-
        cga_sampler_diagnostics(st, case$max_treedepth)
      if (!requireNamespace("posterior", quietly = TRUE))
        stop("The frozen benchmark requires posterior for rank-normalized diagnostics.")
      raw <- rstan::extract(st, permuted = FALSE, inc_warmup = FALSE)
      monitored <- grepl("^(mu_log_a|mu_log_a_free|log_kappa|alpha_structural_market|gamma_structural_market|s0_prob|s0_logit|sigma_b_strat|b_group_centered|sigma_share|sigma_margin|log_sigma_share|log_sigma_margin|rho)",
                        dimnames(raw)[[3L]])
      tab <- as.data.frame(posterior::summarise_draws(
        posterior::as_draws_array(raw[, , monitored, drop = FALSE])))
      sm <- as.matrix(tab[, setdiff(names(tab), "variable"), drop = FALSE])
      rownames(sm) <- tab$variable
      colnames(sm)[colnames(sm) == "rhat"] <- "Rhat"
      if (!is.null(sm) && nrow(sm)) {
        hit <- grepl("^(mu_log_a|mu_log_a_free|log_kappa|alpha_structural_market|gamma_structural_market|s0_prob|s0_logit|sigma_b_strat|b_group_centered|sigma_share|sigma_margin|log_sigma_share|log_sigma_margin|rho)",
                     rownames(sm))
        rh <- if ("Rhat" %in% colnames(sm)) as.numeric(sm[hit, "Rhat"]) else numeric()
        d[["max_rhat"]] <- if (length(rh) && any(is.finite(rh))) max(rh, na.rm = TRUE) else NA_real_
        d[["min_bulk_ess"]] <- cga_first_ess(sm[hit, , drop = FALSE], ".*", "ess_bulk")
        d[["min_tail_ess"]] <- cga_first_ess(sm[hit, , drop = FALSE], ".*", "ess_tail")
        elapsed_safe <- max(as.numeric(elapsed), .Machine$double.eps)
        take_ess <- function(patterns) {
          if (!any(hit)) return(NA_real_)
          z <- cga_first_ess(sm[hit, , drop = FALSE], patterns, "ess_bulk")
          z / elapsed_safe
        }
        if (case$demand_model == "logit") {
          d[["alpha_ess_per_sec"]] <- take_ess(
            "alpha_structural_market")
          d[["alpha_or_gamma_ess_per_sec"]] <- d[["alpha_ess_per_sec"]]
        } else {
          d[["gamma_ess_per_sec"]] <- take_ess(
            "gamma_structural_market")
          d[["alpha_or_gamma_ess_per_sec"]] <- d[["gamma_ess_per_sec"]]
        }
        d[["s0_ess_per_sec"]] <- take_ess("s0_prob")
        d[["tau_ess_per_sec"]] <- take_ess("sigma_b_strat")
        d[["sigma_b_strat_ess_per_sec"]] <- take_ess("sigma_b_strat")
        d[["sigma_share_ess_per_sec"]] <- take_ess("sigma_share")
      }
      sb <- cga_extract(st, "sigma_b_strat")
      ss <- cga_extract(st, "sigma_share")
      if (!is.null(sb) && !is.null(ss) && length(sb) == length(ss) &&
          stats::sd(as.numeric(sb)) > 0 && stats::sd(as.numeric(ss)) > 0)
        d[["tau_sigma_correlation"]] <- stats::cor(as.numeric(sb), as.numeric(ss))
      pp <- cga_predictive_errors(st, panel, fit,
        seed = cga_derived_seed(case$seed, case$case_id, "quality_prediction"))
      p[names(pp)] <- pp
      sp <- rstan::get_sampler_params(st, inc_warmup = TRUE)
      base$gradient_evaluations <- sum(vapply(sp, function(x)
        sum(x[, "n_leapfrog__"]), numeric(1)))
      base$stan_seconds <- sum(rstan::get_elapsed_time(st))
      base$gradient_evaluation_measure <- "sum warmup and sampling n_leapfrog__"
      fields <- c(alpha = "alpha_structural_market", gamma = "gamma_structural_market",
        s0 = "s0_prob", tau = "sigma_b_strat", sigma_share = "sigma_share",
        sigma_margin = "sigma_margin")
      for (name in names(fields)) {
        values <- cga_extract(st, fields[[name]])
        qs <- if (is.null(values)) rep(NA_real_, 3L) else
          stats::quantile(as.numeric(values), c(.05, .5, .95), names = FALSE)
        for (q in seq_along(qs)) base[[paste0(name, c("_q05", "_q50", "_q95")[[q]])]] <- qs[[q]]
      }
    }
  }
  for (nm in names(d)) base[[nm]] <- unname(d[[nm]])
  for (name in c("alpha", "gamma", "s0", "tau", "sigma_share")) {
    base[[paste0(name, "_ess_per_stan_second")]] <-
      if (is.finite(elapsed) && is.finite(base$stan_seconds) && base$stan_seconds > 0)
        base[[paste0(name, "_ess_per_sec")]] * elapsed / base$stan_seconds else NA_real_
  }
  ## Native timings remain usable if the enclosing R timer was lost. Rebuild
  ## each ESS from the already computed normalized diagnostic table.
  if (!is.finite(elapsed) && exists("sm", inherits = FALSE) && base$stan_seconds > 0) {
    native_fields <- c(alpha = "alpha_structural_market", gamma = "gamma_structural_market",
      s0 = "s0_prob", tau = "sigma_b_strat", sigma_share = "sigma_share")
    for (name in names(native_fields)) base[[paste0(name, "_ess_per_stan_second")]] <-
      cga_first_ess(sm, native_fields[[name]], "ess_bulk") / base$stan_seconds
  }
  for (nm in names(p)) base[[nm]] <- unname(p[[nm]])
  gates <- c(base$divergences == 0, base$treedepth_hits == 0,
    base$min_bfmi >= .30, base$max_rhat <= 1.01,
    base$min_bulk_ess >= 100, base$min_tail_ess >= 100)
  base$usable <- length(gates) == 6L && all(!is.na(gates)) && all(gates)
  if (!nzchar(fit_error)) base$status <- if (base$usable) "usable" else "diagnostic_failure"
  base
}

cga_bridge_one <- function(fit) {
  if (is.null(fit) || !requireNamespace("bridgesampling", quietly = TRUE))
    return(list(status = "failed", log_marginal = NA_real_,
                error = "bridgesampling or fit unavailable"))
  st <- cga_unwrap_stanfit(fit)
  if (is.null(st)) return(list(status = "failed", log_marginal = NA_real_,
                               error = "unable to unwrap BayesFit@stanfit"))
  bridge_warnings <- character()
  value <- tryCatch(withCallingHandlers(
    bridgesampling::bridge_sampler(st, repetitions = 3L, maxiter = 1000L,
      method = "normal", cores = 1L, silent = TRUE),
    warning = function(w) {
      bridge_warnings <<- c(bridge_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(value, "condition"))
    return(list(status = "failed", log_marginal = NA_real_,
                error = conditionMessage(value)))
  logml <- if (!is.null(value$logml)) as.numeric(value$logml) else NA_real_
  iterations <- as.numeric(value$niter)
  converged <- length(logml) == 3L && all(is.finite(logml)) &&
    length(iterations) == 3L && all(is.finite(iterations)) && all(iterations < 1000L)
  list(status = if (converged && !length(bridge_warnings)) "complete" else "failed",
       log_marginal = if (all(is.finite(logml))) stats::median(logml) else NA_real_,
       error = paste(c(bridge_warnings,
         if (!converged) "bridge repetition/finite-value/convergence gate failed",
         paste0("repetitions=", length(logml), "; logml_range=",
           if (all(is.finite(logml))) diff(range(logml)) else NA_real_,
           "; niter=", paste(iterations, collapse = "/"))), collapse = " | "))
}

cga_pending_policy <- function(case_id, reason = "pending root policy integration") {
  data.frame(case_id = case_id, policy_status = "pending_root_integration",
             policy_error = reason,
             merger_share_change = NA_real_, merger_price_change = NA_real_,
             cost_share_change = NA_real_, combined_price_change = NA_real_,
             solver_failure_rate = NA_real_, stringsAsFactors = FALSE)
}

cga_run_one <- function(panel, case, out_dir, estimator = NULL,
                        reuse_fit = FALSE, recovered_elapsed = NA_real_,
                        recovered_bridge = NULL) {
  log_dir <- file.path(out_dir, "logs")
  fit_dir <- file.path(out_dir, "stanfit")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0(case$case_id, ".log"))
  fit_path <- file.path(fit_dir, paste0(case$case_id, ".rds"))
  timing_path <- file.path(fit_dir, paste0(case$case_id, ".processing.rds"))
  started <- proc.time()[["elapsed"]]
  fit <- NULL; fit_error <- ""; warns <- character(); output <- character()
  if (reuse_fit) {
    fit <- readRDS(fit_path)
    timing <- if (file.exists(timing_path)) readRDS(timing_path) else NULL
    elapsed <- if (!is.null(timing)) timing$elapsed else recovered_elapsed
    warns <- if (!is.null(timing)) timing$warnings else
      if (file.exists(log_path)) readLines(log_path, warn = FALSE)[startsWith(readLines(log_path, warn = FALSE), "WARNING:")] else character()
  } else {
  if (!cga_api_ready(estimator)) {
    fit_error <- paste("pending production API:",
                       "bayes_est() must expose share_likelihood and",
                       "quality_parameterization")
  } else {
    fn <- estimator
    if (is.null(fn)) fn <- getExportedValue("antitrustBayes", "bayes_est")
    output <- capture.output({
      result <- tryCatch({
        value <- withCallingHandlers(do.call(fn, cga_fit_args(panel, case)),
          warning = function(w) {
            warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
          })
        list(value = value, error = NULL)
      }, error = function(e) list(value = NULL, error = conditionMessage(e)))
      fit <- result$value
      if (!is.null(result$error)) fit_error <- result$error
    }, type = "output")
  }
  elapsed <- proc.time()[["elapsed"]] - started
  writeLines(c(output, if (length(warns)) paste("WARNING:", warns),
               if (nzchar(fit_error)) paste("ERROR:", fit_error)), log_path)
  if (!is.null(fit) && !nzchar(fit_error))
    saveRDS(fit, fit_path)
  saveRDS(list(elapsed = elapsed, warnings = warns, error = fit_error), timing_path)
  }
  diag <- cga_fit_diagnostics(fit, panel, case, elapsed, fit_error,
                              paste(warns, collapse = " | "))
  bridge <- if (!is.null(recovered_bridge)) recovered_bridge else
    if (isTRUE(case$bridge_requested) && !nzchar(fit_error))
    cga_bridge_one(fit) else list(status = "not_requested", log_marginal = NA_real_, error = "")
  list(diagnostics = diag, bridge = bridge,
       policy = cga_pending_policy(case$case_id), log_path = log_path,
       fit = fit)
}

cga_write <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

cga_run <- function(plan = cga_make_plan(), out_dir, estimator = NULL, resume = FALSE) {
  cga_validate_plan(plan)
  out_dir <- cga_require_tmp(out_dir, resume)
  previous_diag <- previous_bridge <- NULL
  if (resume) {
    previous_plan <- read.csv(file.path(out_dir, "benchmark_plan.csv"), stringsAsFactors = FALSE,
      na.strings = c("", "NA"))
    if (!isTRUE(all.equal(previous_plan, plan, check.attributes = FALSE)))
      stop("The saved plan does not equal the frozen resume plan.")
    if (file.exists(file.path(out_dir, "benchmark.csv")))
      previous_diag <- read.csv(file.path(out_dir, "benchmark.csv"), stringsAsFactors = FALSE)
    if (file.exists(file.path(out_dir, "benchmark_bridges.csv")))
      previous_bridge <- read.csv(file.path(out_dir, "benchmark_bridges.csv"), stringsAsFactors = FALSE)
    previous_state <- read.csv(file.path(out_dir, "benchmark_progress.csv"), stringsAsFactors = FALSE)
    uncertain <- previous_state$status != "unstarted" &
      !file.exists(file.path(out_dir, "stanfit", paste0(previous_state$case_id, ".rds")))
    if (any(uncertain)) stop("An interrupted fit has no saved posterior; preserve this run and use a separately labeled new attempt.")
  }
  ## Load common package code before fit timers begin.
  cga_api_ready(estimator)
  cga_write(plan, file.path(out_dir, "benchmark_plan.csv"))
  manifest <- data.frame(
    run_date = as.character(Sys.Date()), run_profile = "primary",
    chains = 4L, iter = 1000L, warmup = 500L, metric = "diag_e", adapt_delta = .95,
    max_treedepth = 13L, use_rho = FALSE, use_market_resid_scale = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, use_hmt = FALSE,
    claim = "no rescue, tuning, or policy claim from pending/failed fits",
    stringsAsFactors = FALSE
  )
  cga_write(manifest, file.path(out_dir, "benchmark_run_manifest.csv"))
  diag_rows <- list(); bridge_rows <- list(); policy_rows <- list(); ledger <- list()
  run_state <- data.frame(case_id = plan$case_id, status = "unstarted",
    artifact = "", stringsAsFactors = FALSE)
  checkpoint <- function() {
    cga_write(run_state, file.path(out_dir, "benchmark_progress.csv"))
    if (length(diag_rows)) cga_write(do.call(rbind, diag_rows), file.path(out_dir, "benchmark.csv"))
    if (length(bridge_rows)) cga_write(do.call(rbind, bridge_rows), file.path(out_dir, "benchmark_bridges.csv"))
    if (length(policy_rows)) cga_write(do.call(rbind, policy_rows), file.path(out_dir, "benchmark_policy_hooks.csv"))
    if (length(ledger)) cga_write(do.call(rbind, ledger), file.path(out_dir, "benchmark_run_ledger.csv"))
  }
  checkpoint()
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    run_state$status[[i]] <- "running"
    checkpoint()
    cat("Starting ", case$case_id, "\n", sep = "")
    flush.console()
    panel <- tryCatch(cga_make_dgp(seed = case$seed,
      demand_model = case$demand_model, supply_model = case$supply_model,
      structure = case$structure, share_scenario = case$share_scenario),
      error = function(e) e)
    if (inherits(panel, "condition")) {
      err <- paste("DGP error:", conditionMessage(panel))
      diag_rows[[length(diag_rows) + 1L]] <- cga_fit_diagnostics(
        NULL, list(data = data.frame()), case, 0, err)
      bridge_rows[[length(bridge_rows) + 1L]] <- data.frame(
        case_id = case$case_id, bridge_status = "not_requested",
        bridge_log_marginal = NA_real_, bridge_error = err, stringsAsFactors = FALSE)
      policy_rows[[length(policy_rows) + 1L]] <- cga_pending_policy(case$case_id, err)
      ledger[[length(ledger) + 1L]] <- data.frame(case_id = case$case_id,
        date = as.character(Sys.Date()), stage = "dgp", status = "dgp_error",
        evidence = err, artifact = "", stringsAsFactors = FALSE)
      run_state$status[[i]] <- "dgp_error"
      checkpoint()
      next
    }
    valid <- tryCatch(cga_validate_dgp(panel), error = function(e) e)
    if (inherits(valid, "condition")) {
      err <- paste("DGP validation failed:", conditionMessage(valid))
      diag_rows[[length(diag_rows) + 1L]] <- cga_fit_diagnostics(
        NULL, panel, case, 0, err)
      bridge_rows[[length(bridge_rows) + 1L]] <- data.frame(
        case_id = case$case_id, bridge_status = "not_requested",
        bridge_log_marginal = NA_real_, bridge_error = err, stringsAsFactors = FALSE)
      policy_rows[[length(policy_rows) + 1L]] <- cga_pending_policy(case$case_id, err)
      ledger[[length(ledger) + 1L]] <- data.frame(case_id = case$case_id,
        date = as.character(Sys.Date()), stage = "dgp", status = "dgp_error",
        evidence = err, artifact = "", stringsAsFactors = FALSE)
      run_state$status[[i]] <- "dgp_error"
      checkpoint()
      next
    }
    reuse <- resume && file.exists(file.path(out_dir, "stanfit", paste0(case$case_id, ".rds")))
    old <- if (!is.null(previous_diag)) previous_diag[previous_diag$case_id == case$case_id, , drop = FALSE] else NULL
    old_bridge <- if (!is.null(previous_bridge)) previous_bridge[previous_bridge$case_id == case$case_id, , drop = FALSE] else NULL
    recovered_bridge <- if (!is.null(old_bridge) && nrow(old_bridge) == 1L) list(
      status = old_bridge$bridge_status[[1]], log_marginal = old_bridge$bridge_log_marginal[[1]], error = old_bridge$bridge_error[[1]]) else NULL
    if (reuse) cat("Reusing saved posterior; no HMC rerun: ", case$case_id, "\n", sep = "")
    result <- cga_run_one(panel, case, out_dir, estimator, reuse_fit = reuse,
      recovered_elapsed = if (!is.null(old) && nrow(old) == 1L) old$elapsed_seconds[[1]] else NA_real_,
      recovered_bridge = recovered_bridge)
    diag_rows[[length(diag_rows) + 1L]] <- result$diagnostics
    bridge_rows[[length(bridge_rows) + 1L]] <- data.frame(
      case_id = case$case_id, bridge_status = result$bridge$status,
      bridge_log_marginal = result$bridge$log_marginal,
      bridge_error = result$bridge$error, stringsAsFactors = FALSE)
    policy_rows[[length(policy_rows) + 1L]] <- result$policy
    ledger[[length(ledger) + 1L]] <- data.frame(case_id = case$case_id,
      date = as.character(Sys.Date()), stage = "fit", status = result$diagnostics$status,
      evidence = cga_or_null(result$diagnostics$error[[1L]], ""),
      artifact = if (file.exists(result$log_path)) result$log_path else "",
      stringsAsFactors = FALSE)
    run_state$status[[i]] <- result$diagnostics$status[[1L]]
    run_state$artifact[[i]] <- result$log_path
    checkpoint()
    cat("Finished ", case$case_id, ": ", result$diagnostics$status[[1L]], "\n", sep = "")
    flush.console()
  }
  diagnostics <- do.call(rbind, diag_rows)
  bridges <- do.call(rbind, bridge_rows)
  policies <- do.call(rbind, policy_rows)
  ledger_df <- do.call(rbind, ledger)
  cga_write(diagnostics, file.path(out_dir, "benchmark.csv"))
  cga_write(bridges, file.path(out_dir, "benchmark_bridges.csv"))
  cga_write(policies, file.path(out_dir, "benchmark_policy_hooks.csv"))
  cga_write(ledger_df, file.path(out_dir, "benchmark_run_ledger.csv"))
  invisible(list(plan = plan, diagnostics = diagnostics, bridges = bridges,
                 policies = policies, ledger = ledger_df, output = out_dir))
}

cga_dry_run <- function(plan = cga_make_plan(), verbose = TRUE) {
  cga_validate_plan(plan)
  dgp_checks <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    panel <- cga_make_dgp(seed = case$seed, demand_model = case$demand_model,
      supply_model = case$supply_model, structure = case$structure,
      share_scenario = case$share_scenario)
    check <- cga_validate_dgp(panel)
    dgp_checks[[i]] <- data.frame(case_id = case$case_id,
      n_observations = check$n_observations,
      market_years = check$n_market_year,
      max_stan_mean_error = check$max_stan_mean_error,
      max_quality_center_error = check$max_quality_center_error,
      max_share_sum_error = check$max_share_sum_error,
      max_foc_error = check$max_foc_error, status = check$status,
      stringsAsFactors = FALSE)
  }
  fixtures <- cga_ladder_fixtures()
  fixture_checks <- lapply(fixtures, function(x) cga_validate_dgp(x))
  ## Exercise the same mode-B argument assembly at the two higher rungs.  No
  ## sampler is called here; this catches accidental one-period/scalar API
  ## assumptions before the production estimator is used.
  ladder_case <- plan[plan$mode == "B" &
    plan$design == "logit_bertrand_noisy_271828", , drop = FALSE][1L, , drop = FALSE]
  ladder_args <- lapply(fixtures, function(x) cga_fit_args(x, ladder_case))
  stopifnot(identical(ladder_case$quality_parameterization, "collapsed"),
            length(ladder_args$four_period$outside_price) == 4L,
            length(ladder_args$multi_market$outside_price) == 4L,
            length(unique(fixtures$multi_market$data$market)) == 2L,
            length(unique(fixtures$multi_market$data$year)) == 2L,
            length(unique(fixtures$multi_market$data$group)) == 48L)
  if (verbose) {
    print(plan[, c("case_id", "mode", "demand_model", "supply_model",
                   "structure", "share_scenario", "k_support", "k_prior")],
          row.names = FALSE)
    cat("DGP VALIDATION: 21/21 primary cells passed; four-period and",
        "2-market x 2-period fixtures passed.\n")
    cat("DRY RUN ONLY: no HMC, fit object, bridge, policy, or recovery claim.\n")
  }
  invisible(list(plan = plan, dgp = do.call(rbind, dgp_checks),
                 fixtures = fixture_checks, fixture_args = ladder_args))
}

if (sys.nframe() == 0L) {
  cli <- cga_parse_cli()
  if (isTRUE(cli$help)) {
    cga_usage()
    quit(save = "no", status = 0L)
  }
  plan <- cga_make_plan()
  cga_validate_plan(plan)
  if (!isTRUE(cli$run)) {
    cga_dry_run(plan)
    quit(save = "no", status = 0L)
  }
  if (is.null(cli$out)) stop("--run requires --out=/tmp/...")
  result <- cga_run(plan, cli$out, resume = cli$resume)
  cat("PRIMARY RUN COMPLETE: compact benchmark CSV, bridge hook results,",
      "policy pending hooks, and ledger written under ",
      normalizePath(result$output), ".\n", sep = "")
}
