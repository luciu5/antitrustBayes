#!/usr/bin/env Rscript
## Independent native target/gradient gate for the market-s0 rotation.
## Run once against the original package (ROT_EXPECT_ROTATED=0) and once
## against the fresh experimental package (ROT_EXPECT_ROTATED=1). No HMC.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("ROT_LIB", unset = "")
if (!nzchar(lib)) stop("ROT_LIB is required", call. = FALSE)
is_rot <- identical(Sys.getenv("ROT_EXPECT_ROTATED", "0"), "1")
chain_mode <- is_rot && identical(Sys.getenv("ROT_CHAIN", "0"), "1")
single_mode <- identical(Sys.getenv("ROT_SINGLE", "0"), "1")
input_mode <- identical(Sys.getenv("ROT_INPUT", "0"), "1")
sign_only <- identical(Sys.getenv("ROT_SIGN_ONLY", "0"), "1")
endpoint <- Sys.getenv("ROT_ENDPOINT", "")
student_mode <- identical(Sys.getenv("ROT_STUDENT", "0"), "1")
centered_alpha_mode <- identical(Sys.getenv("ROT_CENTERED_ALPHA", "0"), "1")
out <- Sys.getenv("ROT_OUTPUT", unset = "/tmp/statistical-redteam-20260916/rotation-native-results.rds")
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = lib))
suppressPackageStartupMessages(library(rstan))

## Four products per market-year, with repeated ownership to exercise the
## quality loading. Two years and two markets distinguish market-year and
## market pooling. A few margins are missing; the remaining margins support
## rho-active cases.
make_fixture <- function() {
  M <- if (single_mode) 1L else 2L
  T <- as.integer(Sys.getenv("ROT_T", "2")); if (T < 1L) stop("ROT_T must be positive")
  J <- 3L; L <- M * T
  list(
    prices = rep(c(1.05, 1.28, 1.52), L) + rep(c(0, .01, -.01), L),
    quantities = rep(c(32, 21, 14), L) + rep(c(0, 1, -1), L),
    margins = rep(c(.30, .27, .33), L),
    market = rep(rep(1:M, each = T), each = J),
    year = rep(rep(1:T, times = M), each = J),
    group = rep(c(1L, 1L, 2L), L),
    margin_observed = rep(c(TRUE, FALSE, TRUE), L),
    outside_price = seq(.85, 1.35, length.out = L),
    loo_cluster = rep(rep(1:M, each = T), each = J)
  )
}
fx <- make_fixture()

set_all <- function(init, name, value) {
  if (!name %in% names(init)) return(init)
  old <- init[[name]]
  val <- if (length(value) == length(old)) as.numeric(value) else
    rep(as.numeric(value), length(old))
  if (!is.null(dim(old))) dim(val) <- dim(old)
  init[[name]] <- val
  init
}
init_one <- function(fit) {
  z <- rstan::get_inits(fit)
  if (is.list(z) && length(z) == 1L && is.list(z[[1L]]) &&
      !is.null(names(z[[1L]]))) z <- z[[1L]]
  if (!is.list(z) || is.null(names(z))) stop("no constrained init", call. = FALSE)
  z[intersect(fit@model_pars, names(z))]
}

model_for <- function(demand, structure) {
  nm <- paste0(if (demand == "ces") "ces_" else "logit_",
               if (structure == "simultaneous") "supply" else "stackelberg")
  antitrustBayes:::stanmodels[[nm]]
}

prepare_case <- function(demand, supply, structure, s0mode, centered, quality,
                         alpha_pool, rho, exact, rotation) {
  support <- if (!exact) NULL else if (structure == "stackelberg") 1:2 else 0:2
  prior <- if (is.null(support)) NULL else if (structure == "stackelberg")
    c(.43, .57) else c(.18, .32, .50)
  args <- list(prices = fx$prices, quantities = fx$quantities,
    margins = fx$margins, market = fx$market, year = fx$year, group = fx$group,
    supply_model = if (supply == "bertrand") 1L else 3L,
    use_cutoff = FALSE, market_side = if (input_mode) "input" else "output", use_rho = rho,
    use_market_resid_scale = FALSE, use_hmt = FALSE,
    use_alpha_market_pool = alpha_pool, use_centered_alpha = centered_alpha_mode,
    use_centered_s0_pool = centered,
    use_s0_market_year_pool = identical(s0mode, "market_year"),
    use_s0_market_pool = identical(s0mode, "market"),
    priors = list(),
    X_s0 = cbind(seq(-.16, .16, length.out = length(fx$outside_price)),
                 rep(c(-.25, .25), length.out = length(fx$outside_price))),
    margin_observed = fx$margin_observed, outside_price = fx$outside_price,
    likelihood_mode = 1L, loo_cluster = fx$loo_cluster,
    estimate_df = FALSE, resid_df_fixed = 5, structure = structure,
    partition_prior = if (exact) "ordered_k" else "cutoff", k_support = support,
    k_prior = prior, k_scope = "market",
    quality_parameterization = quality,
    margin_likelihood = if (student_mode) "student_t" else "gaussian",
    share_likelihood = if (student_mode) "student_t" else "gaussian")
  if (is_rot) args$market_s0_parameterization <- rotation
  fn <- if (demand == "ces") antitrustBayes:::prepare_ces_stan_data
        else antitrustBayes:::prepare_stan_data
  do.call(fn, args)
}

set_common <- function(init, demand, supply, quality, alpha_pool, centered,
                       s0mode, rotation, exact) {
  M <- if (single_mode) 1L else 2L
  T <- max(fx$year)
  b <- if (M == 1L) 0 else c(.22, -.17)
  eta <- seq(-1.16, -.68, length.out = T * M)
  mu <- -1.10; sig <- .37
  init <- set_all(init, "mu_b_market_vec", -.04)
  init <- set_all(init, "sigma_b_market_vec", .43)
  init <- set_all(init, "b_market_raw", (b + .04) / .43)
  init <- set_all(init, "sigma_b_strat", .28)
  init <- set_all(init, "b_group_centered", .11)
  init <- set_all(init, "beta_s0", c(.12, -.08))
  init <- set_all(init, "log_sigma_share", log(.35))
  init <- set_all(init, "log_sigma_margin", log(.22))
  init <- set_all(init, "rho_raw", atanh(.31 / .95))
  market_of_mky <- rep(seq_len(M), each = T)
  rotation_active <- rotation == "rotated" && M > 1L
  if (s0mode == "none") {
    base <- eta
    if (rotation_active) base <- eta + b[market_of_mky]
    init <- set_all(init, if (rotation_active) "c_s0" else "s0_logit", base)
  } else if (s0mode == "market_year") {
    z <- (eta - mu) / sig
    if (rotation_active && !centered)
      init <- set_all(init, "c_s0_raw", z + b[market_of_mky] / sig)
    else if (rotation_active && centered)
      init <- set_all(init, "c_s0", eta + b[market_of_mky])
    else if (centered) init <- set_all(init, "eta_s0_pool", eta)
    else {
      init <- set_all(init, "mu_s0_pool_vec", mu)
      init <- set_all(init, "sigma_s0_pool_vec", sig)
      init <- set_all(init, "s0_pool_raw", z)
    }
  } else {
    em <- vapply(seq_len(M), function(m) {
      idx <- ((m - 1L) * T + 1L):(m * T)
      mean(eta[idx])
    }, numeric(1))
    z <- (em - mu) / sig
    if (rotation_active && !centered)
      init <- set_all(init, "c_s0_raw", z + b / sig)
    else if (rotation_active && centered)
      init <- set_all(init, "c_s0", em + b)
    else if (centered) init <- set_all(init, "eta_s0_market", em)
    else {
      init <- set_all(init, "mu_s0_pool_vec", mu)
      init <- set_all(init, "sigma_s0_pool_vec", sig)
      init <- set_all(init, "s0_market_raw", z)
    }
  }
  init <- set_all(init, "mu_s0_pool_vec", mu)
  init <- set_all(init, "sigma_s0_pool_vec", sig)
  if (demand == "ces") init <- set_all(init, "mu_log_a", log(.82))
  else if (supply == "cournot" && !(alpha_pool && M > 1L))
    init <- set_all(init, "log_kappa_vec", log(.63))
  else if (supply == "cournot" && alpha_pool && M > 1L)
    init <- set_all(init, "log_kappa_market", log(c(.60, .69)))
  else if (!(alpha_pool && M > 1L))
    init <- set_all(init, "mu_log_a_free", log(.82))
  else {
    init <- set_all(init, "sigma_log_a_vec", .25)
    init <- set_all(init, "r_market_a_raw", c(-.2, .3))
    init <- set_all(init, "log_a_market_pool", log(c(.72, .91)))
  }
  init
}

run_case <- function(demand, supply, structure, s0mode, centered, quality,
                     alpha_pool, rho, exact, rotation) {
  dat <- prepare_case(demand, supply, structure, s0mode, centered, quality,
                      alpha_pool, rho, exact, rotation)
  model <- model_for(demand, structure)
  fit0 <- rstan::sampling(model, data = dat, chains = 1L, iter = 1L,
                          warmup = 0L, algorithm = "Fixed_param", init = 0,
                          seed = 19001L, refresh = 0)
  init <- set_common(init_one(fit0), demand, supply, quality, alpha_pool,
                     centered, s0mode, rotation, exact)
  u <- rstan::unconstrain_pars(fit0, init)
  lp <- as.numeric(rstan::log_prob(fit0, u, adjust_transform = TRUE))
  ad <- as.numeric(rstan::grad_log_prob(fit0, u, adjust_transform = TRUE))
  h <- 1e-5
  fd <- vapply(seq_along(u), function(j) {
    up <- u; um <- u
    up[[j]] <- up[[j]] + h
    um[[j]] <- um[[j]] - h
    (as.numeric(rstan::log_prob(fit0, up, adjust_transform = TRUE)) -
       as.numeric(rstan::log_prob(fit0, um, adjust_transform = TRUE))) / (2 * h)
  }, numeric(1))
  fd_rel <- max(abs(ad - fd) / (1 + abs(ad)))
  ## One Fixed_param draw at the checked constrained state exposes economic
  ## transformed/generated outputs without running HMC.
  at_state <- rstan::sampling(model, data = dat, chains = 1L, iter = 1L,
                              warmup = 0L, algorithm = "Fixed_param",
                              init = list(init), seed = 19002L, refresh = 0)
  out_pars <- intersect(c("s0_prob", "share_mean", "pred_marginNormInv",
                          "alpha_structural_market", "gamma_market",
                          "R0_prob", "c_effective_mky",
                          "posterior_component_prob", "prior_component_prob"),
                        at_state@model_pars)
  gq <- if (length(out_pars)) {
    got <- rstan::extract(at_state, pars = out_pars, permuted = TRUE)
    lapply(got, as.numeric)
  } else list()
  list(lp = lp, grad = ad, fd = fd, fd_rel = fd_rel, u = u,
       init = init, fit = fit0, data = dat, gq = gq)
}

## Map an experimental rotated constrained state into the corresponding
## legacy constrained state. The finite-difference Jacobian below includes
## sigma and b cross derivatives of the shifted noncentered coordinate.
rotated_to_legacy_u <- function(rotated, legacy_fit, rotated_fit, s0mode,
                                centered) {
  newc <- rstan::constrain_pars(rotated_fit, rotated)
  oldc <- init_one(legacy_fit)
  for (nm in intersect(names(newc), names(oldc))) {
    if (length(newc[[nm]]) == length(oldc[[nm]])) oldc[[nm]] <- newc[[nm]]
  }
  if (!is.null(newc$b_market_raw) && length(newc$b_market_raw)) {
    b <- as.numeric(newc$mu_b_market_vec) +
      as.numeric(newc$sigma_b_market_vec) * as.numeric(newc$b_market_raw)
  } else b <- 0
  sig <- if (!is.null(newc$sigma_s0_pool_vec) &&
             length(newc$sigma_s0_pool_vec))
    as.numeric(newc$sigma_s0_pool_vec) else NA_real_
  M <- if (length(b)) length(b) else 1L
  mky_market <- rep(seq_len(M), each = max(fx$year))
  if (s0mode == "none") {
    if ("c_s0" %in% names(newc))
      oldc$s0_logit <- as.numeric(newc$c_s0) - b[mky_market]
  } else if (s0mode == "market_year") {
    if (centered && "c_s0" %in% names(newc))
      oldc$eta_s0_pool <- as.numeric(newc$c_s0) - b[mky_market]
    if (!centered && "c_s0_raw" %in% names(newc))
      oldc$s0_pool_raw <- as.numeric(newc$c_s0_raw) - b[mky_market] / sig
  } else {
    if (centered && "c_s0" %in% names(newc))
      oldc$eta_s0_market <- as.numeric(newc$c_s0) - b
    if (!centered && "c_s0_raw" %in% names(newc))
      oldc$s0_market_raw <- as.numeric(newc$c_s0_raw) - b / sig
  }
  as.numeric(rstan::unconstrain_pars(legacy_fit, oldc))
}

chain_check <- function(legacy, rotated, s0mode, centered) {
  ## In the single-market branch rotation is deliberately disabled and the
  ## native parameter vector is identical.  Avoid probing zero-width
  ## declarations through rstan's finite-difference mapper; the direct
  ## target/GQ and AD/FD checks still run for every single-market case.
  if (single_mode)
    return(max(abs(rotated$grad - legacy$grad) /
                 (1 + abs(legacy$grad))))
  f <- function(v) rotated_to_legacy_u(v, legacy$fit, rotated$fit,
                                       s0mode, centered)
  u <- rotated$u
  base <- f(u)
  if (length(base) != length(u)) stop("legacy/rotated dimensions differ")
  h <- 1e-6
  J <- matrix(0, nrow = length(base), ncol = length(u))
  for (j in seq_along(u)) {
    up <- u; um <- u
    up[[j]] <- up[[j]] + h
    um[[j]] <- um[[j]] - h
    J[, j] <- (f(up) - f(um)) / (2 * h)
  }
  expected <- as.numeric(crossprod(J, legacy$grad))
  max(abs(rotated$grad - expected) / (1 + abs(expected)))
}

## Balanced coverage: all four kernels, every s0 mode and centered state,
## sampled/collapsed quality, common/pooled alpha, rho off/on, and exact-k
## fringe/Stackelberg states.
spec <- expand.grid(
  demand = c("logit", "ces"), supply = c("bertrand", "cournot"),
  structure = c("simultaneous", "fringe", "stackelberg"),
  s0mode = c("none", "market_year", "market"), centered = c(FALSE, TRUE),
  quality = c("sampled", "collapsed"), alpha_pool = c(FALSE, TRUE),
  rho = c(FALSE, TRUE), exact = c(FALSE, TRUE), stringsAsFactors = FALSE
)
spec <- spec[spec$exact == (spec$structure != "simultaneous"), , drop = FALSE]
## The declared collapsed-Gaussian rho branch is intentionally supported only
## for one market; this two-market gate records the supported rho-active
## sampled cases and excludes that explicit unsupported combination.
spec <- spec[!(spec$quality == "collapsed" & spec$rho & !single_mode), , drop = FALSE]
## Keep every requested axis in non-exact cells. Exact cells use both quality
## families and both alpha/rho states while avoiding a large redundant grid.
keep <- with(spec, !exact | (exact & s0mode %in% c("market_year", "market") &
  ((quality == "sampled" & !alpha_pool & !rho) |
   (quality == "collapsed" & alpha_pool & rho))))
spec <- spec[keep, , drop = FALSE]
if (sign_only)
  spec <- spec[!spec$exact & spec$demand == "logit" &
    spec$structure == "simultaneous" & spec$s0mode == "market" &
    !spec$centered, , drop = FALSE]
if (endpoint %in% c("oneM1T", "multiM1T"))
  spec <- spec[!spec$exact & spec$s0mode == "market" & !spec$centered &
    spec$quality == "sampled" & !spec$rho &
    spec$structure %in% c("simultaneous", "stackelberg"), , drop = FALSE]
if (student_mode)
  spec <- spec[!spec$exact & spec$demand == "logit" &
    spec$structure == "simultaneous" &
    spec$s0mode %in% c("none", "market") & !spec$centered &
    spec$quality == "sampled" & !spec$alpha_pool & !spec$rho, , drop = FALSE]
if (centered_alpha_mode)
  spec <- spec[!spec$exact & spec$structure == "simultaneous" &
    spec$s0mode == "market" & !spec$centered & spec$quality == "sampled" &
    spec$alpha_pool & !spec$rho, , drop = FALSE]

rows <- list()
details <- list()
for (i in seq_len(nrow(spec))) {
  z <- spec[i, ]
  key <- paste(z$demand, z$supply, z$structure, z$s0mode,
               ifelse(z$centered, "C", "NC"), z$quality,
               ifelse(z$alpha_pool, "pool", "common"),
               ifelse(z$rho, "rho", "norho"),
               ifelse(z$exact, "exact", "full"), sep = "|")
  if (chain_mode) {
    pair <- tryCatch(list(
      legacy = run_case(z$demand, z$supply, z$structure, z$s0mode,
                        z$centered, z$quality, z$alpha_pool, z$rho,
                        z$exact, "legacy"),
      rotated = run_case(z$demand, z$supply, z$structure, z$s0mode,
                         z$centered, z$quality, z$alpha_pool, z$rho,
                         z$exact, "rotated")),
      error = function(e) structure(list(message = conditionMessage(e)),
                                    class = "rotation_error"))
    if (inherits(pair, "rotation_error")) ans <- pair
    else {
      ans <- pair$rotated
      ans$chain_rel <- tryCatch(
        chain_check(pair$legacy, pair$rotated, z$s0mode, z$centered),
        error = function(e) Inf)
    }
  } else {
    ans <- tryCatch(run_case(z$demand, z$supply, z$structure, z$s0mode,
                             z$centered, z$quality, z$alpha_pool, z$rho,
                             z$exact, if (is_rot) "rotated" else "legacy"),
                    error = function(e) structure(list(message = conditionMessage(e)),
                                                    class = "rotation_error"))
  }
  if (inherits(ans, "rotation_error"))
    rows[[length(rows) + 1L]] <- data.frame(key = key, error = ans$message,
                                             lp = NA_real_, fd_rel = NA_real_,
                                             chain_rel = NA_real_,
                                             stringsAsFactors = FALSE)
  else
    rows[[length(rows) + 1L]] <- data.frame(key = key, error = "",
                                             lp = ans$lp, fd_rel = ans$fd_rel,
                                             chain_rel = if (chain_mode) ans$chain_rel else NA_real_,
                                             stringsAsFactors = FALSE)
  if (!inherits(ans, "rotation_error")) {
    details[[key]] <- list(lp = ans$lp, grad = ans$grad, fd = ans$fd,
                           fd_rel = ans$fd_rel, chain_rel = if (chain_mode) ans$chain_rel else NULL,
                           init = ans$init, u = ans$u, gq = ans$gq)
  }
  cat(if (inherits(ans, "rotation_error")) "FAIL " else "OK ", key, "\n", sep = "")
}
saveRDS(list(spec = spec, rows = do.call(rbind, rows),
             details = details, is_rotated = is_rot), out)
