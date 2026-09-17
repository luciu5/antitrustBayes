suppressPackageStartupMessages({
  library(antitrustBayes)
  library(rstan)
})

role <- Sys.getenv("HELMERT_ROLE", "candidate")
out_file <- Sys.getenv("HELMERT_OUT", paste0("/tmp/helmert_", role, ".rds"))
raw_file <- Sys.getenv("HELMERT_RAW", "/tmp/helmert_candidate_raw.rds")
ns <- asNamespace("antitrustBayes")
mods <- get("stanmodels", ns)

## The fixture/helper is deliberately outside either source tree.  It exercises
## the public preparation functions and is shared by both installed versions.
invisible(capture.output(sys.source("doc/performance_audit/native_profile/helmert_panel_fixtures.R", envir = environment())))

kernel_for <- function(demand, structure) {
  if (structure == "simultaneous") paste0(demand, "_supply")
  else paste0(demand, "_stackelberg")
}

param_names <- function(sf) {
  p <- sf@sim$pars_oi
  if (is.null(p)) character() else as.character(p)
}

set_nonzero <- function(cp, sf, nm, value) {
  pnames <- param_names(sf)
  if (!(nm %in% pnames) || is.null(cp[[nm]]) || !length(cp[[nm]])) return(cp)
  cp[[nm]][] <- rep(value, length.out = length(cp[[nm]]))
  cp
}

nonzero_pars <- function(cp, sf) {
  ## Keep all values in a benign interior region.  In particular rho_raw and
  ## bargaining_unit are unconstrained through bounded transforms in Stan.
  for (nm in c("mu_b_market_vec")) cp <- set_nonzero(cp, sf, nm, .15)
  for (nm in c("sigma_b_market_vec", "sigma_b_strat", "sigma_log_a_vec",
               "sigma_s0_pool_vec", "sigma_log_sigma_share_market",
               "sigma_log_sigma_margin_market")) cp <- set_nonzero(cp, sf, nm, .7)
  for (nm in c("b_market_raw", "r_market_a_raw", "s0_logit", "s0_pool_raw",
               "eta_s0_pool", "s0_market_raw", "eta_s0_market", "c_s0",
               "c_s0_raw", "mu_s0_pool_vec", "mu_log_a_free",
               "log_a_market_pool", "log_kappa_vec", "log_kappa_market",
               "log_sigma_share_market_raw", "log_sigma_margin_market_raw",
               "beta_s0")) cp <- set_nonzero(cp, sf, nm, .08)
  for (nm in c("log_sigma_share", "log_sigma_margin"))
    cp <- set_nonzero(cp, sf, nm, -.15)
  cp <- set_nonzero(cp, sf, "rho_raw", .18)
  cp <- set_nonzero(cp, sf, "cutoff_share", .35)
  cp <- set_nonzero(cp, sf, "bargaining_unit", .4)
  cp <- set_nonzero(cp, sf, "nu_minus_2", 3)
  if ("b_group_centered" %in% param_names(sf) &&
      length(cp$b_group_centered))
    cp$b_group_centered[] <- seq(-.3, .3, length.out = length(cp$b_group_centered))
  cp
}

flat_numeric <- function(x) {
  y <- unlist(x, recursive = TRUE, use.names = TRUE)
  y <- y[vapply(y, is.numeric, logical(1))]
  as.numeric(y) |> setNames(names(y))
}

flat_q <- function(q) {
  ## Keep all numeric transformed/generated values.  Names retain variable
  ## and element positions so this catches TP/GQ shape and ordering changes.
  y <- unlist(q, recursive = TRUE, use.names = TRUE)
  y <- y[vapply(y, is.numeric, logical(1))]
  as.numeric(y) |> setNames(names(y))
}

representative_fd <- function(sf, u, lp_fun) {
  n <- length(u)
  ix <- sort(unique(c(1L, ceiling(n / 2), n)))
  h <- 1e-5
  val <- vapply(ix, function(i) {
    up <- u; um <- u
    up[i] <- up[i] + h; um[i] <- um[i] - h
    (lp_fun(up) - lp_fun(um)) / (2 * h)
  }, numeric(1))
  list(index = ix, value = val)
}

one <- function(case_name, spec, supplied_u = NULL) {
  started <- proc.time()[["elapsed"]]
  ans <- tryCatch({
    dat <- do.call(make_data, spec)
    kernel <- kernel_for(spec$demand, spec$structure)
    model <- mods[[kernel]]
    if (is.null(model)) stop("missing stan model ", kernel)
    sf <- suppressWarnings(sampling(model, data = dat, chains = 1L,
      iter = 1L, warmup = 0L, algorithm = "Fixed_param", init = 0,
      refresh = 0L, seed = 48003L))
    nu <- get_num_upars(sf)
    schema <- list(pars_oi = param_names(sf), dims_oi = sf@sim$dims_oi,
                   num_upars = nu)
    if (is.null(supplied_u)) {
      cp <- constrain_pars(sf, rep(0, nu))
      cp <- nonzero_pars(cp, sf)
      u <- unconstrain_pars(sf, cp)
    } else {
      u <- as.numeric(supplied_u)
      if (length(u) != nu)
        stop("raw parameter dimension differs: candidate=", length(u), " current=", nu)
    }
    lp_fun <- function(z, jac = TRUE)
      log_prob(sf, z, adjust_transform = jac, gradient = FALSE)
    lp <- lp_fun(u, TRUE)
    lp_raw <- lp_fun(u, FALSE)
    gr <- grad_log_prob(sf, u, adjust_transform = TRUE)
    gr_raw <- grad_log_prob(sf, u, adjust_transform = FALSE)
    fd <- representative_fd(sf, u, function(z) lp_fun(z, TRUE))
    cp_at <- constrain_pars(sf, u)
    tp <- flat_numeric(cp_at)
    sf_q <- suppressWarnings(sampling(model, data = dat, chains = 1L,
      iter = 1L, warmup = 0L, algorithm = "Fixed_param",
      init = function(...) cp_at, refresh = 0L, seed = 48005L))
    q <- flat_q(extract(sf_q, permuted = TRUE, inc_warmup = FALSE))
    list(status = "ok", case = case_name, kernel = kernel, demand = spec$demand,
      structure = spec$structure, supply_model = spec$supply_model,
      N = dat$N, N_market = dat$N_market, N_market_year = dat$N_market_year,
      N_group = dat$N_group, K_s0 = dat$K_s0, schema = schema, u = u,
      lp = as.numeric(lp), lp_raw = as.numeric(lp_raw), grad = as.numeric(gr),
      grad_raw = as.numeric(gr_raw), fd = fd, tp = tp, q = q,
      elapsed = proc.time()[["elapsed"]] - started)
  }, error = function(e) {
    list(status = "error", case = case_name, error = conditionMessage(e),
      elapsed = proc.time()[["elapsed"]] - started)
  })
  ans
}

res <- list()
raw <- if (role == "old") readRDS(raw_file) else NULL
for (nm in names(cases)) {
  cat(role, nm, "\n", sep = " ")
  res[[nm]] <- one(nm, cases[[nm]], if (role == "old") raw$u[[nm]] else NULL)
}
if (role == "candidate")
  saveRDS(list(u = lapply(res, function(x) if (x$status == "ok") x$u else NULL),
               cases = names(cases)), raw_file)
saveRDS(list(role = role, results = res, cases = names(cases)), out_file)

rows <- lapply(res, function(x) {
  if (x$status != "ok") return(data.frame(case = x$case, status = x$status,
    error = x$error, stringsAsFactors = FALSE))
  data.frame(case = x$case, status = x$status, kernel = x$kernel,
    N = x$N, N_market = x$N_market, N_market_year = x$N_market_year,
    N_group = x$N_group, K_s0 = x$K_s0, num_upars = x$schema$num_upars,
    lp = x$lp, lp_raw = x$lp_raw, grad_l2 = sqrt(sum(x$grad^2)),
    grad_raw_l2 = sqrt(sum(x$grad_raw^2)),
    fd_dim = length(x$fd$index), tp_len = length(x$tp), q_len = length(x$q),
    elapsed = x$elapsed, stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, rows), sub("\\.rds$", ".csv", out_file), row.names = FALSE)
