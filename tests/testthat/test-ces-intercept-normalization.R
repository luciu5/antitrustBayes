## CES counterpart of the Logit inside-market intercept correction
## (model specification 13). Outside-good utility is zero, product qualities
## sum to zero, and the single-market intercept is a free Normal(0, 1) draw.

.ces_r <- c(.25, .20, .16, .13, .11, .09, .06)
.ces_p <- c(1.2, 1.4, 1.6, 1.9, 2.1, 2.4, 2.8)

.ces_helmert <- function(G) {
  sapply(seq_len(G - 1L), function(k) {
    coordinate <- replace(numeric(G - 1L), k, 1)
    b <- numeric(G); suffix <- 0
    for (reverse_j in seq_len(G - 1L)) {
      j <- G + 1L - reverse_j; q <- j - 1L
      v <- coordinate[q] / sqrt(q * (q + 1))
      b[j] <- suffix - q * v; suffix <- suffix + v
    }
    b[1L] <- suffix
    b
  })
}

## Stan data for M markets x T periods of the seven-firm CES panel.
.ces_panel_data <- function(M = 1L, T = 1L, supply = 1L, outside_price = .8,
                            structure = "simultaneous", k = 3L, ...) {
  n <- 7L * M * T
  args <- list(prices = rep(.ces_p, M * T), quantities = rep(1000 * .ces_r, M * T),
    margins = rep(.3, n), market = rep(paste0("m", seq_len(M)), each = 7L * T),
    year = rep(rep(seq_len(T), each = 7L), M),
    group = paste0("f", seq_len(7L), "_", rep(seq_len(M * T), each = 7L)),
    supply_model = supply, use_cutoff = FALSE, use_rho = FALSE, use_hmt = FALSE,
    priors = list(), quality_parameterization = "sampled", quality_scope = "market",
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    outside_price = outside_price, structure = structure, share_basis = "absolute")
  if (!identical(structure, "simultaneous")) {
    args$partition_prior <- "ordered_k"; args$k_support <- k; args$k_prior <- 1
  }
  do.call(prepare_ces_stan_data, c(args, list(...)))
}

test_that("CES conditional shares are recovered at 10, 20 and 40 percent outside with a free intercept", {
  r <- .ces_r; p <- .ces_p
  z <- log(p) - mean(log(p))
  for (p0 in c(.8, 1, 1.4)) for (s0 in c(.10, .20, .40)) {
    d <- .ces_panel_data(outside_price = p0)
    gamma <- 1 + 0.7 / d$log_price_sd
    alpha_std <- gamma * d$log_price_sd # (a + log_price_sd) in the Stan share mean
    intercept <- mean(log(r)) - stats::qlogis(s0) - gamma * log(p0)
    b <- log(r) - mean(log(r)) + gamma * z
    reconstructed <- stats::qlogis(s0) + intercept + b - alpha_std * (z / d$log_price_sd) +
      gamma * log(p0)
    expect_equal(sum(b - gamma * z), 0, tolerance = 1e-13)
    expect_equal(reconstructed, log(r), tolerance = 1e-13)
    ## The removed single-market restriction (intercept = 0) pinned the outside
    ## share to plogis(mean(log(r)) - gamma * log(p0)); it cannot hold for all
    ## three chosen shares at once, whereas the free intercept attains each one.
    expect_equal(exp(reconstructed), r, tolerance = 1e-13)
    expect_gt(abs(intercept), 1e-3)
  }
  expect_lt(stats::plogis(mean(log(r))), 1 / (length(r) + 1))
})

test_that("CES single-market utility normalization is a sampled free intercept in both kernels", {
  for (kernel in c("ces_supply", "ces_stackelberg")) {
    code <- paste(readLines(.ab_test_path("inst", "stan", paste0(kernel, ".stan")),
                           warn = FALSE), collapse = "\n")
    expect_match(code, "vector[share_basis == 2 ? 1 : 0] mu_b_market_vec;", fixed = TRUE)
    expect_match(code, "market_intercept_location = rep_vector(mu_b_market_vec[1], N_market);",
                 fixed = TRUE)
    expect_match(code, "if (share_basis == 1) {\n    market_intercept_location = rep_vector(0.0, N_market);", fixed = TRUE)
    expect_length(gregexpr("normal_lpdf(mu_b_market_vec | 0, 1)", code,
                          fixed = TRUE)[[1L]], 1L)
    expect_match(code, "target += normal_lpdf(mu_b_market_vec | 0, 1);\n  if (is_single_market == 0)",
                 fixed = TRUE)
    expect_match(code, "vector<lower=0>[(share_basis == 2 && is_single_market == 0) ? 1 : 0] sigma_b_market_vec;",
                 fixed = TRUE)
    ## The reported rotated intercept now includes the single-market intercept.
    expect_false(grepl("if (is_single_market == 0) {\n      if (n_mky_per_market[m] > 1)",
                       code, fixed = TRUE))
  }
})

test_that("all four Stan programs still parse after the CES normalization", {
  testthat::skip_if_not_installed("cmdstanr")
  cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  testthat::skip_if(!nzchar(cmdstan_path), "CmdStan must be installed to parse Stan sources")
  stanc <- file.path(cmdstan_path, "bin", if (.Platform$OS.type == "windows") "stanc.exe" else "stanc")
  testthat::skip_if(!file.exists(stanc), "stanc is unavailable")
  for (model in c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg")) {
    source <- normalizePath(.ab_test_path("inst", "stan", paste0(model, ".stan")))
    output <- suppressWarnings(system2(stanc,
      args = c("--allow-undefined", paste0("--name=", model),
               paste0("--o=", tempfile(fileext = ".hpp")), shQuote(source)),
      stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    expect_true(is.null(status) || identical(as.integer(status), 0L),
                info = paste(output, collapse = "\n"))
  }
})

test_that("specification 14 rejects earlier fits, including Logit and CES single-market fits", {
  expect_identical(.antitrust_bayes_model_spec_version, 14L)
  for (demand in c("logit", "ces")) {
    prepared <- if (demand == "ces") .ces_panel_data() else prepare_stan_data(share_basis = "absolute",
      .ces_p, 1000 * .ces_r, rep(.3, 7), rep("m1", 7), rep(1, 7), paste0("f", 1:7),
      supply_model = 1L, market_side = "output", use_rho = FALSE)
    expect_identical(attr(prepared, "model_spec_version"), 14L)
    fit <- methods::new("BayesLogit", stan_data = prepared, supply_model = 1L,
      demand_model = demand, market_side = "output", price_sd = 1,
      call = quote(identity()))
    expect_invisible(validate_fit_spec(fit))
    old <- fit
    attr(old@stan_data, "model_spec_version") <- 12L
    expect_error(validate_fit_spec(old), "current specification is 14.*Refit it with bayes_est")
  }
  stale <- list(prior_normalization = list(mode = "native_stan_normalized",
                                           model_spec_version = 12L))
  fits <- list(methods::new("BayesLogit", stan_data = .ces_panel_data(), supply_model = 1L,
    demand_model = "ces", market_side = "output", price_sd = 1,
    call = quote(identity()), metadata = stale))
  expect_error(antitrustBayes:::.bayes_validate_bridge(fits), "native normalized")
  fits[[1L]]@metadata$prior_normalization$model_spec_version <- 14L
  expect_invisible(antitrustBayes:::.bayes_validate_bridge(fits))
})

test_that("the single-market flag that selects the free intercept reaches both Stan kernels", {
  captured <- NULL
  testthat::local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) list(
    variables = function() list(data = list(is_single_market = 1)),
    sample = function(data, ...) { captured <<- data; stop("captured before HMC") }),
    .package = "antitrustBayes")
  for (demand in c("logit", "ces")) {
    expect_error(bayes_est(share_basis = "absolute", prices = .ces_p, quantities = 1000 * .ces_r, margins = rep(.3, 7),
      market = rep("m1", 7), year = rep(1, 7), group = paste0("f", 1:7),
      supply_model = "bertrand", demand_model = demand, market_side = "output",
      use_rho = FALSE), "captured before HMC")
    expect_identical(as.integer(captured$is_single_market), 1L)
  }
})

## ---- Posterior-to-counterfactual reconstruction with a nonzero intercept ----

.ces_policy_fixture <- function(candidate = "B", outside_price = 1.4, shift = 0) {
  price <- c(1.2, 1.35, 1.5, 1.8)
  conditional <- c(.38, .29, .20, .13)
  outside_share <- .3; gamma <- 7; k <- 2L
  supply <- if (candidate == "B") 1L else 3L
  structure <- if (candidate == "CF") "fringe" else if (candidate == "CS") "stackelberg" else "simultaneous"
  quantities <- 1000 * conditional
  params <- bayes_params(gamma = gamma, s0 = outside_share, supply_model = supply,
    structure = structure, partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k",
    k_support = if (structure == "simultaneous") NULL else k,
    k_prior = if (structure == "simultaneous") NULL else 1)
  market <- bayes_market(price, quantities, rep(.3, 4L), group = letters[1:4],
    outside_price = outside_price, market_side = "output", insideSize = 1000)
  truth <- suppressWarnings(build_conduct_object(params, 1L, market,
    ownerPre = market$group, ownerPost = market$group))
  if (structure != "simultaneous") {
    truth <- suppressWarnings(posterior_policy_models(bayes_specify(params, market),
      draw_plan = posterior_policy_draws(bayes_specify(params, market), n_draws = 1L, seed = 1L),
      insideSize = 1000)$models[[1L]][[1L]])
    truth <- if (methods::is(truth, "AntitrustFit")) truth@model else truth
  }
  lerner <- as.numeric(suppressWarnings(antitrust::calcMargins(truth, preMerger = TRUE)))
  data <- data.frame(prices = price, quantities = quantities, margins = lerner,
    market = "single", year = "period", group = letters[1:4],
    outside_price = outside_price, margin_observed = TRUE)
  stan_data <- prepare_ces_stan_data(share_basis = "absolute", data$prices, data$quantities, data$margins,
    data$market, data$year, data$group, supply_model = supply, use_cutoff = FALSE,
    use_rho = FALSE, use_hmt = FALSE, priors = list(), outside_price = outside_price,
    structure = structure,
    partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k",
    k_support = if (structure == "simultaneous") NULL else k,
    k_prior = if (structure == "simultaneous") NULL else 1,
    quality_parameterization = "sampled", quality_scope = "market",
    share_likelihood = "gaussian", margin_likelihood = "gaussian")
  log_price_sd <- stan_data$log_price_sd
  z <- (log(price) - mean(log(price))) / log_price_sd
  intercept <- mean(log(conditional)) - stats::qlogis(outside_share) - gamma * log(outside_price)
  quality <- log(conditional) - mean(log(conditional)) + (gamma * log_price_sd) * z
  quality <- quality - mean(quality)
  values <- c(`a_market[1]` = (gamma - 1) * log_price_sd, `s0_prob[1]` = outside_share,
    `mu_b_market_vec[1]` = intercept + shift, `b_market[1]` = intercept + shift,
    `market_intercept_location[1]` = intercept + shift,
    `b_market_rot[1]` = intercept + shift, sigma_b_strat = .5, sigma_share = .2,
    sigma_margin = .3, `posterior_component_prob[1]` = 1, `prior_component_prob[1]` = 1)
  values <- c(values,
    stats::setNames(quality, paste0("b_quality_all[", 1:4, "]")),
    stats::setNames(log(conditional) + shift, paste0("share_mean[", 1:4, "]")),
    stats::setNames(1 / lerner, paste0("pred_marginNormInv[", 1:4, "]")))
  posterior <- array(rep(values, each = 2L), c(2L, 1L, length(values)),
    dimnames = list(iteration = 1:2, chain = "1", variable = names(values)))
  fit <- methods::new("BayesLogit", stanfit = structure(list(draws_array = posterior), class = "NutpieFit"),
    stan_data = stan_data, data = data, supply_model = supply, demand_model = "ces",
    market_side = "output", price_sd = log_price_sd,
    metadata = list(backend = "nutpie", inference = "hmc", structure = structure,
      partition_prior = attr(stan_data, "partition_prior"),
      firm_market_ranking = attr(stan_data, "firm_market_ranking"),
      k_support = if (structure == "simultaneous") NULL else k,
      k_prior = if (structure == "simultaneous") NULL else 1,
      quality_parameterization = "sampled", quality_scope = "market",
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      market_intercept_parameterization = "sampled", posterior_draws = 2L,
      model_spec_version = .antitrust_bayes_model_spec_version))
  list(fit = fit, truth = truth, lerner = lerner, price = price, quantities = quantities,
       conditional = conditional, outside_share = outside_share,
       cost = price * (1 - lerner), intercept = intercept + shift, quality = quality)
}

test_that("CES posterior draws with a nonzero intercept reproduce observed equilibria and costs", {
  skip_if_not_installed("coordination", "0.0.0.9006")
  for (outside_price in c(1, 1.4)) for (candidate in c("B", "C", "CF", "CS")) {
    fixture <- .ces_policy_fixture(candidate, outside_price)
    info <- paste(candidate, "outside price", outside_price)
    expect_gt(abs(fixture$intercept), .1)
    expect_true(all(fixture$cost > 0), info = info)
    plan <- posterior_policy_draws(fixture$fit, n_draws = 2L, seed = 915L)
    for (anchor in c("observed", "latent_quality")) {
      route <- suppressWarnings(posterior_policy_models(fixture$fit, draw_plan = plan,
        demand_anchoring = anchor, insideSize = 1000))
      expect_equal(nrow(route$failures), 0L, info = info)
      for (position in seq_along(plan$draw_indices)) {
        result <- route$models[[1L]][[position]]
        model <- if (methods::is(result, "AntitrustFit")) result@model else result
        expect_equal(as.numeric(model@pricePre), fixture$price, tolerance = 1e-6, info = info)
        expect_equal(as.numeric(antitrust::calcShares(model, TRUE, revenue = FALSE)),
          (1 - fixture$outside_share) * fixture$conditional, tolerance = 1e-6, info = info)
        expect_equal(as.numeric(model@mcPre), fixture$cost, tolerance = 1e-6, info = info)
        expect_equal(as.numeric(antitrust::calcMargins(model, TRUE)), fixture$lerner,
          tolerance = 1e-6, info = info)
      }
    }
  }
})

test_that("a CES intercept shift moves the share mean once and leaves the equilibrium unchanged", {
  base <- .ces_policy_fixture("B", 1.4)
  shifted <- .ces_policy_fixture("B", 1.4, shift = .45)
  plan <- posterior_policy_draws(base$fit, n_draws = 2L, seed = 822L)
  other <- posterior_policy_draws(shifted$fit, n_draws = 2L, seed = 822L)
  first <- suppressWarnings(posterior_predict(base$fit, draw_plan = plan, seed = 610L))
  second <- suppressWarnings(posterior_predict(shifted$fit, draw_plan = other, seed = 610L))
  expect_equal(unname(second$replicated$log_shareIn - first$replicated$log_shareIn),
    matrix(.45, 2L, 4L), tolerance = 1e-12)
  expect_equal(second$replicated$share_rep_normalized, first$replicated$share_rep_normalized,
    tolerance = 1e-12)
  expect_equal(second$replicated$log_marginNormInv, first$replicated$log_marginNormInv,
    tolerance = 1e-12)
})

## ---- Opt-in Stan execution at fixed parameters ----

.skip_ces_stan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if(!identical(tolower(Sys.getenv("ANTITRUSTBAYES_RUN_STAN_TESTS")), "true"),
    "set ANTITRUSTBAYES_RUN_STAN_TESTS=true to run Stan execution tests")
}

.ces_fixed_point <- function(kernel, data, init, variables) {
  withr::local_options(cmdstanr_warn_inits = FALSE)
  model <- antitrustBayes:::.cmdstan_model_for_kernel(kernel)
  fit <- suppressWarnings(model$sample(data = data, fixed_param = TRUE, chains = 1L,
    iter_sampling = 1L, init = list(init), refresh = 0L, threads_per_chain = 1L,
    sig_figs = 18L, show_messages = FALSE))
  stats::setNames(lapply(variables, function(v)
    as.numeric(posterior::as_draws_matrix(fit$draws(v)))), variables)
}

## Stan evaluation of lp__ at a point: a near-zero step, no adaptation.
.ces_log_density <- function(kernel, data, init) {
  withr::local_options(cmdstanr_warn_inits = FALSE)
  model <- antitrustBayes:::.cmdstan_model_for_kernel(kernel)
  fit <- suppressWarnings(model$sample(data = data, chains = 1L, iter_warmup = 0L,
    iter_sampling = 1L, init = list(init), refresh = 0L, threads_per_chain = 1L,
    adapt_engaged = FALSE, step_size = 1e-12, sig_figs = 18L, show_messages = FALSE))
  as.numeric(posterior::as_draws_matrix(fit$draws("lp__")))
}

test_that("CES kernels recover the conditional shares and the R margins at fixed parameters", {
  .skip_ces_stan()
  skip_if_not_installed("coordination", "0.0.0.9006")
  r <- .ces_r; p <- .ces_p; H <- .ces_helmert(7L)
  log_a <- .7; k <- 3L
  for (kernel in c("ces_supply", "ces_stackelberg")) for (p0 in c(.8, 1.4)) {
    for (s0 in c(.1, .2, .4)) for (supply in c(1L, 3L)) {
      stack <- identical(kernel, "ces_stackelberg")
      d <- .ces_panel_data(supply = supply, outside_price = p0,
        structure = if (stack) "stackelberg" else "simultaneous", k = k)
      gamma <- 1 + exp(log_a) / d$log_price_sd
      z <- log(p) - mean(log(p))
      intercept <- mean(log(r)) - stats::qlogis(s0) - gamma * log(p0)
      b <- log(r) - mean(log(r)) + gamma * z
      init <- list(mu_b_market_vec = array(intercept, 1L),
        b_group_centered = array(as.numeric(crossprod(H, b)), 6L),
        sigma_b_strat = 1, mu_log_a = log_a, s0_logit = array(stats::qlogis(s0), 1L),
        log_sigma_share = log(.1), log_sigma_margin = log(.1))
      out <- .ces_fixed_point(kernel, d, init,
        c("share_mean", "s0_prob", "pred_marginNormInv"))
      info <- paste(kernel, "p0", p0, "s0", s0, "supply", supply)
      ## Stan maps s0 through eps + (1 - 2 eps) * inv_logit with eps = 1e-6, so the
      ## realized outside share, and with it the share mean, shifts by O(eps).
      bounded <- 1e-6 + (1 - 2e-6) * s0
      expect_equal(out$s0_prob, bounded, tolerance = 1e-12, info = info)
      expect_equal(out$share_mean,
        log(r) + stats::qlogis(bounded) - stats::qlogis(s0), tolerance = 1e-9, info = info)
      params <- if (stack) bayes_params(gamma = gamma, s0 = bounded, supply_model = supply,
        structure = "stackelberg", partition_prior = "ordered_k", k_support = k, k_prior = 1)
      else bayes_params(gamma = gamma, s0 = bounded, supply_model = supply)
      market <- bayes_market(p, 1000 * r, rep(.3, 7L), group = paste0("f", 1:7),
        outside_price = p0, market_side = "output", insideSize = 1000)
      model <- if (stack) {
        fit <- bayes_specify(params, market)
        got <- suppressWarnings(posterior_policy_models(fit,
          draw_plan = posterior_policy_draws(fit, n_draws = 1L, seed = 1L), insideSize = 1000))
        result <- got$models[[1L]][[1L]]
        if (methods::is(result, "AntitrustFit")) result@model else result
      } else suppressWarnings(build_conduct_object(params, 1L, market,
        ownerPre = market$group, ownerPost = market$group))
      expect_equal(out$pred_marginNormInv,
        1 / as.numeric(suppressWarnings(antitrust::calcMargins(model, TRUE))),
        tolerance = 1e-6, info = info)
    }
  }
})

test_that("CES single-period, multi-period and multi-market layouts keep per-market-year outside revenue shares", {
  .skip_ces_stan()
  r <- .ces_r; p <- .ces_p; s0 <- .2
  for (layout in list(c(1L, 4L), c(4L, 1L), c(4L, 4L))) {
    M <- layout[[1L]]; T <- layout[[2L]]
    p0 <- rep(c(.8, 1, 1.2, 1.4), length.out = M * T)
    d <- .ces_panel_data(M = M, T = T, outside_price = p0)
    cells <- M * T
    init <- list(mu_b_market_vec = array(.2, 1L),
      b_group_centered = array(0, d$N_quality_coord), sigma_b_strat = 1, mu_log_a = .3,
      s0_logit = array(stats::qlogis(s0), cells), log_sigma_share = log(.1),
      log_sigma_margin = log(.1))
    if (M > 1L) init <- c(init, list(sigma_b_market_vec = array(.4, 1L), b_market_raw = rep(0, M)))
    out <- .ces_fixed_point("ces_supply", d, init, c("s0_prob", "R0_prob", "share_mean"))
    inside_index <- as.numeric(d$inside_price_index)
    bounded <- 1e-6 + (1 - 2e-6) * s0
    expected <- stats::plogis(stats::qlogis(bounded) + log(p0) - log(inside_index))
    info <- paste0("M", M, "T", T)
    expect_equal(out$s0_prob, rep(bounded, cells), tolerance = 1e-12, info = info)
    expect_equal(out$R0_prob, expected, tolerance = 1e-10, info = info)
    ## CES quantity share and revenue share differ unless p0 equals the index.
    expect_gt(max(abs(out$R0_prob - s0)), 1e-3)
  }
})

test_that("the CES intercept prior is one normalized Normal(0, 1) in every layout", {
  .skip_ces_stan()
  for (M in c(1L, 4L)) for (kernel in c("ces_supply", "ces_stackelberg")) {
    stack <- identical(kernel, "ces_stackelberg")
    d <- .ces_panel_data(M = M, T = 1L, likelihood_mode = 3L,
      structure = if (stack) "stackelberg" else "simultaneous", k = 3L)
    if (stack) next # the ordered-k target adds a mixture over k; text-match test covers its prior
    mu <- c(-1, 0, .7, 1.5)
    density <- vapply(mu, function(x) {
      init <- list(mu_b_market_vec = array(x, 1L),
        b_group_centered = array(0, d$N_quality_coord), sigma_b_strat = 1, mu_log_a = .3,
        s0_logit = array(stats::qlogis(.2), M), log_sigma_share = log(.1),
        log_sigma_margin = log(.1))
      if (M > 1L) init <- c(init, list(sigma_b_market_vec = array(.4, 1L), b_market_raw = rep(0, M)))
      .ces_log_density(kernel, d, init)
    }, numeric(1))
    ## Margin-only mode removes the share equation, so mu enters only through its
    ## prior. The 1e-6 allowance is the bounded outside-share link: the margin term
    ## depends on mu through the market intercept's rotation of s0 (b_market).
    expect_equal(diff(density), diff(stats::dnorm(mu, log = TRUE)), tolerance = 1e-5,
      info = paste("M", M))
  }
  sizes <- antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")$variables()$parameters
  expect_equal(as.integer(sizes$mu_b_market_vec$dimensions), 1L)
})
