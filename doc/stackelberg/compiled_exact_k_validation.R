# Fixed-parameter validation of the compiled exact ordered-k Stan programs.
# Run after installing the package from the current source tree.

suppressPackageStartupMessages(library(antitrustBayes))
source(file.path("doc", "stackelberg", "exact_k_reference.R"))

.lse <- function(x) {
  z <- max(x)
  z + log(sum(exp(x - z)))
}

.draw_vector <- function(x) as.numeric(drop(x[1, , drop = TRUE]))

.student_log_density <- function(y, mu, scale, nu) {
  stats::dt((y - mu) / scale, df = nu, log = TRUE) - log(scale)
}

.fixture <- function() {
  product <- rep(c("a1", "a2", "b", "c", "d"), 3L)
  group <- rep(c("A", "A", "B", "C", "D"), 3L)
  market <- rep(c("m1", "m1", "m2"), each = 5L)
  year <- rep(c("y1", "y2", "y1"), each = 5L)
  prices <- c(2.00, 2.40, 1.80, 2.10, 2.20,
              2.05, 2.35, 1.90, 2.15, 2.25,
              1.70, 2.30, 2.00, 2.45, 1.95)
  quantities <- c(30, 10, 22, 13, 8,
                  16,  8, 18, 30, 12,
                  12,  7, 29, 11, 21)
  margins <- c(.24, .20, .27, .22, .18,
               .23, .19, .25, .21, .17,
               .20, .18, .26, .16, .24)
  margin_observed <- rep(TRUE, length(prices))
  margin_observed[c(2L, 13L)] <- FALSE
  list(product = product, group = group, market = market, year = year,
       prices = prices, quantities = quantities, margins = margins,
       margin_observed = margin_observed)
}

.prepare <- function(x, demand, conduct, structure, likelihood_mode = 1L) {
  support <- if (structure == "fringe") c(3L, 0L, 2L) else c(3L, 1L, 2L)
  prior <- c(.20, .50, .30)
  common <- list(
    prices = x$prices, quantities = x$quantities, margins = x$margins,
    market = x$market, year = x$year, group = x$group,
    supply_model = if (conduct == "bertrand") 1L else 3L,
    use_cutoff = FALSE, use_rho = FALSE, use_market_resid_scale = FALSE,
    use_hmt = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_s0_market_pool = FALSE,
    use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
    priors = list(), margin_observed = x$margin_observed,
    likelihood_mode = likelihood_mode, structure = structure,
    partition_prior = "ordered_k", k_support = support, k_prior = prior,
    k_scope = "market"
  )
  if (demand == "logit") {
    do.call(antitrustBayes:::prepare_stan_data,
            c(common, list(market_side = "output", outside_price = 1)))
  } else {
    do.call(antitrustBayes:::prepare_ces_stan_data,
            c(common, list(market_side = "output", outside_price = 1)))
  }
}

.conditional_inverse <- function(d, ex, demand, conduct, structure, cc, rows) {
  k <- d$partition_component_k[cc]
  ans <- numeric(length(rows))
  for (mky in unique(d$market_year_idx[rows])) {
    ii <- rows[d$market_year_idx[rows] == mky]
    fm <- d$obs_to_firm_mky[ii]
    firm <- as.character(d$firm_mky_group[fm])
    rank <- d$firm_mky_rank[fm]
    market_index <- d$market[ii[[1L]]]
    s0 <- .draw_vector(ex$s0_prob)[mky]
    if (demand == "logit") {
      ans[match(ii, rows)] <- exact_k_logit_inverse_lerner(
        price = d$price[ii], firm = firm,
        firm_share_in = d$firm_mky_shareIn[fm], rank = rank,
        s0 = s0, alpha = .draw_vector(ex$alpha_structural_market)[market_index],
        k = k, structure = structure, conduct = conduct
      )
    } else {
      ans[match(ii, rows)] <- exact_k_ces_inverse_lerner(
        firm = firm, firm_rev_share_in = d$firm_mky_revShare[fm], rank = rank,
        outside_quantity_share = s0, outside_price = d$outside_price[mky],
        inside_price_index = d$inside_price_index[mky],
        gamma = .draw_vector(ex$gamma_structural_market)[market_index],
        k = k, structure = structure, conduct = conduct
      )
    }
  }
  ans
}

.run_case <- function(demand, conduct, structure, likelihood_mode = 1L) {
  x <- .fixture()
  d <- .prepare(x, demand, conduct, structure, likelihood_mode)
  model_name <- paste0(demand, "_stackelberg")
  demand_par <- if (demand == "logit") "alpha_structural_market" else
    "gamma_structural_market"
  fit <- rstan::sampling(
    antitrustBayes:::stanmodels[[model_name]], data = d,
    chains = 1L, iter = 1L, warmup = 0L, algorithm = "Fixed_param",
    init = 0, refresh = 0,
    pars = c("s0_prob", demand_par, "sigma_margin_market", "nu",
             "posterior_component_prob", "prior_component_prob",
             "component_log_score", "pred_marginNormInv", "pred_margin",
             "unit_margin_loglik", "posterior_k_mean", "posterior_k_median",
             "posterior_k_mode", "log_lik_margin_market", "log_lik_market",
             "log_lik", "firm_market_core_prob")
  )
  stopifnot(length(fit@sim$samples) == 1L)
  ex <- rstan::extract(fit)
  ncomp <- d$N_partition_component
  expected_score <- numeric(ncomp)
  expected_post <- numeric(ncomp)
  expected_inv <- expected_margin <- numeric(d$N)
  expected_unit_loglik <- numeric(d$N_partition_unit)
  expected_mean <- expected_median <- expected_mode <- numeric(d$N_partition_unit)
  sigma_margin <- .draw_vector(ex$sigma_margin_market)
  nu <- as.numeric(ex$nu[1])

  for (u in seq_len(d$N_partition_unit)) {
    cc <- seq.int(d$partition_unit_start[u], length.out = d$partition_unit_count[u])
    mm <- d$partition_unit_market[u]
    rows <- which(d$market == mm)
    inv <- vapply(cc, function(j) {
      .conditional_inverse(d, ex, demand, conduct, structure, j, rows)
    }, numeric(length(rows)))
    score <- d$partition_prior_logprob[cc]
    if (likelihood_mode != 2L) {
      observed <- d$margin_obs[rows] == 1L
      scale <- sigma_margin[mm] * sqrt((nu - 2) / nu)
      for (j in seq_along(cc)) {
        score[j] <- score[j] + sum(.student_log_density(
          d$log_marginNormInv[rows][observed],
          log(inv[observed, j] + 1e-9), scale, nu
        ))
      }
    }
    post <- exp(score - .lse(score))
    expected_score[cc] <- score
    expected_post[cc] <- post
    expected_inv[rows] <- as.numeric(inv %*% post)
    expected_margin[rows] <- as.numeric((1 / inv) %*% post)
    expected_unit_loglik[u] <- if (likelihood_mode == 2L) 0 else .lse(score)
    kval <- d$partition_component_k[cc]
    expected_mean[u] <- sum(kval * post)
    ord <- order(kval)
    expected_median[u] <- kval[ord][which(cumsum(post[ord]) >= .5)[1L]]
    expected_mode[u] <- kval[which.max(post)]
  }

  got_score <- .draw_vector(ex$component_log_score)
  got_post <- .draw_vector(ex$posterior_component_prob)
  got_prior <- .draw_vector(ex$prior_component_prob)
  errors <- c(
    component_score = max(abs(got_score - expected_score)),
    component_probability = max(abs(got_post - expected_post)),
    prior_probability = max(abs(got_prior - d$partition_prior_probability)),
    inverse_margin = max(abs(.draw_vector(ex$pred_marginNormInv) - expected_inv)),
    lerner_margin = max(abs(.draw_vector(ex$pred_margin) - expected_margin)),
    mixture_loglik = max(abs(.draw_vector(ex$unit_margin_loglik) - expected_unit_loglik)),
    k_mean = max(abs(.draw_vector(ex$posterior_k_mean) - expected_mean)),
    k_median = max(abs(.draw_vector(ex$posterior_k_median) - expected_median)),
    k_mode = max(abs(.draw_vector(ex$posterior_k_mode) - expected_mode)),
    market_telescope = max(vapply(seq_len(d$N_market), function(mm) {
      abs(.draw_vector(ex$log_lik_market)[mm] -
            sum(.draw_vector(ex$log_lik)[d$market == mm]))
    }, numeric(1)))
  )
  stopifnot(max(errors) < 2e-10)
  data.frame(demand = demand, conduct = conduct, structure = structure,
             likelihood_mode = likelihood_mode, check = names(errors),
             max_abs_error = as.numeric(errors), row.names = NULL)
}

results <- do.call(rbind, c(
  lapply(c("logit", "ces"), function(demand) {
    do.call(rbind, lapply(c("bertrand", "cournot"), function(conduct) {
      do.call(rbind, lapply(c("fringe", "stackelberg"), function(structure) {
        .run_case(demand, conduct, structure)
      }))
    }))
  }),
  list(.run_case("logit", "bertrand", "stackelberg", likelihood_mode = 2L))
))

utils::write.csv(results,
                 file.path("doc", "stackelberg", "compiled_exact_k_validation.csv"),
                 row.names = FALSE)
cat("PASS:", nrow(results), "compiled exact-k checks; maximum absolute error",
    format(max(results$max_abs_error), scientific = TRUE), "\n")
