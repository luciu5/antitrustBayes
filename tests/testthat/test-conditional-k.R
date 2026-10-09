test_that("conditional calls retain observations and priors across rich conduct", {
  for (structure in c("fringe", "stackelberg")) for (supply in c("bertrand", "cournot")) {
    args <- list(prices = rep(1:7, 4), quantities = rep(7:1, 4), margins = rep(.3, 28),
      market = rep("one", 28), year = rep(1:4, each = 7), group = rep(letters[1:7], 4),
      structure = structure, supply_model = supply, k_scope = "market",
      priors = list(prior_alpha_sd = 2), share_likelihood = "gaussian",
      stan_options = list(chains = 4L, iter = 4000L, seed = 12L))
    prior <- (1:6) / sum(1:6)
    calls <- do.call(bayes_conditional_k_arguments, c(args, list(k_support = 1:6, k_prior = prior)))
    expect_equal(unname(attr(calls, "k_prior")), prior)
    for (k in 1:6) {
      x <- calls[[as.character(k)]]
      expect_identical(x$k_support, k); expect_identical(x$k_prior, 1)
      x$k_support <- x$k_prior <- x$partition_prior <- NULL
      expect_identical(x, args)
    }
  }
})

test_that("conditional enumeration rejects unsupported market/state definitions", {
  args <- list(market = "m", group = "g", structure = "fringe", supply_model = "cournot")
  expect_error(do.call(bayes_conditional_k_arguments, c(args, list(k_support = 2))), "feasible")
  args$market <- c("m", "n"); args$group <- c("a", "b")
  expect_error(do.call(bayes_conditional_k_arguments, c(args, list(k_support = 1))), "one market")
  args$market <- c("m", "m"); args$k_scope <- "market_year"
  expect_error(do.call(bayes_conditional_k_arguments, c(args, list(k_support = 1))), "market-scoped")
  args$k_scope <- "market"; args$structure <- "simultaneous"
  expect_error(do.call(bayes_conditional_k_arguments, c(args, list(k_support = 1))), "hypothesis")
})

test_that("conditional evidence matches direct integration with nonuniform k prior", {
  theta_prior <- c(.2, .3, .5)
  likelihood <- matrix(seq(.05, .9, length.out = 18), nrow = 3)
  z <- colSums(likelihood * theta_prior); names(z) <- as.character(1:6)
  prior <- (1:6) / sum(1:6)
  repetitions <- rbind(log(z), log(z), log(z))
  result <- bayes_combine_conditional_k(log(z), k_support = 1:6,
    k_prior = prior, repetitions = repetitions)
  expect_equal(exp(result$logml), sum(z * prior), tolerance = 1e-14)
  expect_equal(result$posterior_k$probability, unname(z * prior / sum(z * prior)))
  expect_equal(result$logml_sd, 0)
  expect_equal(bayes_combine_conditional_k(log(z) + 10000, k_support = 1:6, k_prior = prior)$logml,
    result$logml + 10000)
  expect_error(bayes_combine_conditional_k(log(z[-1]), k_support = 1:6), "Every original")
  bad <- log(z); bad[3] <- NA_real_
  expect_error(bayes_combine_conditional_k(bad, k_support = 1:6), "finite named")
  expect_error(bayes_combine_conditional_k(log(z), k_support = 1:6,
    repetitions = repetitions[, -1]), "all original")
  zero <- prior; zero[1] <- 0; zero <- zero / sum(zero)
  expect_equal(bayes_combine_conditional_k(log(z), k_support = 1:6,
    k_prior = zero)$posterior_k$probability[1], 0)
})

test_that("component scheduling preserves failed fits and explicit numerical options", {
  calls <- bayes_conditional_k_arguments(market = rep("m", 3), group = letters[1:3],
    structure = "stackelberg", supply_model = "bertrand", k_support = 1:3,
    stan_options = list(chains = 4L, seed = 10L), component_stan_options =
      list("1" = list(seed = 11L), "2" = list(seed = 12L), "3" = list(seed = 13L)))
  expect_identical(vapply(calls, function(x) x$stan_options$seed, integer(1)), c("1" = 11L, "2" = 12L, "3" = 13L))
  local_mocked_bindings(bayes_est = function(...) {
    args <- list(...)
    if (args$k_support == 2L) stop("retained failed component")
    args
  })
  result <- bayes_est_conditional_k(market = rep("m", 3), group = letters[1:3],
    structure = "stackelberg", supply_model = "bertrand", k_support = 1:3)
  expect_identical(names(result$fits), as.character(1:3))
  expect_s3_class(result$fits[[2]], "error")
  expect_identical(result$fits[[3]]$k_support, 3L)
})
