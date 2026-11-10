test_that("the public supplied-state merger runs quietly without CmdStan", {
  withr::local_envvar(CMDSTAN = tempfile("absent-cmdstan-"))
  data("pnb", package = "antitrustBayes", envir = environment())
  params <- bayes_params(alpha = 1, s0 = 0.015, supply_model = 1L)
  market <- bayes_market(
    prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
    group = pnb$group, outside_price = 1, market_side = "input",
    insideSize = attr(pnb, "insideSize"), labels = pnb$bank
  )
  fit <- bayes_specify(params, market)
  owner_post <- pnb$group
  owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
  warnings <- character()
  result <- withCallingHandlers(
    simulate(fit, antitrust::counterfactual(ownership = owner_post), mc.cores = 1L),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
  expect_identical(warnings, character())
  expect_s4_class(result, "BayesSim")
  expect_identical(result@metadata$demand_anchoring, "observed")
  pre <- result@price_pre[[1L]]
  post <- result@price_post[[1L]]
  expect_equal(as.numeric(pre[1, ]), pnb$prices, tolerance = 1e-8)
  expect_true(all(is.finite(post)))
  expect_gt(max(abs(post - pre)), 1e-6)
})

test_that("baseline construction retains unrelated validation warnings", {
  market <- bayes_market(prices = c(1, 1.2, .9), quantities = c(20, 15, 10),
    margins = rep(.2, 3), group = c("a", "a", "b"), market_side = "output",
    outside_price = 1, insideSize = 1)
  params <- bayes_params(alpha = 2, s0 = .2, supply_model = 1L)
  original_new <- methods::new
  warned <- FALSE
  local_mocked_bindings(new = function(...) {
    if (!warned) {
      warned <<- TRUE
      warning("material numerical warning fixture")
    }
    original_new(...)
  }, .package = "methods")
  expect_warning(build_conduct_object(params, 1L, market,
    ownerPre = market$group, ownerPost = market$group),
    "material numerical warning fixture")
})

## Material failure: apparently ignored shares also initialize market mass.
## Correct prices alone miss volume and welfare errors from the placeholder.
test_that("supplied Logit and CES baselines preserve independent quantity totals", {
  q <- c(40, 30, 20)
  for (side in c("output", "input")) for (supply in c(1L, 4L)) {
    market <- bayes_market(prices = c(1.1, 1.2, 1.3), quantities = q,
      margins = c(.5, .4, .3), group = letters[1:3], outside_price = 1,
      market_side = side, insideSize = 90)
    params <- bayes_params(alpha = 6, s0 = .2, supply_model = supply)
    prepared <- posterior_policy_models(bayes_specify(params, market),
      n_draws = 1L, seed = 12L)
    expect_equal(nrow(prepared$failures), 0L)
    model <- prepared$models[[1L]][[1L]]@model
    expect_equal(model@mktSize, 90 / .8, tolerance = 1e-10)
    expect_equal(unname(antitrust::calcQuantities(model, TRUE)),
      90 * q / sum(q), tolerance = 1e-10)
  }
  market$market_side <- "output"
  params <- bayes_params(gamma = 10, s0 = .2, supply_model = 1L)
  prepared <- posterior_policy_models(bayes_specify(params, market),
    n_draws = 1L, seed = 12L)
  model <- prepared$models[[1L]][[1L]]@model
  ## Total CES expenditure includes outside-good quantity at its own price.
  expected_budget <- 90 * (sum(market$prices * q / sum(q)) + .2 / .8)
  expect_equal(model@mktSize, expected_budget, tolerance = 1e-10)
  expect_equal(unname(antitrust::calcQuantities(model, TRUE)),
    90 * q / sum(q), tolerance = 1e-10)
})

test_that("the registered baseline adapter retains unrelated constructor warnings", {
  market <- bayes_market(prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), group = letters[1:3], outside_price = 1,
    market_side = "output", insideSize = 90)
  params <- bayes_params(alpha = 6, s0 = .2, supply_model = 1L)
  object <- build_conduct_object(params, 1L, market, market$group, market$group)
  original_specify <- antitrust::specify
  local_mocked_bindings(specify = function(...) {
    warning("material registered-constructor warning fixture")
    original_specify(...)
  }, .package = "antitrust")
  expect_warning(as_antitrust_fit(object),
    "material registered-constructor warning fixture")
})

test_that("the documented estimation call reaches the backend before compilation", {
  readme <- test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme),
    "README extraction is repository-only; public backend dispatch has separate tests.")
  data("pnb", package = "antitrustBayes", envir = environment())
  lines <- readLines(readme, warn = FALSE)
  start <- grep("^estimated_fit <- bayes_est\\(", lines)
  expect_length(start, 1L)
  end <- start + which(lines[(start + 1L):length(lines)] == ")")[[1L]]
  local_mocked_bindings(.cmdstan_model_for_kernel = function(...) {
    stop("README stopped before compilation")
  })
  expect_error(eval(parse(text = lines[start:end])),
    "README stopped before compilation")
})
