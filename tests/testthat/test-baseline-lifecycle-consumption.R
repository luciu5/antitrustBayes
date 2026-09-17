.count_baseline_price_calls <- function(expr) {
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  counter_expr <- substitute(
    assign("n", get("n", envir = COUNTER) + 1L, envir = COUNTER),
    list(COUNTER = counter)
  )
  suppressMessages(trace("calcPrices", tracer = counter_expr,
                         where = asNamespace("antitrust"), print = FALSE))
  on.exit(suppressMessages(untrace("calcPrices",
                                  where = asNamespace("antitrust"))), add = TRUE)
  value <- force(expr)
  list(value = value, calls = counter$n)
}

test_that("posterior realization uses the public observed-baseline lifecycle", {
  skip_if_not("baseline" %in% names(formals(antitrust::specify)),
              "installed antitrust lacks the optional baseline lifecycle")
  market <- bayes_market(
    prices = c(2, 2.2, 2.5), quantities = c(40, 35, 25),
    margins = c(NA_real_, .2, NA_real_), group = c("A", "A", "B"),
    insideSize = 100, outside_price = 1
  )
  for (demand in c("logit", "ces")) for (supply in c(1L, 3L)) {
    pars <- if (demand == "logit") {
      bayes_params(alpha = 1.2, s0 = .2, supply_model = supply)
    } else {
      bayes_params(gamma = 3, s0 = .2, supply_model = supply)
    }
    object <- suppressWarnings(build_conduct_object(
      pars, 1L, market, ownerPre = market$group, ownerPost = market$group
    ))
    observed <- .count_baseline_price_calls(
      suppressWarnings(as_antitrust_fit(object))
    )
    solved <- .count_baseline_price_calls(
      suppressWarnings(as_antitrust_fit(object, baseline = "solve"))
    )
    expect_identical(observed$calls, 0L)
    expect_identical(solved$calls, 1L)
    expect_identical(class(observed$value@model), class(solved$value@model))
    expect_equal(observed$value@model@mcPre, solved$value@model@mcPre)
    expect_equal(observed$value@model@mktSize, solved$value@model@mktSize)
    expect_equal(attr(observed$value@model, "antitrust_cost_state"),
                 attr(solved$value@model, "antitrust_cost_state"))
    for (scenario in c("merger", "cost", "combined")) {
      owners <- if (scenario == "cost") market$group else rep("A", 3L)
      shock <- rep(if (scenario == "merger") 0 else .1, 3L)
      old <- suppressWarnings(antitrust::simulate(
        solved$value, owners, mcDelta = shock
      ))
      new <- suppressWarnings(antitrust::simulate(
        observed$value, owners, mcDelta = shock
      ))
      expect_equal(new@pricePost, old@pricePost, tolerance = 1e-7)
      expect_equal(antitrust::calcShares(new, FALSE),
                   antitrust::calcShares(old, FALSE), tolerance = 1e-8)
      expect_equal(antitrust::CV(new), antitrust::CV(old), tolerance = 1e-6)
      expect_equal(antitrust::calcProducerSurplus(new, FALSE),
                   antitrust::calcProducerSurplus(old, FALSE), tolerance = 1e-6)
    }
  }
})
