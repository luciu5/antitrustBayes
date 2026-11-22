test_that("normal and Warp-III bridges keep independent proposals, score diagnostics and RNG order", {
  set.seed(37)
  draws <- matrix(stats::rt(800, df = 8), ncol = 2L, dimnames = list(NULL, c("a", "b")))
  calls <- 0L
  args <- list(lb = c(a = -Inf, b = -Inf), ub = c(a = Inf, b = Inf), data = NULL,
    log_posterior = function(q, data) {
      calls <<- calls + 1L
      sum(stats::dt(q, df = 8, log = TRUE))
    })
  context <- list(samples = draws, bridge_args = args)
  run <- function(inputs) {
    set.seed(129)
    lapply(1:3, function(i) {
      if (!is.null(inputs$reset)) inputs$reset()
      do.call(bridgesampling::bridge_sampler, c(list(samples = inputs$samples,
        repetitions = 1L, method = method, silent = TRUE), inputs$bridge_args))
    })
  }
  for (method in c("normal", "warp3")) {
    calls <- 0L
    original <- run(context)
    original_seed <- .Random.seed
    expect_equal(calls, if (method == "normal") 1200 else 2400)
    calls <- 0L
    shared <- run(.bayes_shared_bridge_inputs(context, method))
    expect_identical(shared, original)
    expect_identical(.Random.seed, original_seed)
    expect_equal(calls, if (method == "normal") 800 else 1600)
  }
})

test_that("score reuse retries warnings, failed scores and mismatched coordinates", {
  draws <- matrix(1:8 + 0.0, ncol = 2)
  calls <- 0L
  ctx <- list(samples = draws, bridge_args = list(log_posterior = function(q, data) {
    calls <<- calls + 1L
    if (calls == 1) stop("retry")
    if (calls == 2) warning("score warning")
    if (calls == 3) -Inf else sum(q)
  }))
  inputs <- .bayes_shared_bridge_inputs(ctx, "normal")
  at <- function(x) {
    inputs$reset()
    inputs$bridge_args$log_posterior(x, NULL)
  }
  expect_error(at(draws[3, ]), "retry")
  expect_warning(at(draws[3, ]), "score warning")
  expect_identical(at(draws[3, ]), -Inf)
  expect_equal(at(draws[3, ]), 10)
  expect_equal(at(draws[3, ]), 10)
  expect_equal(at(draws[3, ] + 1), 12)
  expect_identical(calls, 5L)
})
