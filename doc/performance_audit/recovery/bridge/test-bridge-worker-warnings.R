test_that("fork bridge warnings remain caller-visible and attached to each job", {
  contexts <- list(list(samples = 1, stanfit_model = 1),
                   list(samples = 2, stanfit_model = 1))
  execute <- function(cores, fail = FALSE) {
    captured <- character()
    value <- withCallingHandlers(.bayes_bridge_dispatch(contexts, rep(NA_character_, 2L),
      c("one", "two"), bridge_repetitions = 2L, bridge_maxiter = 5000L,
      bridge_cores = cores, bridge_seed = 707101L,
      bridge_fun = function(samples, ...) {
        warning(paste("diagnostic", samples), call. = FALSE)
        if (fail && samples == 2) stop("synthetic failure after warning")
        list(logml = runif(1))
      }), warning = function(w) {
        captured <<- c(captured, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(value = value, captured = captured)
  }
  serial <- execute(1L)
  expect_equal(length(serial$captured), 4L)
  expect_equal(nrow(serial$value$warnings), 4L)
  expect_identical(serial$value$warnings$model, c("one", "one", "two", "two"))
  expect_identical(serial$value$warnings$repetition, c(1L, 2L, 1L, 2L))
  expect_equal(nrow(serial$value$failures), 0L)
  if (.Platform$OS.type != "windows") {
    forked <- execute(2L)
    expect_identical(forked$captured, serial$captured)
    expect_identical(forked$value$warnings, serial$value$warnings)
    expect_identical(forked$value$bridge_list, serial$value$bridge_list)
    broken <- execute(2L, fail = TRUE)
    expect_equal(nrow(broken$value$warnings), 4L)
    expect_equal(nrow(broken$value$failures), 2L)
    expect_identical(broken$value$failures$model, c("two", "two"))
    expect_true(all(vapply(broken$value$bridge_list[[2L]], is.null, logical(1))))
    expect_equal(length(broken$captured), 6L)
  }
})
