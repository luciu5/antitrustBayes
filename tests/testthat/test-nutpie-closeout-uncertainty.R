closeout_uncertainty_runner <- function() {
  path <- testthat::test_path("..", "..", "doc", "validation", "nutpie_closeout_uncertainty_audit.R")
  if (!file.exists(path)) skip("Research audit is excluded from built package.")
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment
}

test_that("policy uncertainty sums source errors without assuming independent quality noise", {
  runner <- closeout_uncertainty_runner()
  component <- matrix(c(.01, .02, .03, .04, .05), 5L, dimnames = list(NULL, "mean"))
  bound <- runner$nutpie_uncertainty_bound(component, c(mean = .02))
  expect_equal(bound, c(mean = .17))
  expect_gt(bound[[1L]], sqrt(sum(component^2)) + .02)
  component[1L, 1L] <- NA_real_
  expect_true(is.na(runner$nutpie_uncertainty_bound(component, c(mean = .02))))
  expect_error(runner$nutpie_uncertainty_bound(component[-1L, , drop = FALSE], c(mean = .02)),
    "every original k component")
})

test_that("uncertainty audit perturbs one component while retaining every chain axis", {
  withr::local_preserve_seed()
  runner <- closeout_uncertainty_runner()
  values <- lapply(1:5, function(k) matrix(rep(k, 4000L), 1000L, 4L))
  values[[3L]] <- matrix(seq_len(4000L), 1000L, 4L)
  statistic <- function(changed, weights) {
    stopifnot(all(vapply(c(1L, 2L, 4L, 5L), function(k) identical(values[[k]], changed[[k]]), logical(1))),
      all(vapply(changed, function(x) identical(dim(x), c(1000L, 4L)), logical(1))))
    c(mean = sum(weights * vapply(changed, mean, numeric(1))))
  }
  result <- runner$nutpie_uncertainty_component(values, rep(.2, 5L), 3L, statistic,
    repetitions = 5L, seed = 123L)
  expect_equal(result$point, statistic(values, rep(.2, 5L)))
  expect_gt(result$se[[1L]], 0)
  expect_identical(dim(result$bootstrap), c(5L, 1L))
  expect_identical(result$block_length, 32L)
})

test_that("revised comparison sums backend bounds and keeps exact identities explicit", {
  runner <- closeout_uncertainty_runner()
  a <- list(estimate = c(mean = 1), bound = c(mean = .02), constant = c(mean = FALSE))
  b <- list(estimate = c(mean = 1.03), bound = c(mean = .03), constant = c(mean = FALSE))
  result <- runner$nutpie_uncertainty_difference(a, b, c(mean = .25))
  expect_equal(result$combined_mcse_bound, .05)
  expect_true(result$agreement_passed)
  expect_true(result$precision_passed)
  a$constant[] <- b$constant[] <- TRUE
  a$bound[] <- b$bound[] <- NA_real_
  b$estimate <- a$estimate
  result <- runner$nutpie_uncertainty_difference(a, b, c(mean = .25))
  expect_true(result$exact_constant_identity)
  expect_true(result$precision_passed)
  expect_true(is.na(result$combined_mcse_bound))
  b$estimate[] <- 2
  result <- runner$nutpie_uncertainty_difference(a, b, c(mean = .25))
  expect_false(result$agreement_passed)
  expect_false(result$precision_passed)
})

test_that("saved policy uncertainty restores both chain axes and retains paired failures", {
  runner <- closeout_uncertainty_runner()
  rows <- data.frame(regime = rep(c("merger", "tariff25"), each = 4000L),
    draw_index = rep(1:4000, 2L), status = "ok", price_pct = c(1:4000, -(1:4000)))
  failed <- rows$regime == "merger" & rows$draw_index == 1002L
  rows$status[failed] <- "failed"; rows$price_pct[failed] <- NA_real_
  rows <- rows[nrow(rows):1L, ]
  merger <- runner$nutpie_uncertainty_values(list(list(draws = rows)), "merger")[[1L]]
  tariff <- runner$nutpie_uncertainty_values(list(list(draws = rows)), "tariff25")[[1L]]
  expect_identical(dim(merger), c(1000L, 4L))
  expect_identical(dim(tariff), c(1000L, 4L))
  expect_true(is.na(merger[2L, 2L]))
  expect_true(is.na(tariff[2L, 2L]))
  expect_equal(merger[500L, 4L], 3500)
  expect_equal(tariff[500L, 4L], -3500)
  duplicate <- rbind(rows, rows[1L, ])
  expect_error(runner$nutpie_uncertainty_values(list(list(draws = duplicate)), "merger"), "duplicate")
  incomplete <- rows[!(rows$regime == "tariff25" & rows$draw_index == 2003L), ]
  expect_error(runner$nutpie_uncertainty_values(list(list(draws = incomplete)), "merger"), "Incomplete")
  invalid <- rows
  invalid$price_pct[which(invalid$status == "ok")[1L]] <- NA_real_
  expect_error(runner$nutpie_uncertainty_values(list(list(draws = invalid)), "merger"), "successful saved policy row")
})

test_that("uncertainty audit rejects result files from another producer manifest", {
  runner <- closeout_uncertainty_runner()
  hash <- paste(rep("a", 64L), collapse = "")
  results <- rep(list(list(status = "completed", manifest_hash = hash)), 12L)
  expect_true(runner$nutpie_uncertainty_input_contract(results, hash))
  results[[7L]]$manifest_hash <- paste(rep("b", 64L), collapse = "")
  expect_error(runner$nutpie_uncertainty_input_contract(results, hash), "producer-run manifest")
})
