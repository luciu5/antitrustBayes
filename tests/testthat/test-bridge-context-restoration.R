testthat::test_that("CmdStan bridge evidence requires the exact current source", {
  kernel <- "logit_supply"
  source_file <- .cmdstan_kernel_file(kernel)
  source_hash <- digest::digest(file = source_file, algo = "sha256")
  source_code <- paste0(paste(readLines(source_file, warn = FALSE),
                              collapse = "\n"), "\n")
  metadata <- list(
    backend = "cmdstanr", inference = "hmc", stan_kernel = kernel,
    stan_source_hash = source_hash
  )

  contract <- .bayes_cmdstan_source_contract(metadata, source_code)
  testthat::expect_identical(contract$kernel, kernel)
  testthat::expect_identical(contract$saved_hash, source_hash)
  testthat::expect_identical(contract$current_hash, source_hash)
  testthat::expect_false(contract$migration %||% FALSE)

  stale_metadata <- metadata
  replacement <- if (substr(source_hash, 1L, 1L) == "0") "1" else "0"
  stale_metadata$stan_source_hash <- paste0(replacement, substr(source_hash, 2L, 64L))
  testthat::expect_error(
    .bayes_cmdstan_source_contract(stale_metadata, source_code),
    "source hashes differ"
  )
  changed_code <- paste0(source_code, "// changed source\n")
  testthat::expect_error(
    .bayes_cmdstan_source_contract(metadata, changed_code),
    "does not match the current kernel source"
  )
})

testthat::test_that("LOO log-likelihood extraction preserves chain and unit axes", {
  values <- array(seq_len(24L) / -10, dim = c(4L, 2L, 3L),
                  dimnames = list(
                    iteration = as.character(seq_len(4L)),
                    chain = as.character(seq_len(2L)),
                    variable = paste0("log_lik[", 1:3, "]")
                  ))
  fit <- new.env(parent = emptyenv())
  fit$draws <- function(variables, format) values
  class(fit) <- c("CmdStanMCMC", "CmdStanFit", "R6")

  got <- .bayes_cmdstan_log_lik_array(fit, "log_lik", expected_units = 3L)
  testthat::expect_identical(dim(got), c(4L, 2L, 3L))
  testthat::expect_identical(got, values)
  testthat::expect_error(
    .bayes_cmdstan_log_lik_array(fit, "log_lik", expected_units = 4L),
    "expected predictive-unit count"
  )

  fit$draws <- function(variables, format) values[, , c(1L, 3L), drop = FALSE]
  testthat::expect_error(
    .bayes_cmdstan_log_lik_array(fit, "log_lik"),
    "noncontiguous unit indices"
  )
})
