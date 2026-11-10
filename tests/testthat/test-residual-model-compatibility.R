# Material failure: changing the statistical default can invalidate sealed old
# posterior/evidence artifacts or evaluate them under a different Stan target.
test_that("specification 14 retains its statistical meaning without mutation", {
  data <- list(N = 3L)
  attr(data, "model_spec_version") <- 14L
  fit <- methods::new("BayesLogit", stan_data = data,
    metadata = list(model_spec_version = 14L))
  sealed <- serialize(fit, NULL)
  expect_silent(validate_fit_spec(fit))
  expect_identical(serialize(fit, NULL), sealed)
  fit@stan_data$share_residual_model <- 2L
  expect_error(validate_fit_spec(fit), "cannot contain")
  fit@stan_data$share_residual_model <- NULL
  fit@metadata$model_spec_version <- 15L
  expect_error(validate_fit_spec(fit), "metadata disagree")
  fit@metadata$model_spec_version <- 13L
  attr(fit@stan_data, "model_spec_version") <- 13L
  expect_error(validate_fit_spec(fit), "Refit")
})

test_that("legacy evidence selects only its sealed original kernel", {
  for (kernel in .bayes_bridge_kernel_names) {
    file <- .cmdstan_kernel_file(kernel, 14L)
    source <- paste0(paste(readLines(file, warn = FALSE), collapse = "\n"), "\n")
    metadata <- list(backend = "cmdstanr", inference = "hmc", stan_kernel = kernel,
      model_spec_version = 14L, stan_source_hash = digest::digest(file = file, algo = "sha256"),
      stan_source_identity = .bayes_stan_source_identity(file))
    sealed <- serialize(metadata, NULL)
    contract <- .bayes_cmdstan_source_contract(metadata, source)
    expect_identical(contract$source_file, normalizePath(file))
    expect_identical(contract$current_hash, metadata$stan_source_hash)
    expect_identical(serialize(metadata, NULL), sealed)
    expect_error(.bayes_cmdstan_source_contract(metadata, paste0(source, "// changed\n")), "does not match")
    metadata$model_spec_version <- 15L
    expect_error(.bayes_cmdstan_source_contract(metadata, source), "source")
  }
})
