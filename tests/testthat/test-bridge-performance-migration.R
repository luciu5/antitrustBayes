testthat::test_that("LOO and bridge evidence reject non-HMC inference", {
  kernel <- "ces_supply"
  source_file <- .cmdstan_kernel_file(kernel)
  source_hash <- digest::digest(file = source_file, algo = "sha256")
  source_code <- paste0(paste(readLines(source_file, warn = FALSE),
                              collapse = "\n"), "\n")
  metadata <- list(
    backend = "cmdstanr", inference = "pathfinder", stan_kernel = kernel,
    stan_source_hash = source_hash
  )

  testthat::expect_error(
    .bayes_cmdstan_source_contract(metadata, source_code),
    "require HMC posterior draws"
  )
  testthat::expect_identical(
    .bayes_cmdstan_source_contract(metadata, source_code, require_hmc = FALSE)$inference,
    "pathfinder"
  )
  metadata$backend <- "rstan"
  testthat::expect_error(
    .bayes_cmdstan_source_contract(metadata, source_code, require_hmc = FALSE),
    "CmdStanR backend"
  )
})
