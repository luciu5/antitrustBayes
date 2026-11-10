test_that("line-ending portability preserves raw seals and rejects source edits", {
  kernel <- "logit_supply"
  file <- .cmdstan_kernel_file(kernel)
  raw <- readChar(file, file.info(file)$size, useBytes = TRUE)
  lf <- gsub("\r\n", "\n", raw, fixed = TRUE)
  crlf <- gsub("\n", "\r\n", lf, fixed = TRUE)
  metadata <- list(backend = "cmdstanr", inference = "hmc",
                   stan_kernel = kernel)
  hash <- function(text) digest::digest(text, algo = "sha256", serialize = FALSE)
  identity <- .bayes_stan_source_identity(file)
  expect_identical(identity, list(version = "lf-v1", sha256 = hash(lf)))

  for (source in list(lf, crlf)) {
    metadata$stan_source_hash <- hash(source)
    sealed <- serialize(metadata, NULL)
    for (versioned in c(FALSE, TRUE)) {
      checked <- metadata
      if (versioned) checked$stan_source_identity <- identity
      contract <- .bayes_cmdstan_source_contract(checked, source)
      expect_identical(contract$saved_hash, hash(source))
      expect_identical(contract$current_hash,
                       digest::digest(file = file, algo = "sha256"))
      expect_identical(contract$source_identity, identity)
      expect_error(.bayes_cmdstan_source_contract(checked,
        sub("normal", "student_t", source, fixed = TRUE)), "does not match")
      altered <- checked
      altered$stan_source_hash <- hash(paste0(source, "// edited\n"))
      expect_error(.bayes_cmdstan_source_contract(altered, source),
                   "source hashes differ")
    }
    expect_identical(serialize(metadata, NULL), sealed)
  }
  metadata$stan_source_identity <- list(version = "lf-v2", sha256 = hash(lf))
  expect_error(.bayes_cmdstan_source_contract(metadata, lf),
               "canonical Stan source identities differ")
  metadata$stan_source_identity <- list(version = "lf-v1", sha256 = strrep("0", 64))
  expect_error(.bayes_cmdstan_source_contract(metadata, lf),
               "canonical Stan source identities differ")
})

test_that("package provenance cannot inherit an analysis repository commit", {
  local_mocked_bindings(.bayes_package_description = function() {
    list(Version = "0.3.0.9001", RemoteSha = strrep("a", 40))
  })
  withr::local_envvar(ANTITRUSTBAYES_COMMIT = strrep("b", 40))
  provenance <- package_provenance(data.frame(x = 1), "baseline")
  expect_identical(provenance$commit, strrep("a", 40))
  expect_identical(provenance$package_version, "0.3.0.9001")
  local_mocked_bindings(.bayes_package_description = function() list(Version = "1.0"))
  expect_identical(package_provenance(NULL, "baseline")$commit, NA_character_)
  local_mocked_bindings(.bayes_package_description = function() {
    list(Version = "1.0", RemoteSha = "main")
  })
  expect_identical(package_provenance(NULL, "baseline")$commit, NA_character_)
})
