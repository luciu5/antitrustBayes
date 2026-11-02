.serialization_runner_path <- testthat::test_path(
  "..", "..", "doc", "validation", "nutpie_serialization_validation.R")
testthat::skip_if_not(file.exists(.serialization_runner_path),
  "Saved-fit research audit requires the repository doc/ tree.")

.serialization_runner <- function() {
  environment <- new.env(parent = globalenv())
  source(.serialization_runner_path, local = environment)
  environment
}

test_that("optional-free workers exclude entire native package libraries", {
  environment <- .serialization_runner()
  root <- tempfile("serialization-libraries-")
  directories <- file.path(root, c("compatible", "native", "another native"))
  for (path in directories) dir.create(path, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(directories[[2L]], "nutpieR"))
  file.create(file.path(directories[[2L]], "nutpieR", "DESCRIPTION"))
  dir.create(file.path(directories[[3L]], "bridgestan"))
  file.create(file.path(directories[[3L]], "bridgestan", "DESCRIPTION"))
  actual <- environment$nutpie_serialization_libraries(c(directories, directories[[1L]], "absent-library"))
  expect_identical(actual, normalizePath(directories[[1L]], mustWork = TRUE))
})

test_that("optional-free audit rejects both loaded and discoverable native packages", {
  environment <- .serialization_runner()
  environment$loadedNamespaces <- function() character()
  environment$requireNamespace <- function(package, quietly) FALSE
  expect_invisible(environment$nutpie_serialization_no_native())
  environment$loadedNamespaces <- function() "nutpieR"
  expect_error(environment$nutpie_serialization_no_native(), "unavailable: nutpieR")
  environment$loadedNamespaces <- function() character()
  environment$requireNamespace <- function(package, quietly) package == "bridgestan"
  expect_error(environment$nutpie_serialization_no_native(), "unavailable: bridgestan")
})

test_that("saved-fit audit requires explicit execution and correct rejection reasons", {
  environment <- .serialization_runner()
  old_execute <- Sys.getenv("AB_NUTPIE_SERIALIZATION_EXECUTE", unset = NA_character_)
  on.exit(if (is.na(old_execute)) Sys.unsetenv("AB_NUTPIE_SERIALIZATION_EXECUTE") else
    Sys.setenv(AB_NUTPIE_SERIALIZATION_EXECUTE = old_execute), add = TRUE)
  Sys.setenv(AB_NUTPIE_SERIALIZATION_EXECUTE = "false")
  expect_error(environment$nutpie_serialization_main("not-a-repository"), "requires AB_NUTPIE_SERIALIZATION_EXECUTE=true")
  expect_identical(environment$nutpie_serialization_expect_error(stop("binary hash differs"), "hash differs"),
    "binary hash differs")
  expect_error(environment$nutpie_serialization_expect_error(1, "hash differs"), "did not occur")
  expect_error(environment$nutpie_serialization_expect_error(stop("wrong failure"), "hash differs"), "Unexpected artifact rejection")
})
