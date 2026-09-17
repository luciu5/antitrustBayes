test_that("only validated forward kernel migrations have a coordinate policy", {
  hashes <- lapply(.bayes_bridge_approved_hashes, .bayes_bridge_read_hashes)
  for (kernel in .bayes_bridge_kernel_names) {
    original <- hashes$original[[kernel]]
    rotation <- hashes$rotation[[kernel]]
    helmert <- hashes$helmert[[kernel]]
    geometry <- function(a, b, k = kernel) .bayes_bridge_migration_geometry(
      a, b, k, hashes$original, hashes$rotation, hashes$helmert
    )
    expect_identical(geometry(original, rotation), "legacy")
    expect_identical(geometry(original, helmert), "legacy")
    expect_identical(geometry(rotation, helmert), "preserve")
    expect_null(geometry(helmert, rotation))
    expect_null(geometry(rotation, original))
    expect_null(geometry(helmert, helmert))
    expect_null(geometry(rotation, paste(rep("0", 64L), collapse = "")))
    expect_null(geometry(rotation, helmert, "wrong_kernel"))
    expect_true(.bayes_bridge_is_approved_migration(
      rotation, helmert, kernel, hashes$original, hashes$rotation, hashes$helmert
    ))
    for (other in setdiff(.bayes_bridge_kernel_names, kernel)) {
      expect_null(geometry(rotation, hashes$helmert[[other]]))
    }
  }
})

if (!methods::isClass("bridge_private_data_fixture")) {
  methods::setClass("bridge_private_data_fixture", slots = c(stan_data = "list"))
}

test_that("performance migration preserves legacy and rotated private coordinates", {
  for (flag in c(1L, 2L)) {
    fit <- methods::new("bridge_private_data_fixture",
      stan_data = list(market_s0_parameterization = flag, other = 7L))
    before <- serialize(fit, NULL)
    data <- .bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "preserve"))
    expect_identical(data, fit@stan_data)
    expect_identical(serialize(fit, NULL), before)
  }
  fit <- methods::new("bridge_private_data_fixture", stan_data = list(other = 7L))
  before <- serialize(fit, NULL)
  data <- .bayes_bridge_private_data(fit,
    list(migration = TRUE, migration_geometry = "legacy"))
  expect_identical(data$market_s0_parameterization, 1L)
  expect_identical(serialize(fit, NULL), before)
  for (bad in list(NULL, 0L, 3L, 1.5, NA_real_, Inf, c(1L, 2L), "1")) {
    fit@stan_data <- list(market_s0_parameterization = bad)
    expect_error(.bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "preserve")), "saved private data flag")
  }
  for (bad in list(0L, 2L, 1.5, NA_real_, Inf, c(1L, 2L), "1")) {
    fit@stan_data <- list(market_s0_parameterization = bad)
    expect_error(.bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "legacy")), "legacy private data flag")
  }
  expect_error(.bayes_bridge_private_data(fit, list(migration = TRUE)),
    "no validated coordinate policy")
})
