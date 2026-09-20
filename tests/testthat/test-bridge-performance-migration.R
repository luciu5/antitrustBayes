test_that("only validated forward kernel migrations have a coordinate policy", {
  hashes <- lapply(.bayes_bridge_approved_hashes, .bayes_bridge_read_hashes)
  for (kernel in .bayes_bridge_kernel_names) {
    original <- hashes$original[[kernel]]
    rotation <- hashes$rotation[[kernel]]
    helmert <- hashes$helmert[[kernel]]
    joint <- hashes$joint[[kernel]]
    scale <- hashes$scale[[kernel]]
    margin_first <- hashes$margin_first[[kernel]]
    geometry <- function(a, b, k = kernel) .bayes_bridge_migration_geometry(
      a, b, k, hashes$original, hashes$rotation, hashes$helmert, hashes$joint,
      hashes$scale, hashes$margin_first
    )
    expect_identical(geometry(original, rotation), "legacy")
    expect_identical(geometry(original, helmert), "legacy")
    expect_identical(geometry(rotation, helmert), "preserve")
    expect_identical(geometry(original, joint), "legacy")
    expect_identical(geometry(rotation, joint), "preserve")
    expect_identical(geometry(helmert, joint), "preserve")
    expect_identical(geometry(original, scale), "legacy")
    expect_identical(geometry(rotation, scale), "preserve")
    expect_identical(geometry(helmert, scale), "preserve")
    expect_identical(geometry(joint, scale), "preserve")
    expect_identical(geometry(original, margin_first), "legacy")
    expect_identical(geometry(rotation, margin_first), "preserve")
    expect_identical(geometry(helmert, margin_first), "preserve")
    expect_identical(geometry(joint, margin_first), "preserve")
    expect_identical(geometry(scale, margin_first), "preserve")
    expect_null(geometry(margin_first, scale))
    expect_null(geometry(scale, joint))
    expect_null(geometry(joint, helmert))
    expect_null(geometry(helmert, rotation))
    expect_null(geometry(rotation, original))
    expect_null(geometry(helmert, helmert))
    expect_null(geometry(rotation, paste(rep("0", 64L), collapse = "")))
    expect_null(geometry(rotation, helmert, "wrong_kernel"))
    expect_identical(.bayes_bridge_migration_scale_policy(
      scale, margin_first, kernel, hashes$scale, hashes$margin_first
    ), "preserve")
    expect_identical(.bayes_bridge_migration_scale_policy(
      joint, margin_first, kernel, hashes$scale, hashes$margin_first
    ), "legacy")
    expect_true(.bayes_bridge_is_approved_migration(
      rotation, helmert, kernel, hashes$original, hashes$rotation,
      hashes$helmert, hashes$joint, hashes$scale, hashes$margin_first
    ))
    for (other in setdiff(.bayes_bridge_kernel_names, kernel)) {
      expect_null(geometry(rotation, hashes$helmert[[other]]))
    }
  }
})

test_that("margin-first migration hashes identify the installed Stan kernels", {
  models <- get("stanmodels", envir = asNamespace("antitrustBayes"))
  hashes <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$margin_first)
  for (kernel in .bayes_bridge_kernel_names) {
    current <- digest::digest(models[[kernel]]@model_code,
                              algo = "sha256", serialize = FALSE)
    expect_identical(current, unname(hashes[[kernel]]), info = kernel)
  }
})

if (!methods::isClass("bridge_private_data_fixture")) {
  methods::setClass("bridge_private_data_fixture", slots = c(stan_data = "list"))
}

test_that("performance migration preserves legacy and rotated private coordinates", {
  for (flag in c(1L, 2L, 3L)) {
    fit <- methods::new("bridge_private_data_fixture",
      stan_data = list(market_s0_parameterization = flag, other = 7L))
    before <- serialize(fit, NULL)
    data <- .bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "preserve"))
    expect_identical(data$market_s0_parameterization, flag)
    expect_identical(data$quality_scale_parameterization, 1L)
    expect_identical(data$other, 7L)
    expect_identical(serialize(fit, NULL), before)
  }
  fit <- methods::new("bridge_private_data_fixture", stan_data = list(other = 7L))
  before <- serialize(fit, NULL)
  data <- .bayes_bridge_private_data(fit,
    list(migration = TRUE, migration_geometry = "legacy"))
  expect_identical(data$market_s0_parameterization, 1L)
  expect_identical(data$quality_scale_parameterization, 1L)
  expect_identical(serialize(fit, NULL), before)
  for (bad in list(NULL, 0L, 4L, 1.5, NA_real_, Inf, c(1L, 2L), "1")) {
    fit@stan_data <- list(market_s0_parameterization = bad)
    expect_error(.bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "preserve")), "saved private data flag")
  }
  fit@stan_data <- list(market_s0_parameterization = 3L,
                        quality_scale_parameterization = 2L)
  expect_error(.bayes_bridge_private_data(fit,
    list(migration = TRUE, migration_geometry = "preserve")),
    "invalid quality_scale_parameterization")
  for (scale_flag in c(1L, 2L)) {
    fit@stan_data <- list(market_s0_parameterization = 3L,
                          quality_scale_parameterization = scale_flag)
    data <- .bayes_bridge_private_data(fit, list(
      migration = TRUE, migration_geometry = "preserve",
      migration_scale_policy = "preserve"
    ))
    expect_identical(data$quality_scale_parameterization, scale_flag)
  }
  fit@stan_data <- list(market_s0_parameterization = 3L)
  expect_error(.bayes_bridge_private_data(fit, list(
    migration = TRUE, migration_geometry = "preserve",
    migration_scale_policy = "preserve"
  )), "missing its saved")
  for (scale_flag in c(0L, 3L, 1.5, NA_real_, Inf, c(1L, 2L), "1")) {
    fit@stan_data <- list(market_s0_parameterization = 3L,
                          quality_scale_parameterization = scale_flag)
    expect_error(.bayes_bridge_private_data(fit, list(
      migration = TRUE, migration_geometry = "preserve",
      migration_scale_policy = "preserve"
    )), "invalid quality_scale_parameterization")
  }
  for (bad in list(0L, 2L, 1.5, NA_real_, Inf, c(1L, 2L), "1")) {
    fit@stan_data <- list(market_s0_parameterization = bad)
    expect_error(.bayes_bridge_private_data(fit,
      list(migration = TRUE, migration_geometry = "legacy")), "legacy private data flag")
  }
  expect_error(.bayes_bridge_private_data(fit, list(migration = TRUE)),
    "no validated coordinate policy")
})
