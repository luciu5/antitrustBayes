test_that("bridge migration accepts only the documented kernel hash pairs", {
  original <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$original)
  rotation <- .bayes_bridge_read_hashes(.bayes_bridge_approved_hashes$rotation)
  expect_named(original, .bayes_bridge_kernel_names)
  expect_named(rotation, .bayes_bridge_kernel_names)

  for (kernel in .bayes_bridge_kernel_names) {
    expect_true(.bayes_bridge_is_approved_migration(
      original[[kernel]], rotation[[kernel]], kernel, original, rotation
    ))
    expect_false(.bayes_bridge_is_approved_migration(
      original[[kernel]], rotation[[kernel]], "wrong_kernel", original, rotation
    ))
    altered <- rotation
    altered[[kernel]] <- paste0("0", substr(altered[[kernel]], 2L, 64L))
    expect_false(.bayes_bridge_is_approved_migration(
      original[[kernel]], rotation[[kernel]], kernel, original, altered
    ))
  }
})

test_that("private bridge data canonicalizes optional duplicate values safely", {
  got <- .bayes_bridge_canonicalize_data(list(
    cutoff_scale = NULL,
    cutoff_scale = 0.2,
    case_unit = -1L,
    case_unit = -1L,
    other = 7
  ))
  expect_identical(got$cutoff_scale, 0.2)
  expect_identical(got$case_unit, -1L)
  expect_identical(got$other, 7)
  expect_error(
    .bayes_bridge_canonicalize_data(list(cutoff_scale = 0.2,
                                         cutoff_scale = 0.3)),
    "conflicting duplicate field 'cutoff_scale'"
  )
  expect_error(
    .bayes_bridge_canonicalize_data(list(case_unit = -1L, case_unit = 1L)),
    "conflicting duplicate field 'case_unit'"
  )
  expect_equal(
    .bayes_bridge_canonicalize_data(list(x = 1, x = 1 + 1e-13))$x,
    1
  )
})

test_that("bridge context construction preserves the caller RNG state", {
  set.seed(1729)
  before <- get(".Random.seed", envir = .GlobalEnv)
  value <- .bayes_bridge_preserve_rng({
    set.seed(271828)
    stats::rnorm(4L)
  })
  expect_length(value, 4L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
})

if (!methods::isClass("bridge_context_mock_stanfit")) {
  methods::setClass("bridge_context_mock_stanfit",
                    slots = c(sim = "list", par_dims = "list",
                              inits = "list", model_pars = "character",
                              stan_code = "character"))
}

.bridge_context_mock_stanfit <- function() {
  chain <- list(
    alpha = seq_len(10L) / 10,
    beta = seq_len(10L),
    `tf_generated[1]` = seq_len(10L),
    `tf_generated[2]` = seq_len(10L) + 100,
    `lp__` = seq_len(10L) * -0.5
  )
  methods::new(
    "bridge_context_mock_stanfit",
    sim = list(
      pars_oi = c("alpha", "beta", "tf_generated", "lp__"),
      fnames_oi = c("alpha", "beta", "tf_generated[1]",
                    "tf_generated[2]", "lp__"),
      samples = list(chain, chain)
    ),
    par_dims = list(alpha = numeric(0), beta = 1, `lp__` = numeric(0),
                    tf_generated = 2),
    inits = list(list(alpha = 0.1, beta = 1, tf_generated = c(0, 0))),
    model_pars = c("alpha", "beta", "tf_generated", "lp__"),
    stan_code = paste(
      "parameters {",
      "real alpha;",
      "vector[1] beta;",
      "}",
      "transformed parameters {",
      "vector[2] tf_generated;",
      "}",
      sep = "\n"
    )
  )
}

test_that("positive initialization schemas are exact while zero-width fields may differ", {
  saved <- .bridge_context_mock_stanfit()
  expect_identical(.bayes_bridge_raw_parameter_names(saved), c("alpha", "beta"))
  current <- saved
  current@par_dims$c_s0 <- 0L
  current@inits[[1L]]$c_s0 <- numeric(0L)
  expect_invisible(.bayes_bridge_assert_positive_schema(saved, current))

  current <- saved
  current@par_dims$tf_generated <- 3L
  current@inits[[1L]]$tf_generated <- c(0, 0, 0)
  expect_invisible(.bayes_bridge_assert_positive_schema(saved, current))

  current <- saved
  current@inits[[1L]]$beta <- NULL
  expect_error(
    .bayes_bridge_assert_positive_schema(saved, current),
    "initialization is missing a declared raw parameter"
  )

  current <- saved
  current@par_dims$beta <- 2L
  expect_error(
    .bayes_bridge_assert_positive_schema(saved, current),
    "different nonzero parameter names or dimensions"
  )

  current <- saved
  current@inits[[1L]]$gamma <- 1
  current@par_dims$gamma <- numeric(0)
  current@model_pars <- c("alpha", "gamma", "tf_generated", "lp__")
  current@stan_code <- paste(
    "parameters {", "real alpha;", "real gamma;", "}",
    "transformed parameters {", "vector[2] tf_generated;", "}",
    sep = "\n"
  )
  expect_error(
    .bayes_bridge_assert_positive_schema(saved, current),
    "different nonzero parameter names or dimensions"
  )

  current <- saved
  current@inits[[1L]]$c_s0 <- 0.2
  current@par_dims$c_s0 <- 1L
  current@model_pars <- c("alpha", "beta", "c_s0", "tf_generated", "lp__")
  current@stan_code <- paste(
    "parameters {", "real alpha;", "real beta;", "real c_s0;", "}",
    "transformed parameters {", "vector[2] tf_generated;", "}",
    sep = "\n"
  )
  expect_error(
    .bayes_bridge_assert_positive_schema(saved, current),
    "different nonzero parameter names or dimensions"
  )
})

test_that("saved bridge rows reconstruct the original constrained schema", {
  st <- .bridge_context_mock_stanfit()
  rows <- .bayes_bridge_saved_rows(st, n_rows = 5L)
  expect_equal(nrow(rows), 10L)
  expect_equal(rows$iteration[1:5], c(1L, 3L, 6L, 8L, 10L))
  theta <- .bayes_bridge_theta_at(st, chain = 2L, iteration = 6L)
  expect_named(theta, c("alpha", "beta", "tf_generated", "lp__"))
  expect_equal(theta$alpha, 0.6)
  expect_equal(as.numeric(theta$beta), 6)
  expect_equal(theta$lp__, -3)

  st@sim$warmup2 <- c(5, 5)
  post_rows <- .bayes_bridge_saved_rows(st, n_rows = 5L)
  expect_equal(post_rows$iteration[1:5], 6:10)
})
