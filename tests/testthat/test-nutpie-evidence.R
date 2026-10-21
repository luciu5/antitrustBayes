make_nutpie_evidence_fixture <- function() {
  kernel <- "logit_supply"
  source_file <- .cmdstan_kernel_file(kernel)
  source_hash <- digest::digest(file = source_file, algo = "sha256")
  source_code <- paste0(paste(readLines(source_file, warn = FALSE),
                              collapse = "\n"), "\n")
  values <- array(-seq_len(24L) / 10, c(4L, 2L, 3L), dimnames = list(
    iteration = as.character(1:4), chain = as.character(1:2),
    variable = paste0("log_lik[", c(3L, 1L, 2L), "]")))
  native <- structure(list(draws_array = values, metadata = list(
    source_code = source_code, source_hash = source_hash,
    inference = "hmc")), class = "NutpieFit")
  methods::new("BayesLogit", stanfit = native, metadata = list(
    backend = "nutpie", inference = "hmc", stan_kernel = kernel,
    stan_source_hash = source_hash))
}

test_that("Nutpie predictive units preserve chain identity and explicit indexing", {
  fit <- make_nutpie_evidence_fixture()
  got <- .bayes_log_lik_array(fit@stanfit, "log_lik", expected_units = 3L)
  expect_identical(dim(got), c(4L, 2L, 3L))
  expect_identical(dimnames(got)[[3L]], paste0("log_lik[", 1:3, "]"))
  expect_equal(got[, 2L, 1L], fit@stanfit$draws_array[, 2L, 2L])
  expect_error(.bayes_log_lik_array(fit@stanfit, "log_lik", 4L),
               "expected predictive-unit count")
  fit@stanfit$draws_array[1L, 1L, 1L] <- -Inf
  expect_error(.bayes_log_lik_array(fit@stanfit, "log_lik"), "non-finite")
})

test_that("Nutpie source provenance uses the exact existing kernel and HMC", {
  fit <- make_nutpie_evidence_fixture()
  expect_identical(.bayes_evidence_contract(fit)$backend, "nutpie")
  changed <- fit
  changed@stanfit$metadata$source_hash <- strrep("0", 64L)
  expect_error(.bayes_evidence_contract(changed), "source hashes differ")
  changed <- fit
  changed@stanfit$metadata$source_code <- paste0(
    changed@stanfit$metadata$source_code, "// changed\n")
  expect_error(.bayes_evidence_contract(changed), "does not match")
  changed <- fit
  changed@metadata$inference <- "pathfinder"
  expect_error(.bayes_evidence_contract(changed), "require HMC")
})

test_that("backend choice does not restrict evidence and metadata cannot bypass target contracts", {
  fit <- make_nutpie_evidence_fixture()
  fit@metadata$bridge_validated <- TRUE
  fit@stanfit$metadata$bridge_validated <- TRUE
  expect_identical(.bayes_evidence_contract(fit)$backend, "nutpie")
  changed <- fit
  changed@stanfit$metadata$source_hash <- strrep("0", 64L)
  expect_error(.bayes_evidence_contract(changed), "source hashes differ")
  cmdstan <- fit
  cmdstan@stanfit <- NULL
  cmdstan@metadata$backend <- "cmdstanr"
  expect_identical(unname(.bayes_comparison_backends(list(fit, cmdstan))),
                   c("nutpie", "cmdstanr"))
  expect_identical(unname(.bayes_comparison_backends(list(fit, cmdstan),
    require_evidence = FALSE)), c("nutpie", "cmdstanr"))
  expect_identical(unname(.bayes_comparison_backends(list(fit, fit))),
                   c("nutpie", "nutpie"))
})

test_that("normalized bridge coordinates exclude generated quantities and include one Jacobian", {
  fit <- make_nutpie_evidence_fixture()
  q <- cbind(c(-1, 0, 1, .5), c(-2, 0, 2, -.5))
  theta <- cbind(exp(q[, 1L]), plogis(q[, 2L]))
  ## Retain reversed columns and an unrelated GQ: the model's parameter names
  ## determine the coordinate order, never the position of posterior columns.
  fit@stanfit$draws_array <- array(c(theta[, 2L], rep(999, 4L), theta[, 1L]),
    c(4L, 1L, 3L), dimnames = list(NULL, NULL, c("share", "log_lik[1]", "scale")))
  calls <- new.env(parent = emptyenv())
  calls$flags <- list()
  model <- list(
    param_names = function(include_tp, include_gq) c("scale", "share"),
    param_unc_names = function() c("scale", "share"),
    param_unc_num = function() 2L,
    param_unconstrain = function(theta) c(log(theta[[1L]]), qlogis(theta[[2L]])),
    param_constrain = function(q, include_tp, include_gq)
      c(exp(q[[1L]]), plogis(q[[2L]])),
    log_density = function(q, propto, jacobian) {
      calls$flags[[length(calls$flags) + 1L]] <- c(propto, jacobian)
      ## Exp(1) scale prior and Uniform(0,1) share prior, on q coordinates.
      -exp(q[[1L]]) + q[[1L]] +
        log(plogis(q[[2L]])) + log1p(-plogis(q[[2L]]))
    }
  )
  context <- .bayes_nutpie_bridge_context_from_model(
    fit, model, list(kernel = "logit_supply", saved_hash = "same", current_hash = "same"))
  expect_equal(unname(context$samples), unname(q))
  expect_identical(context$constrained_parameter_names, c("scale", "share"))
  expect_identical(context$bridge_args$lb, c(scale = -Inf, share = -Inf))
  expect_identical(context$bridge_args$ub, c(scale = Inf, share = Inf))
  expect_equal(context$log_posterior(c(0, 0), NULL), -1 - log(4))
  expect_true(all(vapply(calls$flags, identical, logical(1), c(FALSE, TRUE))))
  model$param_constrain <- function(q, include_tp, include_gq) c(2, 2)
  expect_error(.bayes_nutpie_bridge_context_from_model(fit, model, list()),
               "round trip failed")
})

test_that("bridge restoration requires immutable Stan data provenance", {
  fit <- make_nutpie_evidence_fixture()
  fit@stan_data <- list(N = 2L, ignored_metadata = "not Stan data")
  fit@metadata$stan_data_fields <- "N"
  fit@metadata$stan_data_hash <- digest::digest(list(N = 2L), algo = "sha256")
  expect_identical(.bayes_bridge_data(fit), list(N = 2L))
  fit@stan_data$N <- 3L
  expect_error(.bayes_bridge_data(fit), "differ from the data used")
  fit@metadata$stan_data_fields <- c("N", "N")
  expect_error(.bayes_bridge_data(fit), "field contract")
})

test_that("one-dimensional bridge coordinates retain each iteration and chain", {
  fit <- make_nutpie_evidence_fixture()
  q <- seq(-1, 1, length.out = 6L)
  fit@stanfit$draws_array <- array(exp(q), c(3L, 2L, 1L),
    dimnames = list(NULL, NULL, "scale"))
  model <- list(
    param_names = function(include_tp, include_gq) "scale",
    param_unc_names = function() "scale",
    param_unc_num = function() 1L,
    param_unconstrain = log,
    param_constrain = function(q, include_tp, include_gq) exp(q),
    log_density = function(q, propto, jacobian) -exp(q) + q
  )
  context <- .bayes_nutpie_bridge_context_from_model(fit, model, list())
  expect_identical(dim(context$samples), c(6L, 1L))
  expect_equal(as.numeric(context$samples), q)
  expect_identical(context$draw_index$iteration, rep(1:3, 2L))
  expect_identical(context$draw_index$chain, rep(1:2, each = 3L))
})
