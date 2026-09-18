## These tests exercise only the post-sampling existence contract.  They use
## hand-built stanfit objects, so no Stan model is compiled or sampled.

.mock_sampling_chain <- function(n, pars = c("theta[1]", "lp__")) {
  stats::setNames(lapply(pars, function(par) rep(0, n)), pars)
}

.mock_sampling_stanfit <- function(mode = 0L,
                                   sim_chains = 2L,
                                   n_save = c(3L, 3L),
                                   warmup2 = c(1L, 1L),
                                   samples = NULL,
                                   sim_extra = list()) {
  if (is.null(samples)) {
    n_for_chain <- suppressWarnings(as.integer(n_save))
    n_for_chain[!is.finite(n_for_chain) | n_for_chain < 0L] <- 0L
    samples <- lapply(n_for_chain, .mock_sampling_chain)
  }
  sim <- c(
    list(chains = sim_chains, n_save = n_save, warmup2 = warmup2,
         samples = samples),
    sim_extra
  )
  methods::new(
    "stanfit",
    model_name = "sampling-contract-mock",
    model_pars = "lp__",
    par_dims = list(lp__ = 1L),
    mode = mode,
    sim = sim,
    inits = list(),
    stan_args = list(),
    stanmodel = methods::new("stanmodel"),
    date = "",
    .MISC = new.env(parent = emptyenv())
  )
}

test_that("the sampling contract rejects empty and partial fits", {
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(mode = 2L), requested_chains = 2L
    ),
    "mode must be 0"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(sim_chains = 1L), requested_chains = 2L
    ),
    "returned chain count"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(n_save = NULL), requested_chains = 2L
    ),
    "n_save and warmup2"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(n_save = c(2L, 2L), warmup2 = c(2L, 1L)),
      requested_chains = 2L
    ),
    "post-warmup"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(samples = list(.mock_sampling_chain(3L))),
      requested_chains = 2L
    ),
    "one or more chains"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(
        samples = list(.mock_sampling_chain(0L), .mock_sampling_chain(3L))
      ),
      requested_chains = 2L
    ),
    "empty"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(
        n_save = c(2.5, 3L), warmup2 = c(1L, 1L),
        samples = list(.mock_sampling_chain(2L), .mock_sampling_chain(3L))
      ),
      requested_chains = 2L
    ),
    "n_save and warmup2"
  )
  expect_error(
    .assert_stan_sampling_contract(
      .mock_sampling_stanfit(
        n_save = c(3L, 3L), warmup2 = c(1L, 1L),
        samples = list(.mock_sampling_chain(3L), .mock_sampling_chain(1L))
      ),
      requested_chains = 2L
    ),
    "do not match"
  )
})

test_that("the sampling contract rejects incomplete parameter vectors", {
  chain <- .mock_sampling_chain(3L)
  chain[["theta[1]"]] <- numeric(2L)
  expect_error(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(samples = list(chain, .mock_sampling_chain(3L))),
    requested_chains = 2L), "do not match")

  expect_error(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(samples = list(.mock_sampling_chain(3L),
      .mock_sampling_chain(3L, "lp__"))), requested_chains = 2L),
    "parameter names differ")

  expect_error(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(samples = list(list(), .mock_sampling_chain(3L))),
    requested_chains = 2L), "empty")

  chain <- .mock_sampling_chain(3L)
  names(chain) <- NULL
  expect_error(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(samples = list(chain, .mock_sampling_chain(3L))),
    requested_chains = 2L), "do not match")

  ## A sampler-parameter matrix is not a real @sim$samples chain.
  expect_error(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(samples = list(matrix(0, 3L, 2L),
      .mock_sampling_chain(3L))), requested_chains = 2L), "do not match")
})

test_that("the sampling contract accepts real stored draws", {
  ## A one-draw Fixed_param result has zero warmup and is valid.
  expect_invisible(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(
      sim_chains = 1L, n_save = 1L, warmup2 = 0L,
      samples = list(.mock_sampling_chain(1L)),
      sim_extra = list(algorithm = "Fixed_param")
    ),
    requested_chains = 1L
  ))

  ## Ordinary multi-chain draws with warmup are valid.
  expect_invisible(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(
      n_save = c(4L, 4L), warmup2 = c(2L, 2L),
      samples = list(.mock_sampling_chain(4L), .mock_sampling_chain(4L))
    ),
    requested_chains = 2L
  ))

  ## Thinning and save_warmup = FALSE are represented by stored counts.
  expect_invisible(.assert_stan_sampling_contract(
    .mock_sampling_stanfit(
      n_save = c(3L, 3L), warmup2 = c(0L, 0L),
      samples = list(.mock_sampling_chain(3L), .mock_sampling_chain(3L)),
      sim_extra = list(iter = c(15L, 15L), thin = 5L,
                       save_warmup = FALSE)
    ),
    requested_chains = 2L
  ))
})
