panel_test_args <- function(conduct = 1L, structure = "fringe", mask = FALSE,
                            remove = FALSE, duplicate = FALSE,
                            markets = 4L, periods = 4L) {
  x <- expand.grid(firm = 1:7, year = seq_len(periods), market = seq_len(markets))
  # Shuffle rows and give markets different firm rankings.
  x <- x[order((seq_len(nrow(x)) * 53L) %% nrow(x)), ]
  if (remove) x <- x[-1, ]
  if (duplicate) x <- rbind(x, x[1, ])
  margins <- .15 + .015 * x$firm + .002 * x$year
  if (mask) margins[1] <- NA_real_
  list(prices = 1 + .08 * x$firm + .02 * x$year,
    quantities = 4 + (x$firm * x$market) %% 9 + .03 * x$year, margins = margins,
    market = x$market, year = x$year, group = x$firm, market_side = "output",
    margin_observed = is.finite(margins),
    supply_model = conduct, structure = structure, partition_prior = "ordered_k",
    k_support = 1:6, k_prior = (1:6) / 21, share_basis = "relative",
    use_s0_market_pool = TRUE, share_residual_model = "combined", use_rho = TRUE)
}

panel_test_data <- function(...) do.call(prepare_stan_data, panel_test_args(...))

test_that("bayes_est selects the compiled kernel from the actual panel", {
  local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) stop("selected ", kernel, "!"))
  for (mask in c(FALSE, TRUE)) {
    args <- panel_test_args(mask = mask)
    args$supply_model <- "bertrand"
    args$demand_model <- "logit"
    expect_error(do.call(bayes_est, args), if (mask) "selected logit_stackelberg_relative_panel!" else
      "selected logit_ordered_panel!")
  }
})

test_that("only an explicit zero prior selects the sparse ordered-k likelihood", {
  local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) stop("selected ", kernel, "!"))
  args <- panel_test_args()
  args$supply_model <- "bertrand"
  args$k_prior <- NULL
  expect_error(do.call(bayes_est, args), "selected logit_ordered_panel!")
  args$k_support <- 3L
  expect_error(do.call(bayes_est, args), "selected logit_ordered_panel!")
  args$k_support <- 1:6
  args$k_prior <- c(rep(0.2, 5), 1e-100)
  expect_error(do.call(bayes_est, args), "selected logit_ordered_panel!")
  args$k_prior <- as.numeric(args$k_support == 3L)
  expect_error(do.call(bayes_est, args), "selected logit_ordered_panel_sparse!")
  args$demand_model <- "ces"
  expect_error(do.call(bayes_est, args), "selected ces_stackelberg_relative_panel_sparse!")
  args$share_basis <- "absolute"
  expect_error(do.call(bayes_est, args), "selected ces_stackelberg_sparse!")
})

test_that("only complete paired ordered panels enter the balanced kernel", {
  d <- panel_test_data()
  expect_true(.uses_balanced_panel_kernel(d))
  for (bad in list(panel_test_data(mask = TRUE), panel_test_data(remove = TRUE),
                   panel_test_data(duplicate = TRUE))) expect_false(.uses_balanced_panel_kernel(bad))
  for (key in c("share_basis", "likelihood_mode", "quality_parameterization",
                "quality_scope")) {
    bad <- d
    bad[[key]] <- if (key == "use_rho") 0L else if (d[[key]] == 1L) 2L else 1L
    expect_false(.uses_balanced_panel_kernel(bad))
  }
  expect_true(.uses_balanced_panel_kernel(panel_test_data(markets = 2L, periods = 3L)))
  bad <- d
  bad$obs_to_quality[1L] <- bad$obs_to_quality[1L] + 1L
  expect_false(.uses_balanced_panel_kernel(bad))
})

test_that("balanced kernel reuses the original generated quantities and has its own evidence seal", {
  panel <- .cmdstan_kernel_file("logit_stackelberg_panel")
  old <- .cmdstan_kernel_file("logit_stackelberg")
  quantities <- function(path) {
    lines <- readLines(path)
    lines[which(lines == "generated quantities {"):length(lines)]
  }
  expect_identical(quantities(panel), quantities(old))
  meta <- list(backend = "cmdstanr", inference = "hmc", model_spec_version = 15L,
    stan_kernel = "logit_stackelberg_panel", stan_source_hash = digest::digest(file = panel, algo = "sha256"),
    stan_source_identity = .bayes_stan_source_identity(panel))
  expect_identical(.bayes_cmdstan_source_contract(meta, readLines(panel))$current_hash, meta$stan_source_hash)
  expect_error(.bayes_cmdstan_source_contract(meta, readLines(old)), "does not match")
})

test_that("ordinary and ordered Logit and CES use distinct sealed relative kernels", {
  local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) stop("selected ", kernel, "!"))
  for (demand in c("logit", "ces")) for (structure in c("simultaneous", "fringe")) {
    args <- panel_test_args(markets = 2L, periods = 3L)
    args$supply_model <- "bertrand"
    args$demand_model <- demand
    args$structure <- structure
    if (structure == "simultaneous")
      args[c("k_support", "k_prior", "partition_prior")] <- NULL
    base <- paste(demand, if (structure == "fringe") "stackelberg" else "supply", sep = "_")
    kernel <- if (base == "logit_stackelberg") "logit_ordered_panel" else
      paste0(base, "_relative_panel")
    expect_error(do.call(bayes_est, args), paste0("selected ", kernel, "!"))
    source <- .cmdstan_kernel_file(kernel)
    meta <- list(backend = "cmdstanr", inference = "hmc", model_spec_version = 15L,
      stan_kernel = kernel, stan_source_hash = digest::digest(file = source, algo = "sha256"),
      stan_source_identity = .bayes_stan_source_identity(source))
    expect_identical(.bayes_cmdstan_source_contract(meta, readLines(source))$kernel, kernel)
    expect_error(.bayes_cmdstan_source_contract(meta, readLines(.cmdstan_kernel_file(base))),
                 "does not match")
  }
})

test_that("validation manifests verify composed kernels from the requested checkout", {
  repository <- test_path("..", "..")
  collector <- file.path(repository, "doc", "validation", "relative_share_validation_collect.R")
  skip_if_not(file.exists(collector), "Research collector requires a source checkout")
  scripts <- new.env(parent = environment())
  sys.source(collector, scripts)
  sys.source(file.path(repository, "doc", "validation", "relative_share_backend_validation.R"), scripts)
  manifest <- withr::with_dir(repository, scripts$relative_validation_kernel_hashes())
  expect_identical(scripts$.relative_collect_source_hashes(names(manifest), repository), unname(manifest))
  runtime_paths <- c(file.path("R", list.files(file.path(repository, "R"), "[.]R$")),
    "doc/validation/nutpie_backend_validation.R")
  runtime_hashes <- vapply(file.path(repository, runtime_paths), digest::digest, character(1),
    file = TRUE, algo = "sha256")
  sampling <- list(runtime = list(source_hashes = setNames(runtime_hashes, runtime_paths)),
    kernel_hashes = manifest)
  expect_true(all(scripts$.relative_collect_source_table(sampling, repository)$matches))
  missing <- "kernel:ces_supply_relative_panel"
  sampling$kernel_hashes <- manifest[names(manifest) != missing]
  audit <- scripts$.relative_collect_source_table(sampling, repository)
  expect_identical(audit$path[!audit$matches], missing)

  checkout <- tempfile("panel-checkout-")
  on.exit(unlink(checkout, recursive = TRUE))
  stan <- file.path(checkout, "inst", "stan")
  dir.create(stan, recursive = TRUE)
  original <- file.path(repository, "inst", "stan")
  expect_true(all(file.copy(file.path(original, c("ces_supply.stan", "include")), stan, recursive = TRUE)))
  keys <- c("inst/stan/ces_supply.stan", "kernel:ces_supply_relative_panel")
  fragment <- file.path(stan, "include", "relative_panel_functions.stan")
  cat("\n// Changed include\n", file = fragment, append = TRUE)
  current <- scripts$.relative_collect_source_hashes(keys, checkout)
  expect_identical(current[1], unname(manifest[keys[1]]))
  expect_false(is.na(current[2]))
  expect_false(identical(current[2], unname(manifest[keys[2]])))
  expect_identical(scripts$.relative_collect_source_hashes("kernel:unknown", checkout), NA_character_)
})

test_that("compiled M4T4 density, gradient and all ordered-k scores match the generic kernel", {
  # Optional when CmdStan executables have been prepared. Hundreds of fixed
  # parameter evaluations use log_prob only, including its Jacobian.
  executables <- Sys.getenv(c("PANEL_GENERIC_EXE", "PANEL_FAST_EXE"))
  if (any(!nzchar(executables))) {
    skip_if(Sys.getenv("RUN_BALANCED_PANEL_STAN") != "true", "Opt-in compiled fixed-point comparison")
    executables <- vapply(c("logit_stackelberg", "logit_ordered_panel"),
      function(k) .cmdstan_model_for_kernel(k)$exe_file(), character(1))
  }
  expect_true(all(file.exists(executables)))
  work <- tempfile("panel-density-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  data_path <- file.path(work, "data.json")
  points_path <- file.path(work, "points.json")
  base <- c(-.4, .1, -1, .1, -.2, .15, .05, -.7, -.8, .1, -.1, .3, -.4, .45, -2, -1.5)
  value <- function(exe) {
    output <- file.path(work, "output.csv")
    log <- system2(exe, c("log_prob", paste0("unconstrained_params=", shQuote(points_path)),
      "jacobian=1", "data", paste0("file=", shQuote(data_path)), "output",
      paste0("file=", shQuote(output)), "sig_figs=18"), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(log, "status")) && attr(log, "status") != 0L)
      stop(paste(log, collapse = "\n"))
    as.matrix(utils::read.csv(output, comment.char = "#", check.names = FALSE))
  }
  for (supply in c(1L, 3L)) for (structure in c("fringe", "stackelberg")) {
    d <- panel_test_data(supply, structure)
    expect_true(.uses_balanced_panel_kernel(d))
    fields <- .nutpie_stan_data_fields("logit_stackelberg")
    for (k in 0:6) {
      # Isolate each state's normalized score in market 1 while retaining
      # the full, nonuniform mixture in the other three markets.
      data <- d[fields]
      if (k > 0) {
        idx <- d$partition_unit_start[1] + 0:5
        data$partition_prior_probability[idx] <- as.numeric(seq_len(6) == k)
        data$partition_prior_logprob[idx] <- 0 # zero-probability states are guarded in Stan
      }
      cmdstanr::write_stan_json(data, data_path)
      n <- if (k == 0L) 600L else 4L
      points <- t(vapply(seq_len(n), function(i) base + .45 * sin(seq_along(base) * i),
        numeric(length(base))))
      points[, 14L] <- 2 * sin(seq_len(n)) # both signs, including strongly correlated residuals
      cmdstanr::write_stan_json(list(params_r = points), points_path)
      elapsed <- numeric(2L)
      results <- vector("list", 2L)
      for (j in seq_along(executables))
        elapsed[j] <- system.time(results[[j]] <- value(executables[j]))[["elapsed"]]
      delta <- abs(results[[1]] - results[[2]])
      expect_true(all(is.finite(delta)))
      expect_lt(max(delta), 1e-7)
      if (k == 0L) cat(sprintf("\n%s/%s | %d gradients | generic %.3fs panel %.3fs | max abs LP %.3g gradient %.3g\n",
        supply, structure, n, elapsed[1], elapsed[2], max(delta[, 1]), max(delta[, -1])))
    }
  }
})
