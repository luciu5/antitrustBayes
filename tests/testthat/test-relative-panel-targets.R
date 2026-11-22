# Compare native outputs and derivatives at the same fixed parameters. The
# one-row fixed_param run initializes a point and emits all generated quantities.
test_that("Logit and CES panel kernels preserve densities, gradients and generated quantities", {
  directory <- Sys.getenv("RELATIVE_PANEL_EXE_DIR")
  if (!nzchar(directory))
    skip_if(Sys.getenv("RUN_RELATIVE_PANEL_STAN") != "true", "Opt-in compiled panel comparisons")
  executable <- function(kernel) if (nzchar(directory)) file.path(directory, kernel) else
    .cmdstan_model_for_kernel(kernel)$exe_file()
  work <- tempfile("panel-targets-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  data_path <- file.path(work, "data.json")
  run <- function(kernel, args, output) {
    log <- system2(executable(kernel), c(args, "data", paste0("file=", shQuote(data_path)),
      "output", paste0("file=", shQuote(output)), "sig_figs=18"), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(log, "status"))) stop(paste(tail(log, 10), collapse = "\n"))
    as.matrix(utils::read.csv(output, comment.char = "#", check.names = FALSE))
  }
  compare <- function(x, y, case) {
    expect_identical(dimnames(x), dimnames(y), info = case)
    expect_identical(is.finite(x), is.finite(y), info = case)
    expect_identical(x[!is.finite(x)], y[!is.finite(y)], info = case)
    finite <- is.finite(x)
    expect_lt(max(abs(x[finite] - y[finite])), 1e-7, label = case)
  }
  designs <- c("one", "unequal", "missing_period", "partial", "ragged", "global",
               "global_ragged", "global_one", "global_scales", "uncorrelated", "demand_shock", "sampled")
  for (demand in c("logit", "ces")) for (structure in c("simultaneous", "fringe", "stackelberg")) {
    rich <- structure != "simultaneous"
    base <- paste(demand, if (rich) "stackelberg" else "supply", sep = "_")
    for (design in designs) {
      x <- expand.grid(group = 1:5, year = 1:3, market = 1:2)
      if (design == "global_one") x <- subset(x, market == 1)
      if (design == "one") x <- subset(x, market == 1 & year == 1)
      if (design %in% c("unequal", "missing_period", "partial", "ragged", "sampled"))
        x <- subset(x, market == 1 | (year < 3 & group < 5))
      if (design %in% c("ragged", "global_ragged")) x <- x[-1, ]
      x <- x[order((seq_len(nrow(x)) * 17) %% nrow(x)), ]
      observed <- rep(TRUE, nrow(x))
      if (design == "partial") observed[1] <- FALSE
      if (design %in% c("missing_period", "global_scales")) observed[x$market == 1 & x$year == 2] <- FALSE
      global <- startsWith(design, "global")
      args <- list(prices = 1 + .08 * x$group + .02 * x$year,
        quantities = 4 + (x$group * x$market) %% 9 + .03 * x$year,
        margins = ifelse(observed, .15 + .015 * x$group + .002 * x$year, NA_real_),
        market = x$market, year = x$year, group = x$group, margin_observed = observed,
        supply_model = if (structure == "fringe") 1L else 3L,
        share_basis = "relative", use_s0_market_pool = TRUE, structure = structure,
        quality_scope = if (global) "global" else "market",
        share_residual_model = if (design == "demand_shock") "demand_shock" else "combined",
        use_market_resid_scale = design == "global_scales",
        use_rho = !(design == "uncorrelated" || (rich && global && design != "global_one")))
      if (design == "sampled") args$quality_parameterization <- "sampled"
      if (rich) args <- c(args, list(partition_prior = "ordered_k", k_support = 1:3,
                                     k_prior = c(.5, 0, .5)))
      d <- do.call(if (demand == "logit") prepare_stan_data else prepare_ces_stan_data, args)
      fast <- if (demand == "logit" && rich && .uses_balanced_panel_kernel(d))
        "logit_ordered_panel" else paste0(base, "_relative_panel")
      case <- paste(demand, structure, design)
      cmdstanr::write_stan_json(d[.nutpie_stan_data_fields(base)], data_path)
      draw <- file.path(work, "draw.csv")
      alt <- file.path(work, "alt.csv")
      fixed <- c("sample", "num_samples=1", "num_warmup=0", "algorithm=fixed_param",
                 "random", "seed=375", "init=0.4")
      compare(run(base, fixed, draw), run(fast, fixed, alt), case)
      grad <- c("log_prob", paste0("constrained_params=", shQuote(draw)), "jacobian=1")
      compare(run(base, grad, alt), run(fast, grad, alt), case)
    }
  }
})
