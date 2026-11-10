# Extended opt-in integration. Source defines functions only; no HMC on sourcing.
rmi_hash <- function(x, file = FALSE) digest::digest(x, file = file, algo = "sha256")
rmi_internal <- function(name) getFromNamespace(name, "antitrustBayes")
rmi_config <- function() source("doc/validation/residual_model_integration_config.R")$value
rmi_remaining <- function(config = rmi_config()) as.numeric(difftime(
  as.POSIXct(substr(config$deadline_utc, 1L, 19L), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
  Sys.time(), units = "secs"))
rmi_softmax <- function(x) { y <- exp(x - max(x)); y / sum(y) }
rmi_sources <- function() {
  paths <- c("DESCRIPTION", sort(list.files("R", "[.]R$", full.names = TRUE)),
    sort(list.files("inst/stan", "[.]stan$", full.names = TRUE)),
    paste0("doc/validation/", c("residual_model_integration.R",
      "residual_model_integration_config.R", "residual_model_integration_native.R",
      "residual_model_integration_report.R", "residual_model_integration_cli.R",
      "residual_model_integration_collect.R",
      "relative_share_target_checks.R", "coordinated_release_gate.R")))
  setNames(vapply(paths, rmi_hash, character(1), file = TRUE), paths)
}

rmi_fixture <- function(case, config = rmi_config()) {
  ces <- case == "output_ces_panel"; M <- if (ces) 2L else 1L
  T <- if (case == "input_logit_snapshot") 1L else if (ces) 3L else 4L
  owner <- c("a", "a", "b", "b", "c", "d"); J <- length(owner)
  set.seed(config$seed + match(case, config$cases))
  q <- matrix(rnorm(M * 4L, 0, config$quality_sd), M, 4L)
  q <- q - rowMeans(q)
  innovations <- array(rnorm(M * T * J * 3L), c(M, T, J, 3L))
  rows <- truth <- validation <- list(); models <- list()
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    id <- paste0("m", m, "_", 2000L + tt)
    price <- 1.5 + .1 * seq_len(J) + .08 * sin(seq_len(J) * tt + m)
    quality <- q[m, match(owner, unique(owner))]
    shock <- config$demand_shock_sd * innovations[m, tt, , 1L]
    measurement <- config$measurement_sd * innovations[m, tt, , 2L]
    eta <- config$margin_sd * (config$rho_demand * innovations[m, tt, , 1L] +
      sqrt(1 - config$rho_demand^2) * innovations[m, tt, , 3L])
    structural <- if (ces) 6 else 3
    predictor <- quality + shock + if (ces) -structural * log(price) else structural * price
    shares <- .8 * rmi_softmax(predictor)
    native <- suppressWarnings(antitrust::synthetic_market(demand = if (ces) "ces" else "logit",
      supply = "bertrand", mode = "primitives", n_firms = 4L,
      n_products = c(2L, 2L, 1L, 1L), shares = c(shares, .2),
      prices = c(price, 1), reference_price = 1,
      parameters = if (ces) list(gamma = structural) else list(alpha = structural),
      output = ces, seed = config$seed + m * 100L + tt))
    model <- native@model; models[[id]] <- native
    actual <- as.numeric(antitrust::calcShares(model, preMerger = TRUE))[seq_len(J)]
    solved <- as.numeric(antitrust::calcPrices(model, preMerger = TRUE))[seq_len(J)]
    observed <- rmi_softmax(log(shares) + measurement)
    # Match the fitted conditional margin likelihood at its observed firm shares.
    observed_full <- .8 * observed
    firm_share <- as.numeric(outer(owner, owner, "==") %*% observed_full)
    if (ces) {
      revenue <- price * observed_full
      outside_revenue_share <- .2 / (.2 + sum(revenue))
      firm_revenue <- as.numeric(outer(owner, owner, "==") %*% revenue) /
        (.2 + sum(revenue))
      margin_mean_inv <- structural - (structural - 1) * firm_revenue
    } else margin_mean_inv <- price * structural * (1 - firm_share)
    margin <- exp(-eta) / margin_mean_inv
    rows[[id]] <- data.frame(prices = price, quantities = 1000 * observed,
      margins = margin, market = paste0("m", m), year = as.character(2000L + tt), group = owner)
    truth[[id]] <- data.frame(market = paste0("m", m), year = as.character(2000L + tt),
      product = seq_len(J), group = owner, price = price, share = shares,
      conditional_share = shares / sum(shares), quality = quality, demand_shock = shock,
      measurement = measurement, margin_residual = eta, margin_mean_inv = margin_mean_inv,
      structural = structural, s0 = .2, mc = as.numeric(model@mcPre)[seq_len(J)])
    native_foc <- if (ces) {
      rev <- c(price, 1) * c(shares, .2); rev <- rev / sum(rev)
      own <- outer(c(owner, "outside"), c(owner, "outside"), "==")
      coefficient <- structural - (structural - 1) * as.numeric(own %*% rev)
      max(abs((c(price, 1) - model@mcPre) / c(price, 1) * coefficient - 1))
    } else native@diagnostics$synthetic$foc_residual
    validation[[id]] <- data.frame(case, unit = id,
      generator_foc = native_foc,
      share_error = max(abs(actual - shares)), price_error = max(abs(solved - price)),
      repeated_owner = anyDuplicated(owner) > 0L)
  }
  validation <- do.call(rbind, validation)
  stopifnot(all(validation$generator_foc < 1e-8), all(validation$share_error < 1e-8),
    all(validation$price_error < 1e-6), all(validation$repeated_owner))
  list(case = case, data = do.call(rbind, rows), truth = do.call(rbind, truth),
    validation = validation, generating_models = models, quality = q,
    innovations = innovations, demand_model = if (ces) "ces" else "logit",
    market_side = if (ces) "output" else "input", M = M, T = T,
    interpretation = paste("Model-matched working share/margin likelihood with known quality and",
      "demand shocks. Margin means use observed firm shares, as production does.",
      "Variance component recovery is prior dependent, not component identification."))
}

rmi_arguments <- function(fixture, model, directory, config = rmi_config()) {
  c(as.list(fixture$data), list(demand_model = fixture$demand_model,
    market_side = fixture$market_side, supply_model = "bertrand",
    share_basis = "relative", share_residual_model = model,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market",
    use_rho = TRUE, use_alpha_market_pool = fixture$M > 1L,
    use_s0_market_pool = fixture$T > 1L, use_s0_market_year_pool = FALSE,
    priors = list(prior_sigma_measurement = config$measurement_prior_median,
      prior_log_sigma_measurement_sd = config$measurement_prior_log_sd),
    stan_options = list(chains = config$chains, iter = config$warmup + config$retained,
      warmup = config$warmup, cores = config$cores, seed = config$seed,
      sig_figs = 18L, init = 0, refresh = 100L, output_dir = directory,
      control = list(adapt_delta = config$adapt_delta, max_treedepth = config$max_treedepth))))
}

rmi_prepare <- function(output) {
  config <- rmi_config(); stopifnot(rmi_remaining(config) > 0)
  if (file.exists(file.path(output, "fixtures.rds"))) stop("Frozen fixtures already exist.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  fixtures <- setNames(lapply(config$cases, rmi_fixture, config = config), config$cases)
  saveRDS(fixtures, file.path(output, "fixtures.rds"))
  write.csv(do.call(rbind, lapply(fixtures, `[[`, "validation")),
    file.path(output, "fixture_validation.csv"), row.names = FALSE)
  calls <- list()
  for (case in config$cases) for (model in config$models) {
    id <- paste(case, model, sep = "__")
    calls[[id]] <- rmi_arguments(fixtures[[case]], model, file.path(output, id), config)
  }
  m <- list(config = config, fixtures_sha256 = rmi_hash(file.path(output, "fixtures.rds"), TRUE),
    fixture_data_hashes = vapply(fixtures, function(x) rmi_hash(x$data), character(1)),
    fixture_truth_hashes = vapply(fixtures, function(x) rmi_hash(x$truth), character(1)),
    calls = calls, calls_sha256 = rmi_hash(calls), attempts = list(), status = "fixtures_frozen",
    started_utc = config$start_utc, deadline_utc = config$deadline_utc,
    no_bridge_no_bma = TRUE, no_extra_fits = TRUE)
  saveRDS(m, file.path(output, "manifest.rds")); invisible(m)
}

rmi_seal <- function(output) {
  m <- readRDS(file.path(output, "manifest.rds"))
  if (!is.null(m$source_hashes)) stop("Scientific target is already sealed.")
  stopifnot(identical(m$fixtures_sha256, rmi_hash(file.path(output, "fixtures.rds"), TRUE)))
  m$source_hashes <- rmi_sources(); m$priors <- bayes_prior_manifest()
  m$priors_sha256 <- rmi_hash(m$priors); m$sealed_at <- Sys.time()
  m$resolved_call_priors <- lapply(m$calls, function(args) utils::modifyList(m$priors, args$priors))
  m$resolved_call_prior_hashes <- vapply(m$resolved_call_priors, rmi_hash, character(1))
  m$prepared_stan_data <- lapply(m$calls, function(args) {
    helper <- rmi_internal(if (args$demand_model == "ces") "prepare_ces_stan_data" else "prepare_stan_data")
    selected <- args[intersect(names(args), names(formals(helper)))]
    selected$supply_model <- 1L; selected$use_cutoff <- FALSE; selected$use_hmt <- FALSE
    selected$priors <- utils::modifyList(m$priors, args$priors)
    do.call(helper, selected)
  })
  m$prepared_stan_data_hashes <- vapply(m$prepared_stan_data, rmi_hash, character(1))
  m$packages <- setNames(lapply(c("antitrust", "coordination", "cmdstanr", "posterior"),
    function(p) list(version = as.character(packageVersion(p)), path = find.package(p))),
    c("antitrust", "coordination", "cmdstanr", "posterior"))
  m$cmdstan <- list(path = cmdstanr::cmdstan_path(), version = as.character(cmdstanr::cmdstan_version()))
  m$status <- "scientific_target_sealed"; saveRDS(m, file.path(output, "manifest.rds")); invisible(m)
}
rmi_assert <- function(m, output) stopifnot(!is.null(m$source_hashes),
  identical(m$source_hashes, rmi_sources()), identical(m$calls_sha256, rmi_hash(m$calls)),
  identical(m$priors_sha256, rmi_hash(bayes_prior_manifest())),
  identical(m$fixtures_sha256, rmi_hash(file.path(output, "fixtures.rds"), TRUE)))

rmi_diagnose <- function(fit) {
  source("doc/validation/coordinated_release_gate.R", local = TRUE)
  active <- release_gate_active(fit)
  tab <- as.data.frame(posterior::summarise_draws(fit@stanfit$draws(active)))
  native <- fit@stanfit$sampler_diagnostics()
  data.frame(Max_Rhat = max(tab$rhat), Min_Bulk_ESS = min(tab$ess_bulk),
    Divergences = sum(native[, , "divergent__"]),
    Treedepth_Hits = sum(native[, , "treedepth__"] >= rmi_config()$max_treedepth),
    Active_parameters = nrow(tab), Pass = all(is.finite(tab$rhat)) && max(tab$rhat) <= 1.01 &&
      min(tab$ess_bulk) >= 400 && sum(native[, , "divergent__"]) == 0 &&
      sum(native[, , "treedepth__"] >= rmi_config()$max_treedepth) == 0)
}

rmi_worker <- function(output, id) {
  m <- readRDS(file.path(output, "manifest.rds")); rmi_assert(m, output)
  stopifnot(id %in% names(m$calls), rmi_remaining(m$config) > 15)
  dest <- file.path(output, id); dir.create(dest, showWarnings = FALSE)
  if (file.exists(file.path(dest, "fit.rds"))) stop("Completed fit preserved.")
  started <- Sys.time()
  fit <- do.call(bayes_est, m$calls[[id]])
  saveRDS(fit, file.path(dest, "fit.rds"))
  write.csv(rmi_diagnose(fit), file.path(dest, "diagnostics.csv"), row.names = FALSE)
  saveRDS(list(started = started, ended = Sys.time(), fit_sha256 = rmi_hash(file.path(dest, "fit.rds"), TRUE),
    data_sha256 = rmi_hash(fit@stan_data), priors_sha256 = m$priors_sha256,
    source_sha256 = fit@metadata$stan_source_hash), file.path(dest, "seal.rds"))
  source("doc/validation/residual_model_integration_report.R", local = TRUE)
  rmi_report_fit(output, id)
  invisible(fit)
}

rmi_run <- function(output) {
  m <- readRDS(file.path(output, "manifest.rds")); rmi_assert(m, output)
  if (!identical(Sys.getenv("AB_RESIDUAL_HMC_AUTHORIZED"), "root_ack")) stop("Root green acknowledgment required.")
  save <- function() saveRDS(m, file.path(output, "manifest.rds"))
  for (id in names(m$calls)) {
    if (id %in% names(m$attempts)) next
    if (rmi_remaining(m$config) < 60) {m$status <- "absolute_budget_exhausted"; save(); break}
    stopifnot(length(m$attempts) < m$config$maximum_fits)
    m$attempts[[id]] <- list(status = "started", started = Sys.time()); save()
    args <- c("--signal=TERM", "--kill-after=5s", paste0(floor(rmi_remaining(m$config)) - 5L, "s"),
      file.path(R.home("bin"), "Rscript"), "--vanilla", "doc/validation/residual_model_integration_cli.R",
      "worker", output, id)
    result <- processx::run("timeout", args, env = "current", error_on_status = FALSE,
      stdout = file.path(output, paste0(id, ".log")), stderr_to_stdout = TRUE, cleanup_tree = TRUE)
    m$attempts[[id]]$status <- if (result$status == 0) "completed" else if (result$status == 124) "timed_out" else "failed"
    m$attempts[[id]]$exit_status <- result$status; m$attempts[[id]]$ended <- Sys.time(); save()
    cat(id, m$attempts[[id]]$status, format(Sys.time()), "\n"); flush.console()
  }
  invisible(m)
}
