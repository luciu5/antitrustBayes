# Opt-in fixed-point native target checks. Sourcing defines functions only;
# no compilation, installation, HMC, or production writes occur.
.relative_target_internal <- function(name) getFromNamespace(name, "antitrustBayes")
.relative_target_hash <- function(path) digest::digest(file = path, algo = "sha256")
.relative_target_bracket <- function(x) vapply(strsplit(x, ".", fixed = TRUE), function(z)
  if (length(z) == 1L) z else paste0(z[1L], "[", paste(z[-1L], collapse = ","), "]"), character(1))
.relative_target_field <- function(values, name) {
  hit <- names(values) == name | startsWith(names(values), paste0(name, "["))
  as.numeric(values[hit])
}
.relative_target_lse <- function(x) max(x) + log(sum(exp(x - max(x))))
.relative_target_normal <- function(residual, covariance) {
  if (!length(residual)) return(0)
  ch <- chol(covariance)
  white <- forwardsolve(t(ch), residual)
  -.5 * (length(residual) * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(white^2))
}
.relative_target_dense <- function(d, values, margin_mean, rows = seq_len(d$N)) {
  field <- function(name) .relative_target_field(values, name)
  cell <- d$market_year_idx[rows]
  ee <- outer(cell, unique(cell), "==") * 1
  q <- if (nrow(ee) == ncol(ee)) matrix(0, nrow(ee), 0L) else
    qr.Q(qr(ee), complete = TRUE)[, seq.int(ncol(ee) + 1L, nrow(ee)), drop = FALSE]
  ss <- field("sigma_share_market")[d$market[rows]]
  sm <- field("sigma_margin_market")[d$market[rows]]
  rho <- if (d$use_rho == 1L && d$likelihood_mode == 1L && d$N_margin > 0L) field("rho") else 0
  share_covariance <- diag(ss^2, length(rows))
  if (d$quality_parameterization == 2L) {
    p <- diag(d$N_quality)
    for (b in seq_len(d$N_quality_block)) {
      ii <- seq.int(d$quality_block_start[b], length.out = d$quality_block_n[b])
      p[ii, ii] <- p[ii, ii] - 1 / length(ii)
    }
    z <- diag(d$N_quality)[d$obs_to_quality[rows], , drop = FALSE]
    share_covariance <- share_covariance + field("sigma_b_strat_economic")^2 * z %*% p %*% t(z)
  }
  share_residual <- as.numeric(crossprod(q, d$log_shareIn[rows] - field("relative_share_predictor")[rows]))
  projected <- crossprod(q, share_covariance %*% q)
  share <- if (d$likelihood_mode != 3L) .relative_target_normal(share_residual, projected) else 0
  observed <- if (d$likelihood_mode != 2L) which(d$margin_obs[rows] == 1L) else integer()
  if (!length(observed)) return(c(share = share, joint = share))
  mr <- d$log_marginNormInv[rows[observed]] - margin_mean[rows[observed]]
  vm <- diag(sm[observed]^2, length(observed))
  if (d$likelihood_mode == 3L) return(c(share = 0, joint = .relative_target_normal(mr, vm)))
  cross <- matrix(0, length(rows), length(observed))
  cross[cbind(observed, seq_along(observed))] <- rho * ss[observed] * sm[observed]
  cross <- crossprod(q, cross)
  joint <- rbind(cbind(projected, cross), cbind(t(cross), vm))
  c(share = share, joint = .relative_target_normal(c(share_residual, mr), joint))
}

.relative_target_fixture <- function(kernel, markets = 1L, periods = 1L,
    quality = "collapsed", rho = TRUE, scope = "market", conduct = 1L,
    structure = "simultaneous", mixture = FALSE, ragged = FALSE,
    market_side = "output", mode = 1L, smooth = FALSE) {
  ces <- startsWith(kernel, "ces")
  a <- expand.grid(product = 1:4, period = seq_len(periods), market = seq_len(markets))
  if (ragged) a <- a[!(a$product == 4L & a$period %% 2L == 0L), , drop = FALSE]
  a <- a[rev(seq_len(nrow(a))), , drop = FALSE]
  prices <- exp(.1 * a$product + .02 * a$period + .04 * a$market)
  args <- list(prices = prices, quantities = exp(1 - .3 * a$product + .02 * a$period),
    margins = .18 + .01 * a$product, market = as.character(a$market), year = as.character(a$period),
    group = paste0("p", a$product), supply_model = as.integer(conduct), market_side = market_side,
    use_cutoff = smooth, use_rho = rho, use_hmt = FALSE, use_alpha_market_pool = markets > 1L,
    use_s0_market_pool = periods > 1L, use_s0_market_year_pool = FALSE,
    margin_observed = seq_len(nrow(a)) %% 3L != 0L,
    outside_price = .75 + .04 * seq_len(markets * periods), likelihood_mode = mode,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = quality, quality_scale_parameterization = "legacy",
    quality_scope = scope, share_basis = "relative", priors = list())
  if (grepl("stackelberg$", kernel)) args <- c(args, list(structure = structure,
    partition_prior = "ordered_k", k_support = if (mixture) c(2L, 1L) else 2L,
    k_prior = if (mixture) c(.3, .7) else 1))
  d <- do.call(.relative_target_internal(if (ces) "prepare_ces_stan_data" else "prepare_stan_data"), args)
  list(data = d, arguments = args, prices = prices, id = paste(kernel, markets, periods,
    quality, rho, scope, conduct, structure, mixture, ragged, market_side, mode, smooth, sep = "_"))
}

.relative_target_u <- function(model) {
  names <- as.character(unlist(model$param_unc_names()))
  u <- setNames(.12 * sin(seq_along(names)), names)
  u[grepl("^sigma_b_strat$", names)] <- log(.8)
  u[grepl("^log_sigma_share$", names)] <- log(.23)
  u[grepl("^log_sigma_margin$", names)] <- log(.31)
  u[grepl("^rho_raw", names)] <- atanh(.63 / .95)
  u[grepl("^s0_logit|^mu_s0_pool_vec", names)] <- qlogis(.25)
  u[grepl("^sigma_.*_vec", names)] <- log(.2)
  u[grepl("^mu_log_a", names)] <- .1
  u[grepl("^mu_b_fringe_vec", names)] <- .2
  u
}
.relative_target_values <- function(model, u) {
  value <- model$param_constrain(u, include_tp = TRUE, include_gq = TRUE, rng = model$new_rng(137L))
  names(value) <- .relative_target_bracket(unlist(model$param_names(include_tp = TRUE, include_gq = TRUE)))
  value
}
.relative_target_margin <- function(d, values, ces) {
  f <- function(name) .relative_target_field(values, name)
  m <- d$market; t <- d$market_year_idx; s0 <- f("s0_prob")[t]
  if (ces) {
    gamma <- f("gamma_structural_market")[m]; r0 <- f("R0_prob")[t]
    r <- d$firm_revShare * (1 - r0)
    coefficient <- switch(as.character(d$supply_model),
      "1" = gamma - (gamma - 1) * r, "2" = gamma,
      "3" = 1 / (1 / gamma + (gamma - 1) * r / (gamma * (1 + (gamma - 1) * r0))),
      "4" = 1 / (1 - (1 - r)^(1 / (gamma - 1))),
      "5" = (gamma - (gamma - 1) * r) / (1 - f("bargaining_power")))
    weight <- if (d$use_cutoff == 1L) plogis(-(d$firm_revShare - f("cutoff_share")[m]) / d$cutoff_scale) else 0
    pred <- weight * gamma + (1 - weight) * coefficient
  } else {
    alpha <- f("alpha_structural_market")[m]; r <- d$firm_shareIn * (1 - s0)
    coefficient <- switch(as.character(d$supply_model), "1" = 1 - r, "2" = rep(1, d$N),
      "3" = s0 / (s0 + r), "4" = r / -log1p(-r),
      "5" = (f("bargaining_power") / (1 - f("bargaining_power")) * r -
        (1 - r) * log1p(-r)) / -log1p(-r))
    weight <- if (d$use_cutoff == 1L) plogis(-(d$firm_shareIn - f("cutoff_share")[m]) / d$cutoff_scale) else 0
    pred <- d$price * alpha * (weight + (1 - weight) * coefficient)
  }
  log(pred + 1e-9)
}
.relative_target_exact_margin <- function(d, values, ces, structure, k) {
  draws <- list(a_market = matrix(.relative_target_field(values, "a_market"), nrow = 1L),
    s0_prob = matrix(.relative_target_field(values, "s0_prob"), nrow = 1L))
  fit <- methods::new("BayesLogit", stan_data = d, data = data.frame(),
    supply_model = as.integer(d$supply_model), demand_model = if (ces) "ces" else "logit",
    market_side = "output", price_sd = if (ces) d$log_price_sd else d$price_sd,
    metadata = list(structure = structure), stanfit = NULL)
  log(as.numeric(.relative_target_internal(".stack_pure_foc_margin_prediction")(fit,
    parsed = list(y = d$log_shareIn, n_draw = 1L), draws = draws, draw_indices = 1L,
    realized_k = matrix(rep(k, d$N_market), nrow = 1L), realized_k_scope = "market")) + 1e-9)
}

relative_share_target_checks <- function(artifact_dir, output_dir = artifact_dir,
    kernels = c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg"),
    prefix = "relative-native-targets", stop_on_failure = TRUE) {
  stopifnot(requireNamespace("antitrustBayes", quietly = TRUE), requireNamespace("bridgestan", quietly = TRUE))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  results <- list(); warnings <- character(); sources <- character()
  add <- function(kernel, case, test, actual = NA_real_, expected = NA_real_, tolerance = 1e-8, error = NULL) {
    ok <- is.null(error) && is.finite(actual) && is.finite(expected) && abs(actual - expected) <= tolerance
    results[[length(results) + 1L]] <<- data.frame(kernel, case, test, actual, expected,
      absolute_error = abs(actual - expected), tolerance, status = if (ok) "passed" else "failed",
      error = if (is.null(error)) NA_character_ else error, stringsAsFactors = FALSE)
  }
  for (kernel in kernels) {
    artifact <- file.path(artifact_dir, paste0(kernel, "-nutpie-artifact.rds"))
    if (!file.exists(artifact)) {add(kernel, "artifact", "available", error = "Compiled artifact unavailable; no compilation attempted"); next}
    metadata <- readRDS(artifact)
    source <- .relative_target_internal(".cmdstan_kernel_file")(kernel)
    stopifnot(identical(.relative_target_hash(source), metadata$source_hash),
      identical(.relative_target_hash(metadata$lib_path), metadata$lib_hash))
    sources[kernel] <- metadata$source_hash
    alias <- file.path(output_dir, paste0("target_", kernel, "_", substr(metadata$lib_hash, 1L, 16L), ".so"))
    if (!file.exists(alias)) stopifnot(file.copy(metadata$lib_path, alias))
    stopifnot(identical(.relative_target_hash(alias), metadata$lib_hash))
    fixtures <- list(); append_case <- function(...) {x <- .relative_target_fixture(kernel, ...);fixtures[[x$id]] <<- x}
    exact <- grepl("stackelberg$", kernel); ces <- startsWith(kernel, "ces")
    for (M in c(1L, 4L)) for (T in c(1L, 4L)) for (quality in c("sampled", "collapsed"))
      for (rho in c(FALSE, TRUE)) for (scope in c("market", "global")) {
        if (exact && scope == "global" && M > 1L && quality == "collapsed" && rho) next
        if (exact) for (structure in c("fringe", "stackelberg")) for (mixture in c(FALSE, TRUE))
          append_case(markets = M, periods = T, quality = quality, rho = rho, scope = scope,
            structure = structure, mixture = mixture, conduct = if (mixture) 3L else 1L)
        else append_case(markets = M, periods = T, quality = quality, rho = rho,
          scope = scope, conduct = if (T == 1L) 1L else 3L, ragged = T > 1L)
      }
    if (!exact) {
      for (rho in c(FALSE, TRUE)) for (quality in c("sampled", "collapsed"))
        append_case(conduct = 2L, quality = quality, rho = rho)
      for (mode in c(2L, 3L)) append_case(mode = mode, quality = "sampled")
      append_case(smooth = TRUE, conduct = 1L)
      if (!ces) append_case(market_side = "input", conduct = 3L)
    }
    for (fixture in fixtures) {
      case <- fixture$id;d <- fixture$data
      tryCatch(withCallingHandlers({
        data_file <- file.path(output_dir, paste0(prefix, "-data.json"))
        cmdstanr::write_stan_json(d, data_file)
        native <- bridgestan::StanModel$new(alias, data_file, 137L)
        u <- .relative_target_u(native);values <- .relative_target_values(native, u)
        field <- function(name) .relative_target_field(values, name)
        constrained <- native$param_constrain(u)
        roundtrip <- native$param_unconstrain(constrained)
        add(kernel, case, "parameter_roundtrip", max(abs(as.numeric(u) - roundtrip)), 0, 1e-10)
        log_target <- native$log_density(u, propto = FALSE, jacobian = TRUE)
        add(kernel, case, "finite_normalized_target", as.numeric(is.finite(log_target)), 1, 0)
        gradient <- native$log_density_gradient(u, propto = FALSE, jacobian = TRUE)$gradient
        add(kernel, case, "finite_target_gradient", as.numeric(length(gradient) == length(u) && all(is.finite(gradient))), 1, 0)
        if (!exact) {
          margin <- .relative_target_margin(d, values, ces)
          expected <- .relative_target_dense(d, values, margin)
          add(kernel, case, "normalized_share_density", field("relative_share_loglik"), expected[["share"]])
          add(kernel, case, "normalized_joint_density", field("relative_joint_loglik"), expected[["joint"]])
          add(kernel, case, "supply_prediction", max(abs(log(field("pred_marginNormInv") + 1e-9) - margin)), 0, 1e-10)
          if (d$supply_model == 2L && d$N_market_year == 1L) {
            id <- grep("^s0_logit", names(u));u1 <- u;u2 <- u;u1[id] <- -2.1;u2[id] <- -.5
            delta <- native$log_density(u2, propto = FALSE, jacobian = TRUE) - native$log_density(u1, propto = FALSE, jacobian = TRUE)
            prior <- dnorm(-.5, d$s0_prior_mean, d$s0_prior_sd, log = TRUE) - dnorm(-2.1, d$s0_prior_mean, d$s0_prior_sd, log = TRUE)
            add(kernel, case, "MC_outside_prior_factorization", delta, prior)
          }
        } else {
          expected_share <- .relative_target_dense(d, values, rep(0, d$N))[["share"]]
          expected_joint <- expected_share
          for (m in seq_len(d$N_market)) {
            first <- d$market_mky_order[d$market_mky_start[m]];unit <- d$mky_to_partition_unit[first]
            components <- seq.int(d$partition_unit_start[unit], length.out = d$partition_unit_count[unit])
            scores <- vapply(components, function(cc) {
              means <- .relative_target_exact_margin(d, values, ces, fixture$arguments$structure, d$partition_component_k[cc])
              z <- .relative_target_dense(d, values, means, which(d$market == m))
              d$partition_prior_logprob[cc] + z[["joint"]] - z[["share"]]
            }, numeric(1))
            expected_joint <- expected_joint + .relative_target_lse(scores)
            probability <- exp(scores - .relative_target_lse(scores))
            add(kernel, case, paste0("posterior_k_market", m),
              max(abs(field("posterior_component_prob")[components] - probability)), 0, 1e-9)
          }
          add(kernel, case, "normalized_share_density", field("relative_share_loglik"), expected_share)
          add(kernel, case, "normalized_joint_mixture", field("relative_joint_loglik"), expected_joint)
        }
        for (t in seq_len(d$N_market_year)) {
          ii <- which(d$market_year_idx == t);v <- field("relative_share_predictor")[ii]
          prob <- exp(v - .relative_target_lse(v))
          add(kernel, case, paste0("normalized_shares_cell", t), max(abs(field("pred_shareIn")[ii] - prob)), 0, 1e-11)
          if (d$quality_parameterization == 1L)
            add(kernel, case, paste0("derived_intercept_cell", t), field("derived_market_year_intercept")[t],
              -qlogis(field("s0_prob")[t]) - .relative_target_lse(v), 1e-10)
        }
      }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w));invokeRestart("muffleWarning")}),
      error = function(e) add(kernel, case, "evaluation_error", error = conditionMessage(e)))
    }
  }
  table <- do.call(rbind, results)
  write.csv(table, file.path(output_dir, paste0(prefix, ".csv")), row.names = FALSE)
  result <- list(created = Sys.time(), results = table, warnings = unique(warnings),
    source_hashes = sources, fixture_kind = "fixed unconstrained points; no HMC or compilation",
    passed = !any(table$status != "passed"))
  saveRDS(result, file.path(output_dir, paste0(prefix, ".rds")))
  if (stop_on_failure && !result$passed) stop(sum(table$status != "passed"), " native target checks failed; see retained result table.")
  invisible(result)
}

# Numerically delicate invariants use analytic oracles, rather than subtracting
# almost singular dense covariance matrices. Still fixed-point work only.
relative_share_target_invariants <- function(artifact_dir, output_dir = artifact_dir,
    prefix = "relative-native-invariants", stop_on_failure = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  add <- function(kernel, test, actual, expected, tolerance) {
    rows[[length(rows) + 1L]] <<- data.frame(kernel, test, actual, expected,
      absolute_error = abs(actual - expected), tolerance,
      status = if (is.finite(actual) && abs(actual - expected) <= tolerance) "passed" else "failed")
  }
  for (kernel in c("logit_supply", "ces_supply")) {
    metadata <- readRDS(file.path(artifact_dir, paste0(kernel, "-nutpie-artifact.rds")))
    stopifnot(identical(.relative_target_hash(.relative_target_internal(".cmdstan_kernel_file")(kernel)), metadata$source_hash),
      identical(.relative_target_hash(metadata$lib_path), metadata$lib_hash))
    alias <- file.path(output_dir, paste0("invariant_", kernel, "_", substr(metadata$lib_hash, 1L, 16L), ".so"))
    if (!file.exists(alias)) stopifnot(file.copy(metadata$lib_path, alias))
    ces <- startsWith(kernel, "ces")
    args <- list(prices = rep(1, 4L), quantities = rep(1, 4L), margins = rep(.2, 4L),
      market = rep("m", 4L), year = rep(c("t1", "t2"), each = 2L), group = rep(c("g1", "g2"), 2L),
      supply_model = 2L, use_cutoff = FALSE, use_rho = FALSE, use_hmt = FALSE, priors = list(),
      use_alpha_market_pool = FALSE, outside_price = 1, likelihood_mode = 2L,
      quality_scope = "market", quality_parameterization = "collapsed", quality_scale_parameterization = "legacy",
      share_basis = "relative", share_likelihood = "gaussian", margin_likelihood = "gaussian")
    d <- do.call(.relative_target_internal(if (ces) "prepare_ces_stan_data" else "prepare_stan_data"), args)
    delta <- 2^-34
    # Any common cell level is irrelevant. This exact binary representative
    # avoids rounding the tiny contrast increment through a softmax/log cycle.
    d$log_shareIn <- c(.25, -.25, .25 + delta, -.25 - delta)
    file <- file.path(output_dir, paste0(prefix, "-thin-", kernel, ".json"))
    cmdstanr::write_stan_json(d, file)
    # The default JSON writer retains 15 significant digits. That changes a
    # 2^-34 increment on a .25 baseline enough to confound this exact oracle.
    # Preserve the binary fixture's doubles, without altering production I/O.
    json <- paste(readLines(file, warn = FALSE), collapse = "\n")
    json <- sub('"log_shareIn": \\[[^]]*\\]',
      paste0('"log_shareIn": [', paste(sprintf("%.17g", d$log_shareIn), collapse = ", "), "]"), json)
    writeLines(json, file)
    serialized <- as.numeric(jsonlite::fromJSON(file)$log_shareIn)
    stopifnot(identical(serialized, as.numeric(d$log_shareIn)))
    increment <- serialized[3L] - serialized[1L]
    add(kernel, "thin_noise_serialized_increment", increment / delta, 1, 0)
    native <- suppressWarnings(bridgestan::StanModel$new(alias, file, 137L))
    u <- .relative_target_u(native);u[grepl("^log_sigma_share$", names(u))] <- log(delta)
    u[grepl("^sigma_b_strat$", names(u))] <- 0
    values <- .relative_target_values(native, u)
    expected <- dnorm(.5 + increment, 0, sqrt(2 + delta^2), log = TRUE) +
      dnorm(-increment, 0, delta, log = TRUE)
    add(kernel, "thin_noise_normalized_density", .relative_target_field(values, "relative_share_loglik"), expected, 1e-7)
    gradient <- native$log_density_gradient(u, propto = FALSE, jacobian = TRUE)$gradient
    add(kernel, "thin_noise_finite_gradient", as.numeric(length(gradient) == length(u) && all(is.finite(gradient))), 1, 0)
    if (ces) {
      f <- .relative_target_fixture(kernel, quality = "sampled", rho = TRUE, conduct = 1L)
      old <- f$data;lambda <- 2.7;new <- old
      new$price <- old$price * lambda;new$outside_price <- old$outside_price * lambda
      new$inside_price_index <- old$inside_price_index * lambda
      path1 <- file.path(output_dir, paste0(prefix, "-ces-scale-original.json"))
      path2 <- file.path(output_dir, paste0(prefix, "-ces-scale-transformed.json"))
      cmdstanr::write_stan_json(old, path1);cmdstanr::write_stan_json(new, path2)
      first <- suppressWarnings(bridgestan::StanModel$new(alias, path1, 137L))
      second <- suppressWarnings(bridgestan::StanModel$new(alias, path2, 137L))
      u <- .relative_target_u(first);a <- .relative_target_values(first, u);b <- .relative_target_values(second, u)
      add(kernel, "CES_price_unit_target", second$log_density(u, propto = FALSE, jacobian = TRUE),
        first$log_density(u, propto = FALSE, jacobian = TRUE), 1e-9)
      gamma <- .relative_target_field(a, "gamma_structural_market")
      add(kernel, "CES_price_unit_revenue_outside", .relative_target_field(b, "R0_prob"),
        .relative_target_field(a, "R0_prob"), 1e-12)
      add(kernel, "CES_price_unit_derived_intercept", .relative_target_field(b, "derived_market_year_intercept"),
        .relative_target_field(a, "derived_market_year_intercept") - gamma * log(lambda), 1e-10)
      q <- .relative_target_field(a, "b_quality_all")[old$obs_to_quality]
      meanval1 <- .relative_target_field(a, "derived_market_year_intercept") + q + gamma * mean(log(old$price))
      meanval2 <- .relative_target_field(b, "derived_market_year_intercept") + q + gamma * mean(log(new$price))
      add(kernel, "CES_price_unit_physical_meanval", max(abs(meanval1 - meanval2)), 0, 1e-10)
    }
  }
  table <- do.call(rbind, rows)
  write.csv(table, file.path(output_dir, paste0(prefix, ".csv")), row.names = FALSE)
  saveRDS(list(created = Sys.time(), results = table, source = "fixed points; no HMC or compilation"),
    file.path(output_dir, paste0(prefix, ".rds")))
  if (stop_on_failure && any(table$status != "passed")) stop("Native invariant checks failed; see retained table.")
  invisible(table)
}

# Old native libraries are density references only. Their posterior draws are
# never reused as relative-share fits. No historical artifact is mutated.
relative_share_absolute_preservation <- function(artifact_dir, reference_fits,
    output_dir = artifact_dir, prefix = "absolute-native-preservation",
    points = 10L, stop_on_failure = TRUE) {
  stopifnot(is.character(reference_fits), !is.null(names(reference_fits)), points > 0L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  for (kernel in names(reference_fits)) {
    saved <- readRDS(reference_fits[[kernel]])
    old <- saved@stanfit$metadata
    current <- readRDS(file.path(artifact_dir, paste0(kernel, "-nutpie-artifact.rds")))
    stopifnot(identical(.relative_target_hash(old$lib_path), old$lib_hash),
      identical(.relative_target_hash(current$lib_path), current$lib_hash),
      identical(.relative_target_hash(.relative_target_internal(".cmdstan_kernel_file")(kernel)), current$source_hash))
    alias <- function(metadata, role) {
      path <- file.path(output_dir, paste0("absolute_", role, "_", kernel, "_", substr(metadata$lib_hash, 1L, 16L), ".so"))
      if (!file.exists(path)) stopifnot(file.copy(metadata$lib_path, path))
      stopifnot(identical(.relative_target_hash(path), metadata$lib_hash))
      path
    }
    data_old <- saved@stan_data
    if (!is.null(data_old$share_basis) && data_old$share_basis != 2L)
      stop("A reference fit must have the absolute-share target")
    data_new <- data_old;data_new$share_basis <- 2L
    path_old <- file.path(output_dir, paste0(prefix, "-", kernel, "-old.json"))
    path_new <- file.path(output_dir, paste0(prefix, "-", kernel, "-new.json"))
    cmdstanr::write_stan_json(data_old, path_old);cmdstanr::write_stan_json(data_new, path_new)
    model_old <- suppressWarnings(bridgestan::StanModel$new(alias(old, "old"), path_old, 137L))
    model_new <- suppressWarnings(bridgestan::StanModel$new(alias(current, "new"), path_new, 137L))
    unames_old <- unlist(model_old$param_unc_names());unames_new <- unlist(model_new$param_unc_names())
    stopifnot(identical(unames_old, unames_new),
      identical(model_old$param_names(), model_new$param_names()))
    for (point in seq_len(points)) {
      u <- setNames(.2 * sin(seq_along(unames_old) * point), unames_old)
      u[grepl("^rho_raw", names(u))] <- atanh(.61 / .95)
      expected <- model_old$log_density(u, propto = FALSE, jacobian = TRUE)
      actual <- model_new$log_density(u, propto = FALSE, jacobian = TRUE)
      results[[length(results) + 1L]] <- data.frame(kernel, point, actual, expected,
        absolute_error = abs(actual - expected), tolerance = 1e-8,
        status = if (is.finite(actual) && abs(actual - expected) <= 1e-8) "passed" else "failed",
        old_source_hash = old$source_hash, new_source_hash = current$source_hash)
    }
  }
  table <- do.call(rbind, results)
  write.csv(table, file.path(output_dir, paste0(prefix, ".csv")), row.names = FALSE)
  saveRDS(list(created = Sys.time(), references = reference_fits, results = table,
    source = "old native normalized target versus explicit absolute-share branch; no HMC"),
    file.path(output_dir, paste0(prefix, ".rds")))
  if (stop_on_failure && any(table$status != "passed")) stop("Absolute-share preservation failed; retained table available.")
  invisible(table)
}

# Focused native flow checks for the additional allowed relative-share scale
# configurations. These complement, rather than replace, the baseline matrix.
relative_share_target_scale_checks <- function(artifact_dir, output_dir = artifact_dir,
    prefix = "relative-native-scales", stop_on_failure = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list(); serial <- 0L
  add <- function(kernel, case, test, actual = NA_real_, expected = NA_real_, tolerance = 1e-8, error = NULL) {
    passed <- is.null(error) && is.finite(actual) && is.finite(expected) && abs(actual - expected) <= tolerance
    rows[[length(rows) + 1L]] <<- data.frame(kernel, case, test, actual, expected,
      absolute_error = abs(actual - expected), tolerance,
      status = if (passed) "passed" else "failed", error = if (is.null(error)) NA_character_ else error)
  }
  remake <- function(fixture, kernel, ...) {
    fixture$arguments <- utils::modifyList(fixture$arguments, list(...))
    fixture$data <- do.call(.relative_target_internal(if (startsWith(kernel, "ces"))
      "prepare_ces_stan_data" else "prepare_stan_data"), fixture$arguments)
    fixture
  }
  for (kernel in c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg")) {
    metadata <- readRDS(file.path(artifact_dir, paste0(kernel, "-nutpie-artifact.rds")))
    stopifnot(identical(.relative_target_hash(.relative_target_internal(".cmdstan_kernel_file")(kernel)), metadata$source_hash),
      identical(.relative_target_hash(metadata$lib_path), metadata$lib_hash))
    alias <- file.path(output_dir, paste0("scales_", kernel, "_", substr(metadata$lib_hash, 1L, 16L), ".so"))
    if (!file.exists(alias)) stopifnot(file.copy(metadata$lib_path, alias))
    stopifnot(identical(.relative_target_hash(alias), metadata$lib_hash))
    model <- function(d) {
      serial <<- serial + 1L
      path <- file.path(output_dir, paste0(prefix, "-data-", serial, ".json"))
      cmdstanr::write_stan_json(d, path)
      suppressWarnings(bridgestan::StanModel$new(alias, path, 137L))
    }
    exact <- grepl("stackelberg$", kernel); ces <- startsWith(kernel, "ces")
    fixtures <- list()
    for (rho in c(FALSE, TRUE)) for (quality in c("sampled", "collapsed")) {
      scopes <- if (exact) "market" else c("market", "global")
      for (scope in scopes) for (structure in if (exact) c("fringe", "stackelberg") else "simultaneous")
        for (mixture in if (exact) c(FALSE, TRUE) else FALSE) {
          fixture <- .relative_target_fixture(kernel, markets = 4L, periods = 4L,
            quality = quality, rho = rho, scope = scope, conduct = 3L,
            structure = structure, mixture = mixture, ragged = TRUE)
          fixtures[[length(fixtures) + 1L]] <- remake(fixture, kernel, use_market_resid_scale = TRUE)
        }
    }
    for (fixture in fixtures) {
      case <- paste0("market_residual_", fixture$id); d <- fixture$data
      tryCatch({
        native <- model(d); u <- .relative_target_u(native)
        u[grepl("^sigma_log_sigma_share_market", names(u))] <- log(.4)
        u[grepl("^sigma_log_sigma_margin_market", names(u))] <- log(.25)
        u[grepl("^log_sigma_share_market_raw", names(u))] <- c(-.7, -.1, .4, 1.1)
        u[grepl("^log_sigma_margin_market_raw", names(u))] <- c(.8, -.4, .3, -.6)
        values <- .relative_target_values(native, u); field <- function(name) .relative_target_field(values, name)
        gradient <- native$log_density_gradient(u, propto = FALSE, jacobian = TRUE)$gradient
        add(kernel, case, "finite_normalized_target", as.numeric(is.finite(native$log_density(u,
          propto = FALSE, jacobian = TRUE))), 1, 0)
        add(kernel, case, "finite_gradient", as.numeric(length(gradient) == length(u) && all(is.finite(gradient))), 1, 0)
        add(kernel, case, "heterogeneous_share_scales", max(abs(field("sigma_share_market") - .23 * exp(.4 * c(-.7, -.1, .4, 1.1)))), 0, 1e-12)
        add(kernel, case, "heterogeneous_margin_scales", max(abs(field("sigma_margin_market") - .31 * exp(.25 * c(.8, -.4, .3, -.6)))), 0, 1e-12)
        if (!exact) {
          margin <- .relative_target_margin(d, values, ces)
          expected <- .relative_target_dense(d, values, margin)
          add(kernel, case, "normalized_share_density", field("relative_share_loglik"), expected[["share"]])
          add(kernel, case, "normalized_joint_density", field("relative_joint_loglik"), expected[["joint"]])
        } else {
          expected_share <- .relative_target_dense(d, values, rep(0, d$N))[["share"]]
          expected_joint <- expected_share
          for (m in seq_len(d$N_market)) {
            first <- d$market_mky_order[d$market_mky_start[m]]; unit <- d$mky_to_partition_unit[first]
            components <- seq.int(d$partition_unit_start[unit], length.out = d$partition_unit_count[unit])
            scores <- vapply(components, function(cc) {
              means <- .relative_target_exact_margin(d, values, ces, fixture$arguments$structure,
                d$partition_component_k[cc])
              z <- .relative_target_dense(d, values, means, which(d$market == m))
              d$partition_prior_logprob[cc] + z[["joint"]] - z[["share"]]
            }, numeric(1))
            expected_joint <- expected_joint + .relative_target_lse(scores)
            add(kernel, case, paste0("posterior_k_market", m), max(abs(field("posterior_component_prob")[components] -
              exp(scores - .relative_target_lse(scores)))), 0, 1e-9)
          }
          add(kernel, case, "normalized_share_density", field("relative_share_loglik"), expected_share)
          add(kernel, case, "normalized_joint_mixture", field("relative_joint_loglik"), expected_joint)
        }
      }, error = function(e) add(kernel, case, "evaluation_error", error = conditionMessage(e)))
    }
    # Both coordinates represent the same physical sigma and tau. With
    # unconstrained coordinates, |d(log tau,log sigma)/d(log T,logit lambda)|=1/2.
    # Test the normalized target and analytic gradient chain rule explicitly.
    for (layout in c(1L, 4L)) for (rho in c(FALSE, TRUE))
      for (scope in if (exact) "market" else c("market", "global"))
        for (structure in if (exact) c("fringe", "stackelberg") else "simultaneous")
          for (mixture in if (exact) c(FALSE, TRUE) else FALSE) {
            fixture <- .relative_target_fixture(kernel, markets = layout, periods = layout,
              quality = "collapsed", rho = rho, scope = scope, conduct = 3L,
              structure = structure, mixture = mixture, ragged = layout > 1L)
            case <- paste0("total_fraction_", fixture$id)
            tryCatch({
              legacy <- model(fixture$data)
              transformed <- model(remake(fixture, kernel, quality_scale_parameterization = "total_fraction")$data)
              u1 <- .relative_target_u(legacy); u2 <- u1
              stopifnot(identical(legacy$param_unc_names(), transformed$param_unc_names()))
              iq <- match("sigma_b_strat", names(u1)); is <- match("log_sigma_share", names(u1))
              tau <- exp(u1[iq]); sigma <- exp(u1[is]); lambda <- sigma^2 / (sigma^2 + tau^2)
              u2[iq] <- log(sqrt(sigma^2 + tau^2)); u2[is] <- qlogis(lambda)
              v1 <- .relative_target_values(legacy, u1); v2 <- .relative_target_values(transformed, u2)
              g1 <- legacy$log_density_gradient(u1, propto = FALSE, jacobian = TRUE)
              g2 <- transformed$log_density_gradient(u2, propto = FALSE, jacobian = TRUE)
              add(kernel, case, "physical_share_scale", .relative_target_field(v2, "sigma_share"), sigma, 1e-12)
              add(kernel, case, "physical_quality_scale", .relative_target_field(v2, "sigma_b_strat_economic"), tau, 1e-12)
              for (name in c("relative_share_loglik", "relative_joint_loglik"))
                add(kernel, case, paste0(name, "_coordinate_invariant"), .relative_target_field(v2, name),
                  .relative_target_field(v1, name))
              add(kernel, case, "normalized_unconstrained_density_jacobian",
                transformed$log_density(u2, propto = FALSE, jacobian = TRUE),
                legacy$log_density(u1, propto = FALSE, jacobian = TRUE) - log(2))
              expected_gradient <- g1$gradient
              expected_gradient[iq] <- g1$gradient[iq] + g1$gradient[is]
              expected_gradient[is] <- -.5 * lambda * g1$gradient[iq] + .5 * (1 - lambda) * g1$gradient[is]
              add(kernel, case, "gradient_coordinate_chain_rule", max(abs(g2$gradient - expected_gradient)), 0, 1e-7)
              add(kernel, case, "finite_transformed_gradient", as.numeric(length(g2$gradient) == length(u2) && all(is.finite(g2$gradient))), 1, 0)
            }, error = function(e) add(kernel, case, "evaluation_error", error = conditionMessage(e)))
          }
  }
  table <- do.call(rbind, rows)
  write.csv(table, file.path(output_dir, paste0(prefix, ".csv")), row.names = FALSE)
  saveRDS(list(created = Sys.time(), results = table, source = "fixed native points; no HMC or compilation"),
    file.path(output_dir, paste0(prefix, ".rds")))
  if (stop_on_failure && any(table$status != "passed")) stop("Native scale checks failed; see retained table.")
  invisible(table)
}
