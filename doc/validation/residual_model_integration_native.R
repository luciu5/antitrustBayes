# Independent Gaussian observation-covariance oracle, with native densities and
# conditional role weights. No HMC, bridge evidence, or model probabilities.
# Frozen relative-share helpers supply fixture geometry and structural FOC means;
# this file independently builds/project the new full covariance.
rmi_native_dense <- function(d, v, margin_mean, rows = seq_len(d$N)) {
  field <- function(name) as.numeric(v[names(v) == name |
    startsWith(names(v), paste0(name, "["))])
  log_normal <- function(residual, covariance) {
    factor <- chol(covariance)
    white <- forwardsolve(t(factor), residual)
    -.5 * (length(residual) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(white^2))
  }
  cell <- d$market_year_idx[rows]
  incidence <- outer(cell, unique(cell), "==") * 1
  basis <- qr.Q(qr(incidence), complete = TRUE)[,
    seq.int(ncol(incidence) + 1L, nrow(incidence)), drop = FALSE]
  m <- d$market[rows]; omega <- field("sigma_demand_shock_market")[m]
  measurement <- field("sigma_measurement")
  share_cov <- diag(omega^2 + measurement^2, length(rows))
  if (d$quality_parameterization == 2L) {
    p <- diag(d$N_quality)
    for (b in seq_len(d$N_quality_block)) {
      ii <- seq.int(d$quality_block_start[b], length.out = d$quality_block_n[b])
      p[ii, ii] <- p[ii, ii] - 1 / length(ii)
    }
    z <- diag(d$N_quality)[d$obs_to_quality[rows], , drop = FALSE]
    share_cov <- share_cov + field("sigma_b_strat_economic")^2 * z %*% p %*% t(z)
  }
  residual <- as.numeric(crossprod(basis,
    d$log_shareIn[rows] - field("relative_share_predictor")[rows]))
  share_cov <- crossprod(basis, share_cov %*% basis)
  share <- log_normal(residual, share_cov)
  observed <- which(d$margin_obs[rows] == 1L)
  if (!length(observed)) return(c(share = share, joint = share))
  sigma_margin <- field("sigma_margin_market")[m[observed]]
  margin_residual <- d$log_marginNormInv[rows[observed]] - margin_mean[rows[observed]]
  cross <- matrix(0, length(rows), length(observed))
  cross[cbind(observed, seq_along(observed))] <-
    field("rho_demand_shock") * omega[observed] * sigma_margin
  cross <- crossprod(basis, cross)
  covariance <- rbind(cbind(share_cov, cross),
    cbind(t(cross), diag(sigma_margin^2, length(observed))))
  c(share = share, joint = log_normal(c(residual, margin_residual), covariance))
}

rmi_native_checks <- function(output, compile = TRUE) {
  source("doc/validation/relative_share_target_checks.R", local = TRUE)
  kernels <- c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg")
  dest <- file.path(output, "native"); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  rows <- list(); artifacts <- list()
  add <- function(kernel, case, test, actual, expected = 0, tolerance = 1e-8, reason = "") {
    rows[[length(rows) + 1L]] <<- data.frame(kernel, case, test, actual, expected, tolerance,
      absolute_error = abs(actual - expected),
      passed = is.finite(actual) && is.finite(expected) && abs(actual - expected) <= tolerance, reason)
    write.csv(do.call(rbind, rows), file.path(dest, "checks.csv"), row.names = FALSE)
  }
  for (kernel in kernels) {
    artifact <- file.path(dest, paste0(kernel, "-artifact.rds"))
    hash <- rmi_hash(rmi_internal(".cmdstan_kernel_file")(kernel), TRUE)
    if (file.exists(artifact)) {
      metadata <- readRDS(artifact)
      stopifnot(identical(metadata$source_hash, hash),
        identical(metadata$lib_hash, rmi_hash(metadata$lib_path, TRUE)))
    } else {
      if (!compile) stop("No compiled native artifact.")
      cat("COMPILE", kernel, format(Sys.time()), "\n"); flush.console()
      metadata <- rmi_internal(".nutpie_model_for_kernel")(kernel)$metadata
      stopifnot(identical(metadata$source_hash, hash))
      saveRDS(metadata, artifact)
    }
    artifacts[[kernel]] <- metadata
    exact <- grepl("stackelberg$", kernel); ces <- startsWith(kernel, "ces")
    cases <- expand.grid(quality = c("sampled", "collapsed"), measurement = c(.04, .25),
      market_scale = c(FALSE, TRUE), residual_model = c("combined", "demand_shock"),
      stringsAsFactors = FALSE)
    cases <- cases[!(cases$residual_model == "combined" & cases$measurement == .25), ]
    for (j in seq_len(nrow(cases))) {
      spec <- cases[j, ]; id <- paste(spec, collapse = "_")
      tryCatch({
        fx <- .relative_target_fixture(kernel, markets = 2L, periods = 2L,
          quality = spec$quality, rho = TRUE, scope = "market",
          conduct = if (exact) 3L else c(1L, 2L, 3L, 4L, 5L)[(j - 1L) %% 5L + 1L],
          structure = if (j %% 2L) "fringe" else "stackelberg", mixture = exact,
          ragged = TRUE, market_side = if (!ces && !exact) "input" else "output")
        args <- fx$arguments; args$share_residual_model <- spec$residual_model
        args$use_market_resid_scale <- spec$market_scale
        d <- do.call(rmi_internal(if (ces) "prepare_ces_stan_data" else "prepare_stan_data"), args)
        file <- file.path(dest, paste0(kernel, "-data.json")); cmdstanr::write_stan_json(d, file)
        native <- bridgestan::StanModel$new(metadata$lib_path, file, 137L)
        u <- .relative_target_u(native)
        u[grepl("^log_sigma_measurement", names(u))] <- log(spec$measurement)
        u[grepl("^rho_raw", names(u))] <- atanh(c(-.5, 0, .6, .25)[(j - 1L) %% 4L + 1L] / .95)
        u[grepl("^log_sigma_share_market_raw", names(u))] <- c(-.7, .7)
        u[grepl("^log_sigma_margin_market_raw", names(u))] <- c(.4, -.4)
        values <- .relative_target_values(native, u)
        field <- function(name) .relative_target_field(values, name)
        add(kernel, id, "measurement_coordinate_dimension",
          sum(grepl("^log_sigma_measurement", names(u))),
          as.numeric(spec$residual_model == "demand_shock"), 0)
        if (spec$residual_model == "combined")
          add(kernel, id, "combined_measurement_exact_zero", field("sigma_measurement"), 0, 0)
        add(kernel, id, "native_finite_normalized_density", as.numeric(is.finite(
          native$log_density(u, propto = FALSE, jacobian = TRUE))), 1, 0)
        grad <- native$log_density_gradient(u, propto = FALSE, jacobian = TRUE)$gradient
        add(kernel, id, "native_finite_gradient", as.numeric(all(is.finite(grad))), 1, 0)
        add(kernel, id, "total_variance", max(abs(field("sigma_share_market")^2 -
          field("sigma_demand_shock_market")^2 - field("sigma_measurement")^2)), 0, 1e-12)
        add(kernel, id, "correlated_demand_covariance", max(abs(field("rho_market") *
          field("sigma_share_market") - field("rho_demand_shock") *
          field("sigma_demand_shock_market"))), 0, 1e-12)
        if (!exact) {
          margin <- .relative_target_margin(d, values, ces)
          expected <- rmi_native_dense(d, values, margin)
          add(kernel, id, "dense_share_density", field("relative_share_loglik"), expected[["share"]], 2e-8)
          add(kernel, id, "dense_joint_density", field("relative_joint_loglik"), expected[["joint"]], 2e-8)
        } else {
          expected_share <- rmi_native_dense(d, values, rep(0, d$N))[["share"]]
          expected_joint <- expected_share
          for (m in seq_len(d$N_market)) {
            first <- d$market_mky_order[d$market_mky_start[m]]
            unit <- d$mky_to_partition_unit[first]
            components <- seq.int(d$partition_unit_start[unit], length.out = d$partition_unit_count[unit])
            scores <- vapply(components, function(cc) {
              means <- .relative_target_exact_margin(d, values, ces, args$structure, d$partition_component_k[cc])
              z <- rmi_native_dense(d, values, means, which(d$market == m))
              d$partition_prior_logprob[cc] + z[["joint"]] - z[["share"]]
            }, numeric(1))
            expected_joint <- expected_joint + .relative_target_lse(scores)
            probability <- exp(scores - .relative_target_lse(scores))
            add(kernel, id, paste0("dense_conditional_role_weights_market", m),
              max(abs(field("posterior_component_prob")[components] - probability)), 0, 2e-8)
          }
          add(kernel, id, "dense_share_density", field("relative_share_loglik"), expected_share, 2e-8)
          add(kernel, id, "dense_joint_role_mixture", field("relative_joint_loglik"), expected_joint, 2e-8)
        }
      }, error = function(e) add(kernel, id, "evaluation_error", Inf, reason = conditionMessage(e)))
    }
  }
  result <- list(checks = do.call(rbind, rows), artifacts = artifacts, at = Sys.time(),
    sampling = FALSE, passed = all(vapply(rows, function(x) x$passed, logical(1))))
  saveRDS(result, file.path(dest, "result.rds")); invisible(result)
}
