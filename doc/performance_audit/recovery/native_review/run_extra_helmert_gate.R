suppressPackageStartupMessages({
  library(antitrustBayes)
  library(rstan)
  library(digest)
})

role <- Sys.getenv("HELMERT_ROLE", "candidate")
if (!role %in% c("candidate", "old"))
  stop("HELMERT_ROLE must be candidate or old")
root <- Sys.getenv("HELMERT_REVIEW_ROOT",
  "doc/performance_audit/recovery/native_review")
artifact_dir <- Sys.getenv("HELMERT_ARTIFACT_DIR",
  file.path(root, ".artifacts"))
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(artifact_dir, paste0("extra_helmert_", role, ".rds"))
raw_file <- file.path(artifact_dir, "extra_helmert_candidate_raw.rds")
failure_file <- file.path(root, paste0("extra_helmert_", role, "_failures.csv"))
case_root <- file.path(artifact_dir, paste0("cases_", role))
dir.create(case_root, recursive = TRUE, showWarnings = FALSE)
force <- tolower(Sys.getenv("HELMERT_FORCE", "false")) %in%
  c("1", "true", "yes")

source(file.path(root, "cache_helpers.R"), local = TRUE)
ns <- asNamespace("antitrustBayes")
mods <- get("stanmodels", ns)
invisible(capture.output(sys.source(
  file.path(root, "extra_helmert_fixtures.R"), envir = environment())))

kernel_for <- function(demand, structure, partition_prior = "cutoff") {
  exact <- identical(structure, "stackelberg") ||
    (identical(structure, "fringe") && identical(partition_prior, "ordered_k"))
  paste0(demand, if (exact) "_stackelberg" else "_supply")
}

param_names <- function(sf) {
  p <- sf@sim$pars_oi
  if (is.null(p)) character() else as.character(p)
}

set_nonzero <- function(cp, sf, nm, value) {
  pnames <- param_names(sf)
  if (!(nm %in% pnames) || is.null(cp[[nm]]) || !length(cp[[nm]])) return(cp)
  cp[[nm]][] <- rep(value, length.out = length(cp[[nm]]))
  cp
}

nonzero_pars <- function(cp, sf) {
  for (nm in c("mu_b_market_vec")) cp <- set_nonzero(cp, sf, nm, .15)
  for (nm in c("sigma_b_market_vec", "sigma_b_strat", "sigma_log_a_vec",
               "sigma_s0_pool_vec", "sigma_log_sigma_share_market",
               "sigma_log_sigma_margin_market"))
    cp <- set_nonzero(cp, sf, nm, .7)
  for (nm in c("b_market_raw", "r_market_a_raw", "s0_logit", "s0_pool_raw",
               "eta_s0_pool", "s0_market_raw", "eta_s0_market", "c_s0",
               "c_s0_raw", "mu_s0_pool_vec", "mu_log_a_free",
               "log_a_market_pool", "log_kappa_vec", "log_kappa_market",
               "log_sigma_share_market_raw", "log_sigma_margin_market_raw",
               "beta_s0")) cp <- set_nonzero(cp, sf, nm, .08)
  for (nm in c("log_sigma_share", "log_sigma_margin"))
    cp <- set_nonzero(cp, sf, nm, -.15)
  cp <- set_nonzero(cp, sf, "rho_raw", .18)
  cp <- set_nonzero(cp, sf, "cutoff_share", .35)
  cp <- set_nonzero(cp, sf, "bargaining_unit", .4)
  cp <- set_nonzero(cp, sf, "nu_minus_2", 3)
  if ("b_group_centered" %in% param_names(sf) &&
      length(cp$b_group_centered))
    cp$b_group_centered[] <- seq(-.3, .3, length.out = length(cp$b_group_centered))
  cp
}

flat_numeric <- function(x) {
  y <- unlist(x, recursive = TRUE, use.names = TRUE)
  y <- y[vapply(y, is.numeric, logical(1))]
  as.numeric(y) |> setNames(names(y))
}

flat_q <- function(q) {
  y <- unlist(q, recursive = TRUE, use.names = TRUE)
  y <- y[vapply(y, is.numeric, logical(1))]
  as.numeric(y) |> setNames(names(y))
}

representative_fd <- function(u, lp_fun) {
  n <- length(u)
  ix <- sort(unique(c(1L, ceiling(n / 2), n)))
  h <- 1e-5
  val <- vapply(ix, function(i) {
    up <- u; um <- u
    up[i] <- up[i] + h; um[i] <- um[i] - h
    (lp_fun(up) - lp_fun(um)) / (2 * h)
  }, numeric(1))
  list(index = ix, value = val)
}

hash_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

source_hashes <- list(
  runner = hash_file(file.path(root, "run_extra_helmert_gate.R")),
  fixture = hash_file(file.path(root, "extra_helmert_fixtures.R")),
  cache_helper = hash_file(file.path(root, "cache_helpers.R"))
)
settings <- list(chains = 1L, iter = 1L, warmup = 0L,
                 algorithm = "Fixed_param", init = 0,
                 refresh = 0L, fd_step = 1e-5)
seeds <- list(target = 58003L, generated_quantities = 58005L)

package_dll_manifest <- local({
  pkg <- find.package("antitrustBayes")
  files <- list.files(file.path(pkg, "libs"), full.names = TRUE,
                      pattern = "\\.(so|dll|dylib)$", ignore.case = TRUE)
  if (!length(files)) {
    list(files = character(), md5 = character(), hash = NA_character_)
  } else {
    md5 <- as.character(tools::md5sum(files))
    names(md5) <- basename(files)
    list(files = names(md5), md5 = unname(md5),
         hash = native_review_hash(md5))
  }
})

prepare_case <- function(spec) {
  dat <- do.call(make_extra_data, spec)
  prior <- if (is.null(spec$partition_prior)) "cutoff" else spec$partition_prior
  kernel <- kernel_for(spec$demand, spec$structure, prior)
  model <- mods[[kernel]]
  if (is.null(model)) stop("missing stan model ", kernel)
  loaded_dso <- tryCatch(model@dso@dso_filename,
                          error = function(e) NA_character_)
  loaded_dso_hash <- if (length(loaded_dso) && file.exists(loaded_dso))
    hash_file(loaded_dso) else NA_character_
  list(dat = dat, model = model, kernel = kernel,
       metadata = list(
         input_hash = native_review_hash(dat),
         kernel_code_hash = digest::digest(model@model_code, algo = "sha256",
                                           serialize = FALSE),
         package_dll_hash = package_dll_manifest$hash,
         loaded_dso_hash = loaded_dso_hash,
         source_hashes = source_hashes,
         settings = settings, seeds = seeds,
         package_version = as.character(packageVersion("antitrustBayes"))))
}

one <- function(case_name, spec, prepared, supplied_u = NULL) {
  started <- proc.time()[["elapsed"]]
  dat <- prepared$dat
  model <- prepared$model
  kernel <- prepared$kernel
  meta <- prepared$metadata
  tryCatch({
    sf <- suppressWarnings(sampling(model, data = dat, chains = settings$chains,
      iter = settings$iter, warmup = settings$warmup,
      algorithm = settings$algorithm, init = settings$init,
      refresh = settings$refresh, seed = seeds$target))
    nu <- get_num_upars(sf)
    schema <- list(pars_oi = param_names(sf), dims_oi = sf@sim$dims_oi,
                   num_upars = nu)
    if (is.null(supplied_u)) {
      cp <- constrain_pars(sf, rep(0, nu))
      cp <- nonzero_pars(cp, sf)
      u <- unconstrain_pars(sf, cp)
    } else {
      u <- as.numeric(supplied_u)
      if (length(u) != nu)
        stop("raw parameter dimension differs: stored=", length(u),
             " current=", nu)
    }
    lp_fun <- function(z, jac = TRUE)
      log_prob(sf, z, adjust_transform = jac, gradient = FALSE)
    lp <- lp_fun(u, TRUE)
    lp_raw <- lp_fun(u, FALSE)
    gr <- grad_log_prob(sf, u, adjust_transform = TRUE)
    gr_raw <- grad_log_prob(sf, u, adjust_transform = FALSE)
    fd <- representative_fd(u, function(z) lp_fun(z, TRUE))
    cp_at <- constrain_pars(sf, u)
    tp <- flat_numeric(cp_at)
    sf_q <- suppressWarnings(sampling(model, data = dat, chains = settings$chains,
      iter = settings$iter, warmup = settings$warmup,
      algorithm = settings$algorithm,
      init = function(...) cp_at, refresh = settings$refresh,
      seed = seeds$generated_quantities))
    extracted <- extract(sf_q, permuted = TRUE, inc_warmup = FALSE)
    q <- flat_q(extracted)
    quality <- as.numeric(extracted[["b_quality_all"]])
    expected_quality_length <- if (dat$quality_parameterization == 1L)
      dat$N_group else 0L
    if (length(quality) != expected_quality_length)
      stop("quality output length differs from the selected parameterization")
    c(list(status = "ok", case = case_name, kernel = kernel,
      demand = spec$demand, structure = spec$structure,
      supply_model = spec$supply_model, N = dat$N,
      N_market = dat$N_market, N_market_year = dat$N_market_year,
      N_group = dat$N_group, K_s0 = dat$K_s0, schema = schema, u = u,
      lp = as.numeric(lp), lp_raw = as.numeric(lp_raw), grad = as.numeric(gr),
      grad_raw = as.numeric(gr_raw), fd = fd, tp = tp, q = q,
      quality = as.numeric(quality), elapsed = proc.time()[["elapsed"]] - started),
      meta)
  }, error = function(e) c(list(status = "error", case = case_name,
                                kernel = kernel, error = conditionMessage(e),
                                elapsed = proc.time()[["elapsed"]] - started),
                            meta))
}

raw <- if (identical(role, "old")) {
  if (!file.exists(raw_file)) stop("missing candidate raw coordinates: ", raw_file)
  readRDS(raw_file)
} else if (file.exists(raw_file) && !force) {
  readRDS(raw_file)
} else {
  list(schema = 2L, role = "candidate", entries = list())
}
if (identical(role, "old") &&
    (!identical(raw$schema, 2L) || !identical(raw$role, "candidate") ||
     !is.list(raw$entries)))
  stop("candidate raw checkpoint lacks the required provenance schema")

res <- list()
for (nm in names(extra_cases)) {
  cat(role, nm, "\n", sep = " ")
  spec <- extra_cases[[nm]]
  prepared <- tryCatch(prepare_case(spec), error = function(e)
    list(error = conditionMessage(e), kernel = kernel_for(spec$demand,
         spec$structure), metadata = list(input_hash = NA_character_,
         kernel_code_hash = NA_character_, package_dll_hash = NA_character_,
         loaded_dso_hash = NA_character_, source_hashes = source_hashes,
         settings = settings, seeds = seeds,
         package_version = as.character(packageVersion("antitrustBayes")))))
  meta <- prepared$metadata
  entry <- if (identical(role, "old")) raw$entries[[nm]] else NULL
  supplied <- if (identical(role, "old") && is.list(entry)) entry$u else NULL
  raw_hash <- if (identical(role, "old") && !is.null(supplied))
    native_review_hash(supplied) else NA_character_
  cache_spec <- list(schema = 2L, role = role, case = nm,
    fixture = spec, kernel = prepared$kernel, input_hash = meta$input_hash,
    kernel_code_hash = meta$kernel_code_hash,
    package_dll_hash = meta$package_dll_hash,
    loaded_dso_hash = meta$loaded_dso_hash, source_hashes = source_hashes,
    raw_points_hash = raw_hash, settings = settings, seeds = seeds,
    coordinate_initializer = "nonzero_pars:v2")
  callback <- function() {
    if (!is.null(prepared$error))
      return(list(status = "error", case = nm, kernel = prepared$kernel,
                  error = prepared$error, input_hash = meta$input_hash,
                  kernel_hash = meta$kernel_code_hash,
                  package_version = meta$package_version, elapsed = 0))
    if (identical(role, "old")) {
      guard <- tryCatch({
        native_review_validate_raw(entry, nm, meta$input_hash, source_hashes,
                                   settings, seeds)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(guard))
        return(list(status = "error", case = nm, kernel = prepared$kernel,
                    error = guard, input_hash = meta$input_hash,
                    kernel_hash = meta$kernel_code_hash,
                    package_version = meta$package_version, elapsed = 0))
    }
    one(nm, spec, prepared, supplied)
  }
  cached <- native_review_cache(file.path(case_root, nm), cache_spec,
                                callback, force = force)
  ans <- cached$result
  if (is.null(ans$case)) ans$case <- nm
  if (is.null(ans$kernel)) ans$kernel <- prepared$kernel
  if (is.null(ans$elapsed)) ans$elapsed <- 0
  if (is.null(ans$error) && !identical(ans$status, "ok"))
    ans$error <- "native review attempt failed without an error message"
  ans$cache_hash <- cached$cache_hash
  ans$cache_generation <- cached$generation
  ans$cache_attempt <- cached$attempt
  ans$cache_reused <- cached$reused
  ans$cache_state <- cached$state
  res[[nm]] <- ans
  if (identical(role, "candidate") && identical(ans$status, "ok")) {
    raw$entries[[nm]] <- list(schema = 2L, role = "candidate", case = nm,
      u = ans$u, u_hash = native_review_hash(ans$u),
      input_hash = meta$input_hash, kernel_code_hash = meta$kernel_code_hash,
      source_hashes = source_hashes, settings = settings, seeds = seeds)
    native_review_atomic_save(raw, raw_file)
  }
}

if (identical(role, "candidate")) native_review_atomic_save(raw, raw_file)
native_review_atomic_save(list(role = role, results = res,
                               cases = names(extra_cases)), out_file)

rows <- lapply(res, function(x) {
  if (!identical(x$status, "ok"))
    return(data.frame(case = x$case, status = x$status, error = x$error,
                      cache_hash = x$cache_hash,
                      stringsAsFactors = FALSE))
  data.frame(case = x$case, status = x$status, kernel = x$kernel, N = x$N,
    N_market = x$N_market, N_market_year = x$N_market_year,
    N_group = x$N_group, K_s0 = x$K_s0, num_upars = x$schema$num_upars,
    input_hash = x$input_hash, kernel_code_hash = x$kernel_code_hash,
    package_dll_hash = x$package_dll_hash,
    loaded_dso_hash = x$loaded_dso_hash,
    lp = x$lp, lp_raw = x$lp_raw, grad_l2 = sqrt(sum(x$grad^2)),
    grad_raw_l2 = sqrt(sum(x$grad_raw^2)), fd_dim = length(x$fd$index),
    tp_len = length(x$tp), q_len = length(x$q),
    quality_len = length(x$quality), cache_hash = x$cache_hash,
    cache_reused = x$cache_reused, elapsed = x$elapsed,
    stringsAsFactors = FALSE)
})
bind_rows_complete <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    for (column in setdiff(columns, names(row))) row[[column]] <- NA
    row[columns]
  })
  do.call(rbind, rows)
}
write.csv(bind_rows_complete(rows), file.path(root,
  paste0("extra_helmert_", role, ".csv")), row.names = FALSE)
fails <- Filter(function(x) !identical(x$status, "ok"), res)
failure_rows <- lapply(fails, function(x)
  data.frame(case = x$case, status = x$status, error = x$error,
             cache_hash = x$cache_hash, elapsed = x$elapsed,
             stringsAsFactors = FALSE))
if (length(failure_rows)) {
  write.csv(do.call(rbind, failure_rows), failure_file, row.names = FALSE)
} else {
  write.csv(data.frame(case = character(), status = character(),
                        error = character(), cache_hash = character(),
                        elapsed = numeric()), failure_file, row.names = FALSE)
}
