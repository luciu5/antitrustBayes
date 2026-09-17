## Preparation-only manifest for the Stan performance audit.
##
## This file deliberately does not compile, sample, call gqs(), or mutate a
## package/library. It records the design and checks that the already-built
## package exposes the four kernels. A later measurement runner may source
## these helpers after adding fixed-parameter target/gradient calls.

stan_performance_design <- function() {
  expand.grid(
    N = c(100L, 1000L, 10000L),
    quality = c("sampled", "collapsed"),
    share_likelihood = c("gaussian", "student_t"),
    kernel = c("logit_supply", "ces_supply",
               "logit_stackelberg", "ces_stackelberg"),
    operation = c("model_target_gradient", "generated_quantities"),
    stringsAsFactors = FALSE
  )
}

stan_performance_environment <- function() {
  package_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(pkg))
  }
  data.frame(
    field = c("R", "rstan", "StanHeaders", "RcppParallel", "logical_cpus"),
    value = c(
      as.character(getRversion()),
      package_version("rstan"),
      package_version("StanHeaders"),
      package_version("RcppParallel"),
      as.character(parallel::detectCores(logical = TRUE))
    ),
    stringsAsFactors = FALSE
  )
}

stan_performance_preflight <- function() {
  if (!requireNamespace("antitrustBayes", quietly = TRUE)) {
    stop("Load an installed antitrustBayes package before running preflight.",
         call. = FALSE)
  }
  models <- get("stanmodels", envir = asNamespace("antitrustBayes"))
  expected <- c("ces_stackelberg", "ces_supply",
                "logit_stackelberg", "logit_supply")
  if (!identical(sort(names(models)), sort(expected))) {
    stop("The installed package does not expose the expected four kernels.",
         call. = FALSE)
  }
  design <- stan_performance_design()
  ## Collapsed quality is implemented only for Gaussian shares. Keeping the
  ## rows in the manifest makes unsupported cells explicit for the later
  ## runner instead of silently dropping them here.
  design$supported <- with(design,
                           quality == "sampled" |
                             share_likelihood == "gaussian")
  list(
    design = design,
    environment = stan_performance_environment(),
    models = expected,
    status = "preflight_only; no compile, target call, gqs, or sampling"
  )
}
