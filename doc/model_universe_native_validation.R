#!/usr/bin/env Rscript
# Public-interface integration check using a frozen, previously fitted BF draw.
# This is not a new recovery experiment or a conduct-universe comparison.
suppressPackageStartupMessages(library(antitrustBayes))
fit_path <- Sys.getenv("ANTITRUSTBAYES_UNIVERSE_FIXTURE",
  "doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_BF_noisy_271828__D.rds")
if (!file.exists(fit_path)) stop("Missing archived native fixture: ", fit_path)
fit <- readRDS(fit_path)
universe <- bayes_model_universe(
  list(small_core = list(fit = "BF", k = 1:2),
       larger_core = list(fit = "BF", k = 3:5)),
  prior_weights = c(small_core = .4, larger_core = .6))
set.seed(20260916)
comparison <- bayes_compare(list(BF = fit), universe = universe,
                            bridge = TRUE, bridge_repetitions = 3)
stopifnot(identical(attr(comparison, "status"), "complete"))
k <- posterior_k(fit)
post <- c(sum(k$posterior[k$k %in% 1:2]),
          sum(k$posterior[k$k %in% 3:5]))
prior <- c(sum(k$prior[k$k %in% 1:2]),
           sum(k$prior[k$k %in% 3:5]))
expected <- c(.4, .6) * post / prior
expected <- expected / sum(expected)
stopifnot(max(abs(comparison$posterior_probability - expected)) < 1e-10,
          max(abs(comparison$posterior_mass - post)) < 1e-10,
          max(abs(comparison$prior_mass - prior)) < 1e-10)
dimensions <- bayes_dimensions(comparison, dimension = "conduct")
stopifnot(abs(sum(dimensions$posterior_probability) - 1) < 1e-10)
out <- Sys.getenv("ANTITRUSTBAYES_UNIVERSE_OUT",
                  "doc/model_universe_native_validation.csv")
comparison$independent_expected_probability <- expected
comparison$probability_error <- abs(comparison$posterior_probability - expected)
utils::write.csv(comparison, out, row.names = FALSE, na = "")
archive <- Sys.getenv("ANTITRUSTBAYES_UNIVERSE_ARCHIVE",
  "doc/statistical_toolkit_runs/20260916/model-universe-native")
dir.create(archive, recursive = TRUE, showWarnings = FALSE)
saveRDS(comparison, file.path(archive, "comparison.rds"))
writeLines(c(paste("fixture_sha256", digest::digest(file = fit_path, algo = "sha256")),
             capture.output(sessionInfo())), file.path(archive, "session.txt"))
cat("Public model-universe native bridge and independent mass checks passed.\n")
