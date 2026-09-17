suppressPackageStartupMessages(library(antitrustBayes))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("baseline", "candidate"))
  stop("supply baseline or candidate")
ns <- asNamespace("antitrustBayes")
approved <- get(".bayes_bridge_approved_hashes", ns)
generation <- if (args == "candidate") "helmert" else "rotation"
expected <- approved[[generation]]
if (is.null(expected)) stop("missing approved generation ", generation)
mods <- get("stanmodels", ns)
actual <- vapply(mods, function(model)
  digest::digest(model@model_code, algo = "sha256", serialize = FALSE),
  character(1))
print(actual)
if (!identical(sort(names(actual)), sort(names(expected))) ||
    !all(actual[names(expected)] == expected))
  stop("loaded kernel generation differs from expected ", generation)
cat("loaded generation PASS:", generation, "\n")
