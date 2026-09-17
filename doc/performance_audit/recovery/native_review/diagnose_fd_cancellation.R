# Secondary numerical diagnostic. Original central-difference gate is retained.
suppressPackageStartupMessages({library(antitrustBayes); library(rstan)})
repo <- "/mnt/d/projects/antitrustbayes"
recovery <- file.path(repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery")
role <- Sys.getenv("HELMERT_ROLE", "candidate")
artifact <- readRDS(file.path(recovery, "extra-native",
  paste0("extra_helmert_", role, ".rds")))
case <- "ces_second_score_auction"
point <- artifact$results[[case]]
stopifnot(identical(point$status, "ok"))
source(file.path(repo,
  "doc/performance_audit/recovery/native_review/extra_helmert_fixtures.R"))
dat <- do.call(make_extra_data, extra_cases[[case]])
model <- get("stanmodels", asNamespace("antitrustBayes"))[[point$kernel]]
stopifnot(identical(digest::digest(model@model_code, algo="sha256",
  serialize=FALSE), point$kernel_code_hash))
sf <- sampling(model, data=dat, chains=1, iter=1, warmup=0,
  algorithm="Fixed_param", init=0, refresh=0, seed=58003)
u <- point$u
ad <- grad_log_prob(sf, u, adjust_transform=TRUE)
lp <- function(z) log_prob(sf, z, adjust_transform=TRUE)
indices <- point$fd$index
# Frozen numerical-sensitivity grid. No HMC, prior or economic settings change.
steps <- c(1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
rows <- list()
for (i in indices) for (h in steps) {
  value <- function(offset) {
    z <- u; z[i] <- z[i]+offset
    lp(z)
  }
  minus <- value(-h); plus <- value(h)
  minus_half <- value(-h/2); plus_half <- value(h/2)
  central <- (plus-minus)/(2*h)
  half <- (plus_half-minus_half)/h
  # Richardson removes the central difference's leading h^2 error.
  fourth <- (4*half-central)/3
  rows[[length(rows)+1L]] <- data.frame(role=role, case=case,
    coordinate=i, step=h, minus=minus, plus=plus,
    minus_half=minus_half, plus_half=plus_half,
    central=central, richardson=fourth, analytic=ad[i],
    central_error=abs(central-ad[i]), richardson_error=abs(fourth-ad[i]))
}
write.csv(do.call(rbind, rows), file.path(repo,
  "doc/performance_audit/recovery/native_review",
  paste0("fd_cancellation_", role, ".csv")), row.names=FALSE)
