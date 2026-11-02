## Prepare fixed, collapsed-quality stress inputs.  This is a benchmark
## fixture generator only: it does not fit a posterior or change package
## source.  Each N uses G=min(N,1000) persistent groups and repeated years.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript ... N output.rds")
N <- as.integer(args[[1L]]); out <- args[[2L]]
if (!N %in% c(1000L, 10000L)) stop("N must be 1000 or 10000")
G <- min(N, 1000L)
T <- N %/% G
if (G * T != N) stop("N must be divisible by G")
e <- new.env(parent = globalenv())
sys.source("doc/collapsed_gaussian_benchmark.R", envir = e)
panel <- e$cga_make_dgp(
  seed = 20260916L, demand_model = "logit", supply_model = "bertrand",
  structure = "simultaneous", share_scenario = "noisy",
  n_markets = 1L, n_periods = T, n_products = G
)
dat <- panel$data
stopifnot(nrow(dat) == N, length(unique(dat$group)) == G,
          all(is.finite(dat$prices)), all(is.finite(dat$quantities)),
          all(is.finite(dat$margins)), all(dat$margins > 0),
          all(dat$margins < 1))
saveRDS(list(
  N = N, G = G, T = T, seed = 20260916L,
  data = dat, outside_price = as.numeric(panel$parameters$outside_price),
  description = "Fixed_param collapsed-quality stress fixture; not posterior inference."
), out)
cat("prepared N=", N, " G=", G, " T=", T, " rows=", nrow(dat), "\n", sep = "")
