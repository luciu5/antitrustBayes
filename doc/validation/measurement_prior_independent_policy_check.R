# Independent input-Logit demand and Broyden equilibrium audit. No fitting.
# Detects solver/index errors that a replay of the original evaluator would miss.
args <- commandArgs(TRUE)
output <- normalizePath(if (length(args)) args[[1L]] else
  ".scratch/measurement_prior_diagnostic_20261109", mustWork = TRUE)
.libPaths(c("/mnt/d/Projects/mergerBayes/.scratch/paper-library", .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
inputs <- readRDS(file.path(output, "inputs.rds"))
extract <- get(".bayes_backend_extract", asNamespace("antitrustBayes"))
results <- list()
for (cell in c("snapshot_estimated", "history_estimated")) {
  policy <- readRDS(file.path(output, cell, "policy_draws.rds"))
  saved <- readRDS(file.path(output, cell, "fit.rds"))
  input <- inputs[[sub("_.*", "", cell)]]
  target <- which(input$rows$year == "1960"); rows <- input$rows[target, ]
  owner <- as.character(rows$group); post <- owner
  post[post == rows$target[1L]] <- rows$acquirer[1L]
  parties <- owner %in% c(rows$target[1L], rows$acquirer[1L])
  own_pre <- outer(owner, owner, "==") * 1
  own_post <- outer(post, post, "==") * 1
  a <- extract(saved$native_fit, variables = "a_market")$a_market
  pre_foc <- post_foc <- index_error_pp <- 0
  for (j in seq_len(10L)) {
    d <- policy$diagnostics[[j]]
    alpha <- as.numeric(a[policy$indices[j], 1L]) / input$data$price_sd
    share <- function(price) {
      w <- d$reference_quantity_shares * exp(alpha * (price - d$reference_price))
      w / (1 - sum(d$reference_quantity_shares) + sum(w))
    }
    foc <- function(price, ownership) d$mc_pre - price -
      1 / (alpha * (1 - as.numeric(ownership %*% share(price))))
    pre_foc <- max(pre_foc, max(abs(foc(d$latent_pre_price, own_pre))))
    solved <- nleqslv::nleqslv(d$latent_pre_price, function(price) foc(price, own_post),
      method = "Broyden", control = list(xtol = 1e-12, ftol = 1e-12, maxit = 1000L))
    post_foc <- max(post_foc, max(abs(solved$fvec)))
    s <- share(solved$x)
    market <- sum(solved$x * s) / sum(d$latent_pre_price * s) - 1
    party <- sum(solved$x[parties] * s[parties]) /
      sum(d$latent_pre_price[parties] * s[parties]) - 1
    index_error_pp <- max(index_error_pp,
      100 * abs(market - policy$effects$Market_Paasche[j]),
      100 * abs(party - policy$effects$Party_Paasche[j]))
  }
  stopifnot(pre_foc < 2e-8, post_foc < 1e-8, index_error_pp < 1e-5)
  results[[cell]] <- data.frame(Cell = cell, Draws = 10L, Max_Pre_FOC = pre_foc,
    Max_Post_FOC = post_foc, Max_Paasche_Error_PP = index_error_pp,
    Posterior_Rows = paste(policy$indices[seq_len(10L)], collapse = ","))
}
results <- do.call(rbind, results)
write.csv(results, file.path(output, "independent_policy_audit.csv"), row.names = FALSE)
print(results)
