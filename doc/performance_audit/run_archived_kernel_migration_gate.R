# Validate the composed original-to-Helmert path on actual archived posteriors.
library(antitrustBayes)
ns <- asNamespace("antitrustBayes")
e <- new.env(parent = globalenv())
for (n in ls(ns, all.names = TRUE)) assign(n, get(n, ns), e)
sys.source("R/bma.R", e)
archive <- Sys.getenv("MIGRATION_ARCHIVE_DIR")
out <- Sys.getenv("MIGRATION_ARCHIVE_OUTPUT")
stopifnot(nzchar(archive), nzchar(out))
files <- c("logit_bertrand_noisy_271828__D.rds", "ces_bertrand_noisy_271828__D.rds",
           "logit_BF_noisy_271828__D.rds", "logit_BS_noisy_271828__D.rds",
           "logit_cournot_noisy_271828__D.rds")
rows <- lapply(files, function(file) {
  cat(file, "\n"); flush.console()
  fit <- readRDS(file.path(archive, file))
  before <- serialize(fit, NULL); set.seed(9931L); before_rng <- .Random.seed
  ctx <- e$.bayes_bridge_restore_context(fit)
  contract <- e$.bayes_bridge_code_contract(fit@stanfit)
  data <- e$.bayes_bridge_private_data(fit, contract)
  stopifnot(identical(contract$migration_geometry, "legacy"),
            identical(data$market_s0_parameterization, 1L),
            identical(.Random.seed, before_rng), identical(serialize(fit, NULL), before))
  data.frame(file = file, kernel = ctx$kernel, checked_rows = nrow(ctx$target_check),
    max_normalized_target_difference = max(ctx$target_check$abs_difference),
    restored_flag = data$market_s0_parameterization,
    added_zero_width_fields = paste(ctx$zero_width_fields_added, collapse = ";"),
    rng_preserved = TRUE, fit_immutable = TRUE, raw_schema = "PASS",
    scope = "archived normalized target restoration; not a bridge evidence estimate")
})
write.csv(do.call(rbind, rows), out, row.names = FALSE)
