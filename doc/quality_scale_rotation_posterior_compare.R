#!/usr/bin/env Rscript

## Compare the economic posterior from the legacy scale coordinates with the
## exact total-scale/variance-fraction coordinates.  Both inputs are compact
## posterior summaries produced under seed 29091 by the frozen synthetic HMC
## benchmark.

legacy_path <- "doc/s0_joint_rotation_posterior_equivalence.csv"
rotated_path <- "doc/quality_scale_rotation_posteriors.csv"
output_path <- "doc/quality_scale_rotation_posterior_equivalence_checks.csv"

legacy <- read.csv(legacy_path, stringsAsFactors = FALSE)
legacy <- legacy[
  legacy$case == "market_joint_nc" & legacy$seed == 29091,
  c("variable", "mean", "sd", "mcse_mean")
]
rotated <- read.csv(rotated_path, stringsAsFactors = FALSE)
rotated <- rotated[
  rotated$case == "market_joint_nc" & rotated$seed == 29091,
  c("variable", "mean", "sd", "mcse_mean")
]

names(legacy)[-1L] <- paste0("legacy_", names(legacy)[-1L])
names(rotated)[-1L] <- paste0("total_fraction_", names(rotated)[-1L])
out <- merge(legacy, rotated, by = "variable", all = FALSE, sort = FALSE)
if (!nrow(out)) stop("No common economic posterior variables were found.")

out$difference <- out$total_fraction_mean - out$legacy_mean
out$combined_mcse <- sqrt(
  out$legacy_mcse_mean^2 + out$total_fraction_mcse_mean^2
)
out$abs_difference_over_combined_mcse <-
  abs(out$difference) / out$combined_mcse
out$pooled_sd <- sqrt((out$legacy_sd^2 + out$total_fraction_sd^2) / 2)
out$abs_difference_over_pooled_sd <- abs(out$difference) / out$pooled_sd
out$passed <- out$abs_difference_over_combined_mcse <= 3 |
  out$abs_difference_over_pooled_sd <= .10

write.csv(out, output_path, row.names = FALSE)
cat(sum(out$passed), "/", nrow(out), "economic posterior checks passed\n",
    sep = "")
if (!all(out$passed)) quit(status = 1L)
