#!/usr/bin/env Rscript

## Compare economic posterior summaries before and after replacing the exact
## conditional factorization of the collapsed joint Gaussian. The target is
## identical; differences should be Monte Carlo error.

conditional <- read.csv(
  "doc/quality_scale_rotation_posteriors.csv", stringsAsFactors = FALSE
)
margin_first <- read.csv(
  "doc/margin_first_posteriors.csv", stringsAsFactors = FALSE
)
select <- function(x) x[
  x$case == "market_joint_nc" & x$seed == 29091,
  c("variable", "mean", "sd", "mcse_mean")
]
conditional <- select(conditional)
margin_first <- select(margin_first)
names(conditional)[-1L] <- paste0("conditional_", names(conditional)[-1L])
names(margin_first)[-1L] <- paste0("margin_first_", names(margin_first)[-1L])
out <- merge(conditional, margin_first, by = "variable", all = FALSE,
             sort = FALSE)
if (!nrow(out)) stop("No common economic posterior variables were found.")

out$difference <- out$margin_first_mean - out$conditional_mean
out$combined_mcse <- sqrt(
  out$conditional_mcse_mean^2 + out$margin_first_mcse_mean^2
)
out$abs_difference_over_combined_mcse <-
  abs(out$difference) / out$combined_mcse
out$pooled_sd <- sqrt((out$conditional_sd^2 + out$margin_first_sd^2) / 2)
out$abs_difference_over_pooled_sd <- abs(out$difference) / out$pooled_sd
out$passed <- out$abs_difference_over_combined_mcse <= 3 |
  out$abs_difference_over_pooled_sd <= .10

write.csv(out, "doc/margin_first_posterior_equivalence_checks.csv",
          row.names = FALSE)
cat(sprintf("%d/%d economic posterior checks passed\n",
            sum(out$passed), nrow(out)))
if (!all(out$passed)) quit(status = 1L)
