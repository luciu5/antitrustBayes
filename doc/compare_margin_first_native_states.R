old <- read.csv("doc/collapsed_joint_conditional_native_state_values.csv",
                stringsAsFactors = FALSE)
new <- read.csv("doc/collapsed_joint_margin_first_native_state_values.csv",
                stringsAsFactors = FALSE)

keys <- c("kernel", "scenario", "coordinate", "field", "element")
comparison <- merge(old, new, by = keys, all = TRUE,
                    suffixes = c("_conditional", "_margin_first"))
comparison$abs_error <- abs(comparison$value_conditional -
                              comparison$value_margin_first)
comparison$tolerance <- 3e-8
comparison$passed <- is.finite(comparison$abs_error) &
  comparison$abs_error <= comparison$tolerance

write.csv(comparison,
          "doc/collapsed_joint_margin_first_native_gq_checks.csv",
          row.names = FALSE)

cat(sprintf("native generated quantities %d/%d; maximum error %.12g\n",
            sum(comparison$passed), nrow(comparison),
            max(comparison$abs_error)))
if (any(!comparison$passed)) {
  print(comparison[!comparison$passed, , drop = FALSE])
  stop("conditional and margin-first generated quantities differ")
}
