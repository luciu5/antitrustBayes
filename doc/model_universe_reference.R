#!/usr/bin/env Rscript
# Independent finite-state evidence/conditioning calculation. No fitted model,
# bridge estimate, or generated-probability shortcut is used in this reference.
model_universe_reference <- function() {
  theta_prior <- c(.5, .5)
  k_prior <- c(.5, .5)
  # Rows are two shared continuous states; columns are k=0 and k=1.
  margin_likelihood <- matrix(c(1, 4, 4, 1), nrow = 2L, byrow = TRUE)
  state_likelihood <- as.numeric(margin_likelihood %*% k_prior)
  z_full <- sum(theta_prior * state_likelihood^2)
  z_restricted <- sum(theta_prior * margin_likelihood[, 2L]^2)
  theta_full <- theta_prior * state_likelihood^2 / z_full
  conditional_k1 <- margin_likelihood[, 2L] * k_prior[2L] /
    state_likelihood
  prior_event <- k_prior[2L]^2
  posterior_event <- sum(theta_full * conditional_k1^2)
  z_identity <- z_full * posterior_event / prior_event
  theta_restricted <- theta_full * conditional_k1^2 / posterior_event
  wrong_product_of_means <- sum(theta_full * conditional_k1)^2

  stopifnot(abs(z_full - 6.25) < 1e-14,
            abs(z_restricted - 8.5) < 1e-14,
            abs(z_identity - z_restricted) < 1e-14,
            abs(posterior_event - .34) < 1e-14,
            abs(wrong_product_of_means - .25) < 1e-14,
            max(abs(theta_restricted - c(16, 1) / 17)) < 1e-14)
  data.frame(
    quantity = c("full_evidence", "restricted_evidence_direct",
                 "restricted_evidence_identity", "prior_joint_event",
                 "posterior_joint_event", "incorrect_product_of_means",
                 "restricted_theta_state_1", "restricted_theta_state_2"),
    value = c(z_full, z_restricted, z_identity, prior_event, posterior_event,
              wrong_product_of_means, theta_restricted),
    status = c(rep("verified", 5L), "deliberately_incorrect_control",
               rep("verified", 2L))
  )
}

if (sys.nframe() == 0L) {
  result <- model_universe_reference()
  print(result, row.names = FALSE)
  write.csv(result, "doc/model_universe_reference_validation.csv", row.names = FALSE)
}
