# Preserve the original failed central-difference gate. Resolve only the
# diagnosed subtraction issue using a separate, more accurate derivative test.
root <- "/mnt/d/projects/antitrustbayes/doc/performance_audit/recovery/native_review"
gate <- read.csv(file.path(root, "extra_helmert_comparison.csv"))
tol <- 1e-10
stopifnot(nrow(gate) == 7L, !anyDuplicated(gate$case),
  all(gate$status_candidate == "ok"), all(gate$status_old == "ok"))
required_flags <- c("schema_equal", "raw_u_same", "input_hash_equal",
  "finite_targets", "finite_gradients", "finite_fd", "finite_transformed",
  "finite_gq", "gq_nonfinite_types", "quality_centering_ok")
stopifnot(all(vapply(gate[required_flags], function(x) all(x), logical(1))),
  all(gate$quality_len == gate$N_group))
required_deltas <- c("normalized_target_abs", "raw_target_abs",
  "normalized_gradient_abs", "raw_gradient_abs", "transformed_abs", "gq_abs")
stopifnot(all(vapply(gate[required_deltas], function(x)
  all(is.finite(x) & x <= tol), logical(1))))
failed <- which(!gate$pass)
stopifnot(length(failed) == 1L,
  gate$case[failed] == "ces_second_score_auction",
  gate$failure[failed] == "numeric_tolerance", gate$fd_abs[failed] > tol)
candidate <- read.csv(file.path(root, "fd_cancellation_candidate.csv"))
old <- read.csv(file.path(root, "fd_cancellation_old.csv"))
candidate <- candidate[candidate$step == 1e-3, ]
old <- old[old$step == 1e-3, ]
stopifnot(identical(candidate$coordinate, old$coordinate),
  identical(candidate$coordinate, c(1L, 10L, 19L)),
  all(candidate$case == gate$case[failed]),
  all(old$case == gate$case[failed]),
  all(is.finite(candidate$richardson)), all(is.finite(old$richardson)))
error <- max(candidate$richardson_error, old$richardson_error)
delta <- max(abs(candidate$richardson-old$richardson))
stopifnot(error <= tol, delta <= tol)
out <- data.frame(case=gate$case, central_v1_pass=gate$pass,
  resolved_pass=TRUE, fd_validation=ifelse(gate$pass, "central 1e-5",
    "fourth-order Richardson 1e-3; both analytic errors <=1e-10"),
  normalized_target_abs=gate$normalized_target_abs,
  full_gradient_abs=gate$normalized_gradient_abs,
  tp_gq_abs=pmax(gate$transformed_abs, gate$gq_abs),
  central_fd_abs=gate$fd_abs,
  refined_fd_abs=ifelse(gate$pass, NA_real_, delta),
  refined_analytic_error_max=ifelse(gate$pass, NA_real_, error),
  quality_len=gate$quality_len,
  quality_centering_ok=gate$quality_centering_ok)
write.csv(out, file.path(root, "extra_helmert_resolved_gate.csv"), row.names=FALSE)
cat("seven resolved cases; original six-of-seven central gate retained\n")
