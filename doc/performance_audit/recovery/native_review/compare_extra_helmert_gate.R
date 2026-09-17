root <- Sys.getenv("HELMERT_REVIEW_ROOT",
  "doc/performance_audit/recovery/native_review")
artifact_dir <- Sys.getenv("HELMERT_ARTIFACT_DIR",
  file.path(root, ".artifacts"))
candidate_file <- file.path(artifact_dir, "extra_helmert_candidate.rds")
old_file <- file.path(artifact_dir, "extra_helmert_old.rds")
out_file <- file.path(root, "extra_helmert_comparison.csv")
failure_file <- file.path(root, "extra_helmert_comparison_failures.csv")
if (!file.exists(candidate_file) || !file.exists(old_file))
  stop("both candidate and old RDS files are required")

candidate <- readRDS(candidate_file)
old <- readRDS(old_file)
tol <- 1e-10

vector_metrics <- function(a, b) {
  nonfinite_signature <- function(x) {
    vapply(x, function(z) {
      if (is.na(z)) return(if (is.nan(z)) "NaN" else "NA")
      if (is.infinite(z)) return(if (z > 0) "+Inf" else "-Inf")
      "finite"
    }, character(1))
  }
  if (length(a) != length(b))
    return(list(same_length = FALSE, same_names = FALSE, finite_pattern = FALSE,
                nonfinite_types = FALSE, all_finite_a = FALSE,
                all_finite_b = FALSE, max_abs = Inf, max_rel = Inf))
  same_names <- identical(names(a), names(b))
  finite_a <- is.finite(a); finite_b <- is.finite(b)
  both <- finite_a & finite_b
  d <- if (any(both)) abs(a[both] - b[both]) else numeric()
  scale <- if (any(both)) pmax(1, abs(a[both]), abs(b[both])) else numeric()
  list(same_length = TRUE, same_names = same_names,
       finite_pattern = identical(finite_a, finite_b),
       nonfinite_types = identical(nonfinite_signature(a), nonfinite_signature(b)),
       all_finite_a = all(finite_a), all_finite_b = all(finite_b),
       max_abs = if (length(d)) max(d) else 0,
       max_rel = if (length(d)) max(d / scale) else 0)
}

quality_metrics <- function(a, b) {
  va <- a$quality; vb <- b$quality
  ma <- vector_metrics(va, vb)
  list(metrics = ma,
       candidate_centering = if (length(va)) sum(va) else NA_real_,
       old_centering = if (length(vb)) sum(vb) else NA_real_,
       candidate_first = if (length(va)) va[1] else NA_real_,
       old_first = if (length(vb)) vb[1] else NA_real_,
       candidate_last = if (length(va)) va[length(va)] else NA_real_,
       old_last = if (length(vb)) vb[length(vb)] else NA_real_,
       candidate_centering_ok = length(va) == 0L ||
         (all(is.finite(va)) && abs(sum(va)) <= tol),
       old_centering_ok = length(vb) == 0L ||
         (all(is.finite(vb)) && abs(sum(vb)) <= tol))
}

kernel_hash <- function(x) {
  if (!is.null(x$kernel_code_hash)) x$kernel_code_hash else x$kernel_hash
}

one <- function(nm) {
  a <- candidate$results[[nm]]
  b <- old$results[[nm]]
  if (!identical(a$status, "ok") || !identical(b$status, "ok")) {
    return(data.frame(case = nm, status_candidate = a$status,
      status_old = b$status, pass = FALSE,
      failure = paste(c(if (!identical(a$status, "ok")) a$error else NULL,
                        if (!identical(b$status, "ok")) b$error else NULL),
                      collapse = " || "), stringsAsFactors = FALSE))
  }
  lp <- vector_metrics(a$lp, b$lp)
  lp_raw <- vector_metrics(a$lp_raw, b$lp_raw)
  gr <- vector_metrics(a$grad, b$grad)
  gr_raw <- vector_metrics(a$grad_raw, b$grad_raw)
  fd <- vector_metrics(a$fd$value, b$fd$value)
  tp <- vector_metrics(a$tp, b$tp)
  q <- vector_metrics(a$q, b$q)
  quality <- quality_metrics(a, b)
  schema_equal <- identical(a$schema, b$schema)
  raw_u <- vector_metrics(a$u, b$u)
  metrics <- list(lp = lp, lp_raw = lp_raw, grad = gr, grad_raw = gr_raw,
                  fd = fd, tp = tp, q = q, quality = quality$metrics,
                  raw_u = raw_u)
  max_abs <- vapply(metrics, function(x) x$max_abs, numeric(1))
  max_rel <- vapply(metrics, function(x) x$max_rel, numeric(1))
  numeric_ok <- all(is.finite(max_abs), is.finite(max_rel)) &&
    all(max_abs <= tol, max_rel <= tol)
  structure_ok <- schema_equal && raw_u$same_length && raw_u$same_names &&
    tp$same_length && tp$same_names && q$same_length && q$same_names
  required_finite <- c("lp", "lp_raw", "grad", "grad_raw", "fd", "tp", "raw_u")
  finite_ok <- all(vapply(metrics[required_finite], function(x)
    isTRUE(x$all_finite_a) && isTRUE(x$all_finite_b), logical(1))) &&
    isTRUE(q$finite_pattern) && isTRUE(q$nonfinite_types)
  quality_centering_ok <- isTRUE(quality$candidate_centering_ok) &&
    isTRUE(quality$old_centering_ok) &&
    length(a$quality) == a$N_group && length(b$quality) == b$N_group
  input_hash_equal <- identical(a$input_hash, b$input_hash)
  pass <- numeric_ok && structure_ok && finite_ok && quality_centering_ok &&
    input_hash_equal
  data.frame(
    case = nm, status_candidate = a$status, status_old = b$status, pass = pass,
    kernel = a$kernel, N = a$N, N_market = a$N_market,
    N_market_year = a$N_market_year, N_group = a$N_group,
    num_upars = a$schema$num_upars, schema_equal = schema_equal,
    raw_u_same = raw_u$same_length && raw_u$same_names,
    input_hash_equal = input_hash_equal,
    input_hash_candidate = a$input_hash, input_hash_old = b$input_hash,
    kernel_hash_candidate = kernel_hash(a), kernel_hash_old = kernel_hash(b),
    normalized_target_abs = lp$max_abs, normalized_target_rel = lp$max_rel,
    raw_target_abs = lp_raw$max_abs, raw_target_rel = lp_raw$max_rel,
    normalized_gradient_abs = gr$max_abs, normalized_gradient_rel = gr$max_rel,
    raw_gradient_abs = gr_raw$max_abs, raw_gradient_rel = gr_raw$max_rel,
    fd_abs = fd$max_abs, fd_rel = fd$max_rel,
    transformed_abs = tp$max_abs, transformed_rel = tp$max_rel,
    gq_abs = q$max_abs, gq_rel = q$max_rel,
    finite_targets = lp$finite_pattern && lp_raw$finite_pattern,
    finite_gradients = gr$finite_pattern && gr_raw$finite_pattern,
    finite_fd = fd$finite_pattern, finite_transformed = tp$finite_pattern,
    finite_gq = q$finite_pattern, gq_nonfinite_types = q$nonfinite_types,
    quality_centering_ok = quality_centering_ok,
    quality_centering_bound = tol, quality_len = length(a$quality),
    candidate_centering = quality$candidate_centering,
    old_centering = quality$old_centering,
    candidate_first_quality = quality$candidate_first,
    old_first_quality = quality$old_first,
    candidate_last_quality = quality$candidate_last,
    old_last_quality = quality$old_last,
    elapsed_candidate = a$elapsed, elapsed_old = b$elapsed,
    failure = if (pass) "" else paste(c(
      if (!schema_equal) "schema", if (!structure_ok) "raw_schema",
      if (!finite_ok) "finite_pattern", if (!numeric_ok) "numeric_tolerance",
      if (!quality_centering_ok) "quality_centering",
      if (!input_hash_equal) "input_hash"
    ), collapse = ";"), stringsAsFactors = FALSE)
}

expected_cases <- c(
  "ces_g1000_repeated_panel", "ces_sampled_margin_only_panel",
  "ces_second_score_auction", "logit_bargaining",
  "logit_g1000_single_product", "logit_sampled_share_only_panel",
  "logit_smooth_cutoff_hmt")
validate_manifest <- function(x, label) {
  if (!is.character(x$cases) || anyDuplicated(x$cases) ||
      !identical(sort(x$cases), expected_cases))
    stop(label, " case manifest must contain exactly the seven expected cases")
  result_names <- names(x$results)
  if (!is.list(x$results) || is.null(result_names) || anyDuplicated(result_names) ||
      !identical(sort(result_names), expected_cases))
    stop(label, " result manifest must contain exactly the seven expected cases")
}
validate_manifest(candidate, "candidate")
validate_manifest(old, "old")
case_names <- expected_cases
rows <- lapply(case_names, one)
columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
rows <- lapply(rows, function(row) {
  for (column in setdiff(columns, names(row))) row[[column]] <- NA
  row[columns]
})
out <- do.call(rbind, rows)
write.csv(out, out_file, row.names = FALSE)
fail <- out[!out$pass, , drop = FALSE]
if (!nrow(fail)) fail <- out[FALSE, , drop = FALSE]
write.csv(fail, failure_file, row.names = FALSE)
pass_count <- sum(out$pass)
cat("extra Helmert native review:", pass_count, "/", nrow(out), "cases passed\n")
if (pass_count != nrow(out)) quit(status = 1L)
