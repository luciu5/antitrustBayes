# Rebuild the compact table from retained executed measurements.
# Memory ratios explicitly distinguish process memory from copied array sizes.
out_dir <- "doc/performance_audit"
rows <- list()
add <- function(change, benchmark, old, new, memory, difference, diagnostics,
                recommendation, memory_basis, evidence_basis, seed = NA_integer_) {
  rows[[length(rows) + 1L]] <<- data.frame(change = change, benchmark = benchmark,
    seed = seed, old_wall_time = old, new_wall_time = new,
    speed_ratio = if (new > 0) old / new else NA_real_, memory_ratio = memory,
    numerical_difference = difference, diagnostics_failures = diagnostics,
    recommendation = recommendation, memory_basis = memory_basis,
    evidence_basis = evidence_basis)
}
c <- read.csv(file.path(out_dir, "chain_results.csv"))
for (i in which(c$cores == 1L)) {
  j <- which(c$cores == 4L & c$fixture == c$fixture[i] & c$seed == c$seed[i])
  stopifnot(length(j) == 1L, c$posterior_draw_hash[i] == c$posterior_draw_hash[j])
  add("parallel chains, cores 1 to 4", c$fixture[i], c$wall_seconds[i],
    c$wall_seconds[j], c$process_tree_peak_pss_bytes[j] / c$process_tree_peak_pss_bytes[i],
    0, sprintf("divergences %d/%d; depth hits %d/%d", c$divergences[i],
      c$divergences[j], c$treedepth_hits[i], c$treedepth_hits[j]),
    "option-aware default; explicit 4 processes when memory permits",
    "process-tree peak PSS, 0.25-second polling", "paired proper posterior HMC", c$seed[i])
}
sel <- file.path(out_dir, "selective_extraction")
for (fixture in c("sim-logit", "sim-ces", "exact-BF", "exact-BS")) {
  a <- read.csv(file.path(sel, paste0("benchmark-policies-full-", fixture, ".csv")))
  b <- read.csv(file.path(sel, paste0("benchmark-policies-selected-", fixture, ".csv")))
  stopifnot(identical(a[c("fixture", "policy", "repeat.")],
                      b[c("fixture", "policy", "repeat.")]),
            all(a$valid == 200L), all(b$valid == 200L))
  for (policy in unique(a$policy)) {
    i <- which(a$policy == policy)
    add("selective extraction, total counterfactual", paste(fixture, policy, "200 draws x 3 repeats"),
      median(a$policy_seconds[i]), median(b$policy_seconds[i]),
      median(b$selected_or_full_bytes[i]) / median(a$selected_or_full_bytes[i]),
      0, "all 200 draws valid; prices/shares/welfare/identities/failures exact",
      "retain for copied-array reduction; whole-policy gains inconsistent",
      "copied extracted arrays; not process peak memory", "paired archived posterior simulation")
  }
}
for (n in c(1000L, 10000L)) {
  a <- read.csv(file.path(sel, paste0("stress-full-N", n, ".csv")))
  b <- read.csv(file.path(sel, paste0("stress-selected-N", n, ".csv")))
  stopifnot(identical(a$canonical_hashes, b$canonical_hashes), all(a$finite), all(b$finite))
  add("selective extraction, extraction only", paste0("N=", n, "; 200 Fixed_param rows"),
    a$retained_mode_seconds[1L], b$retained_mode_seconds[1L],
    b$retained_bytes[1L] / a$retained_bytes[1L], 0,
    "required economic hashes exact; finite; no HMC inference",
    "retain selective extraction; these are serialization measurements",
    "copied extracted arrays; not process peak memory", "Fixed_param serialization stress")
}
b <- read.csv(file.path(out_dir, "baseline_consumer_summary.csv"))
for (i in seq_len(nrow(b))) {
  add("public observed-baseline construction", paste(b$fixture[i], b$policy[i], "200 draws x 3 repeats"),
    b$old_wall_seconds[i], b$new_wall_seconds[i], b$process_peak_memory_ratio[i],
    b$max_numerical_difference[i],
    "all 200 valid; baseline calls 200 to 0; policy iterations/evaluations identical",
    "architecture candidate; no general wall-time gain; cross-package integration held",
    "whole benchmark process peak RSS, not individual policy", "paired archived posterior simulation")
}
h <- read.csv(file.path(out_dir, "native_profile/helmert_matched_large_comparison.csv"))
for (phase in c("setup", "target", "gradient")) {
  add("linear-time sampled-quality Helmert candidate", paste0("N=G=10000; ", phase),
    h[[paste0("baseline_", phase, "_seconds")]][1L],
    h[[paste0("optimized_", phase, "_seconds")]][1L],
    h$optimized_peak_rss_kib[1L] / h$baseline_peak_rss_kib[1L], h$target_abs_diff[1L],
    "finite; same target and gradient norm at fixed point; not HMC",
    "hold integration for independent panel/rho and saved-fit compatibility gates",
    "whole probe process peak RSS, not isolated phase", "matched-compiler fixed-point native probe")
}
result <- do.call(rbind, rows)
write.csv(result, file.path(out_dir, "summary.csv"), row.names = FALSE, na = "")
cat("summary rows:", nrow(result), "\n")
