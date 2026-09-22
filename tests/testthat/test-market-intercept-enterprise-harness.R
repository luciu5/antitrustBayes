.mic_enterprise_env <- local({
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "inst", "benchmarks",
                                "market_intercept_enterprise",
                                "common.R"), envir = env)
  env
})

test_that("enterprise manifest is deterministic, unique, and staged", {
  one <- .mic_enterprise_env$mic_make_manifest()
  two <- .mic_enterprise_env$mic_make_manifest()
  .mic_enterprise_env$mic_validate_manifest(one)
  expect_identical(one, two)
  expect_identical(nrow(one), 124L)
  expect_identical(anyDuplicated(one$run_key), 0L)
  expect_equal(as.vector(table(one$tier)[c("gate", "replicate", "comparison", "sc",
                                           "sc_comparison", "stress", "scale")]),
               c(32L, 32L, 16L, 16L, 4L, 8L, 16L))
  expect_true(all(c("B", "C", "MC", "A", "BF", "CF", "BS", "CS") %in%
                    unique(one$model)))
  expect_true(all(c("logit", "ces") %in% unique(one$demand)))
  expect_true(all(one$benchmark_kind[one$tier == "scale"] == "timing"))
  expect_true(all(one$chains[one$tier == "scale"] == 1L))
  expect_true(all(one$true_k[one$model %in% c("BF", "CF", "BS", "CS") &
                               one$source == "synthetic"] == 2L))
  frozen <- one
  frozen$manifest_hash <- .mic_enterprise_env$mic_manifest_hash(frozen)
  expect_identical(.mic_enterprise_env$mic_manifest_hash(frozen),
                   unique(frozen$manifest_hash))
  frozen$model[[1L]] <- "C"
  expect_false(identical(.mic_enterprise_env$mic_manifest_hash(frozen),
                         unique(frozen$manifest_hash)))
})

test_that("all synthetic enterprise cells satisfy the no-sampling contract", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "doc",
                                             "market_quality_benchmark.R")),
              "source-repository DGP helpers are excluded from package builds")
  plan <- .mic_enterprise_env$mic_make_manifest()
  plan <- plan[plan$source == "synthetic", , drop = FALSE]
  checked <- lapply(seq_len(nrow(plan)), function(i)
    .mic_enterprise_env$mic_smoke_task(plan[i, , drop = FALSE]))
  checked <- do.call(rbind, checked)
  expect_identical(nrow(checked), nrow(plan))
  expect_true(all(checked$status == "smoke_passed"))
  expect_identical(as.integer(checked$N_market), as.integer(plan$M))
  expect_true(all(checked$observed_margins > 0L))
  expect_true(all(checked$observed_margins <= checked$N))
})

test_that("ordered-k support counts firms and parameterizations remain paired", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "doc",
                                             "market_quality_benchmark.R")),
              "source-repository DGP helpers are excluded from package builds")
  plan <- .mic_enterprise_env$mic_make_manifest()
  rich <- plan[plan$source == "synthetic" & plan$model == "BS" &
                 plan$demand == "logit" & plan$T == 4L, , drop = FALSE]
  collapsed <- rich[rich$intercept == "collapsed", , drop = FALSE][1L, ]
  sampled <- rich[rich$intercept == "sampled", , drop = FALSE][1L, ]
  fixture <- .mic_enterprise_env$mic_load_fixture(collapsed)
  cargs <- .mic_enterprise_env$mic_fit_args(collapsed, fixture)
  sargs <- .mic_enterprise_env$mic_fit_args(sampled, fixture)
  expect_identical(cargs$k_support, 1:5)
  expect_equal(sum(cargs$k_prior), 1)
  expect_identical(cargs$market_intercept_parameterization, "collapsed")
  expect_identical(cargs$market_s0_parameterization, "legacy")
  expect_identical(sargs$market_intercept_parameterization, "sampled")
  expect_identical(sargs$market_s0_parameterization, "joint")
  expect_identical(cargs$use_rho, TRUE)
})

test_that("stress fixtures implement the declared statistical regimes", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "doc",
                                             "market_quality_benchmark.R")),
              "source-repository DGP helpers are excluded from package builds")
  plan <- .mic_enterprise_env$mic_make_manifest()
  high <- plan[plan$tier == "stress" & plan$scenario == "high_rho", , drop = FALSE][1L, ]
  low <- plan[plan$tier == "stress" & plan$scenario == "low_share_noise", , drop = FALSE][1L, ]
  partial <- plan[plan$tier == "stress" & plan$scenario == "partial_margins", , drop = FALSE][1L, ]
  high_fx <- .mic_enterprise_env$mic_load_fixture(high)
  low_fx <- .mic_enterprise_env$mic_load_fixture(low)
  base_fx <- .mic_enterprise_env$mic_make_synthetic(low$model, low$demand,
    low$M, low$T, low$F, low$dgp_seed, "compositional")
  partial_fx <- .mic_enterprise_env$mic_load_fixture(partial)
  expect_equal(high_fx$truth$rho, .9)
  expect_lt(low_fx$truth$share_residual_sd,
            stats::sd(log(base_fx$data$quantities /
              ave(base_fx$data$quantities,
                  interaction(base_fx$data$market, base_fx$data$year), FUN = sum)) -
              base_fx$truth$latent_log_share))
  first <- unique(partial_fx$data$market)[[1L]]
  expect_false(any(partial_fx$data$margin_observed[partial_fx$data$market == first]))
  expect_true(any(partial_fx$data$margin_observed[partial_fx$data$market != first]))
})

test_that("atomic DGP conduct schedules match the public package equations", {
  theory <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "R", "priors.R"), envir = theory)
  sys.source(testthat::test_path("..", "..", "R", "theory.R"), envir = theory)
  q <- c(.03, .11, .24)
  s0 <- .27
  for (model in c("B", "C", "MC", "A")) {
    supply <- .mic_enterprise_env$mic_model_contract(model)$supply
    expect_equal(
      .mic_enterprise_env$mic_conduct_schedule(q, model, "logit", s0),
      theory$conduct_schedule(q, supply, s0 = s0, demand_model = "logit"),
      tolerance = 1e-14
    )
    expect_equal(
      .mic_enterprise_env$mic_conduct_schedule(q, model, "ces", s0,
                                               gamma = 3.5),
      theory$conduct_schedule(q, supply, s0 = s0,
                              demand_model = "ces", gamma = 3.5),
      tolerance = 1e-14
    )
  }
})

test_that("task paths isolate concurrent workers", {
  root <- tempfile("mic-enterprise-test-")
  paths <- vapply(1:20, function(i)
    .mic_enterprise_env$mic_task_result_path(root, i), character(1))
  expect_identical(anyDuplicated(paths), 0L)
  expect_true(all(mapply(grepl, sprintf("task-%04d", 1:20), paths,
                         MoreArgs = list(fixed = TRUE))))
})
