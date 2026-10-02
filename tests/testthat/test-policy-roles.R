test_that("rich-game roles follow post-merger ownership without reranking", {
  for (conduct in c("BF", "BS")) {
    market <- suppressWarnings(bayes_calibrated_market(
      conduct = conduct, J = 5L, k = 2L, seed = 9103L,
      dirichlet_concentration = 1,
      mode = "parameter", alpha = 12,
      costs = c(0.8, 1.1, 0.9, 1.2, 1.0)
    ))
    model <- market$truth$model
    ranking <- as.character(market$design$role_ranking)
    owner_post <- as.character(model@firmOwnerPre)
    owner_post[owner_post == ranking[[2L]]] <- ranking[[1L]]
    expected <- if (identical(conduct, "BS"))
      list(leadersPost = ranking[[1L]]) else
      list(corePost = ranking[[1L]])
    expect_identical(bayes_policy_roles(model, owner_post), expected)
    expect_identical(bayes_policy_roles(market, owner_post), expected)

    # Posterior structural models may prefix pre-owner IDs while a policy
    # manifest uses raw product IDs after the merger.
    prefixed <- model
    prefixed@firmOwnerPre <- paste0("market_1_", model@firmOwnerPre)
    if (identical(conduct, "BS"))
      prefixed@leadersPre <- paste0("market_1_", model@leadersPre) else
      prefixed@corePre <- paste0("market_1_", model@corePre)
    expect_identical(bayes_policy_roles(prefixed, owner_post), expected)

    bad <- model
    if (identical(conduct, "BS")) bad@leadersPre <- "unknown" else
      bad@corePre <- "unknown"
    expect_error(bayes_policy_roles(bad, owner_post),
                 "Pre-policy firm roles")
  }
  expect_error(bayes_policy_roles(new.env(), c("a", "b")),
               "native Stackelberg or exact core-fringe")
})
