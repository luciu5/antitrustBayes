suppressPackageStartupMessages({
  library(antitrustBayes)
  library(rstan)
})

ns <- asNamespace("antitrustBayes")
mods <- get("stanmodels", ns)

make_data <- function(demand, structure = "simultaneous", quality = "sampled",
                      G = 3L, supply_model = 1L) {
  N <- 3L
  args <- list(
    prices = seq(1, 2, length.out = N), quantities = rep(1, N),
    margins = rep(.2, N), market = rep("M", N), year = rep("Y", N),
    group = paste0("G", ((seq_len(N) - 1L) %% G) + 1L),
    supply_model = supply_model,
    use_cutoff = FALSE, use_rho = FALSE, use_market_resid_scale = FALSE,
    use_hmt = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_s0_market_pool = FALSE,
    use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
    market_s0_parameterization = "legacy",
    priors = get("resolve_prior_manifest", ns)("baseline", list()),
    grainsize = 1L, margin_observed = NULL, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = rep("M", N), estimate_df = FALSE,
    resid_df_fixed = 5, structure = structure, cutoff_scale = NULL,
    partition_prior = "cutoff", k_support = NULL, k_prior = NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian", quality_parameterization = quality
  )
  if (structure != "simultaneous") {
    args$partition_prior <- "ordered_k"
    args$k_support <- if (structure == "fringe") 0:3 else 1:3
    args$k_prior <- rep(1 / length(args$k_support), length(args$k_support))
  }
  fun <- if (demand == "logit") get("prepare_stan_data", ns)
    else get("prepare_ces_stan_data", ns)
  if (demand == "ces") args$market_side <- "output"
  dat <- do.call(fun, args)
  dat$structure <- switch(structure, simultaneous = 0L, fringe = 1L,
                          stackelberg = 2L)
  dat
}

one <- function(kernel, demand, structure = "simultaneous", G = 3L,
                beta = c(1, 2), supply_model = 1L) {
  dat <- make_data(demand, structure, G = G, supply_model = supply_model)
  model <- mods[[kernel]]
  sf0 <- suppressWarnings(sampling(model, data = dat, chains = 1L,
    iter = 1L, warmup = 0L, algorithm = "Fixed_param", init = 0,
    refresh = 0L, seed = 48003L))
  u0 <- rep(0, get_num_upars(sf0))
  pars <- constrain_pars(sf0, u0)
  pars$b_group_centered <- beta
  u <- unconstrain_pars(sf0, pars)
  lp <- log_prob(sf0, u, adjust_transform = TRUE, gradient = FALSE)
  gr <- grad_log_prob(sf0, u, adjust_transform = TRUE)
  gate_id <- paste(kernel, structure, supply_model, G, sep="_")
  write.csv(data.frame(index=seq_along(gr), gradient=gr), file.path(Sys.getenv("GATE_DETAIL"),paste0(gate_id,"_gradient.csv")),row.names=FALSE)
  h <- 1e-5
  fd <- vapply(seq_along(u), function(i) {
    up <- u
    um <- u
    up[i] <- up[i] + h
    um[i] <- um[i] - h
    (log_prob(sf0, up, adjust_transform = TRUE, gradient = FALSE) -
       log_prob(sf0, um, adjust_transform = TRUE, gradient = FALSE)) / (2 * h)
  }, numeric(1))
  sf <- suppressWarnings(sampling(model, data = dat, chains = 1L,
    iter = 1L, warmup = 0L, algorithm = "Fixed_param",
    init = function(...) pars, refresh = 0L, seed = 48005L))
  q <- as.numeric(extract(sf, permuted = TRUE)$b_quality_all)
  saveRDS(extract(sf,permuted=FALSE,inc_warmup=FALSE), file.path(Sys.getenv("GATE_DETAIL"),paste0(gate_id,"_output.rds")))
  expected <- if (G == 1L) 0 else c(1 / sqrt(2) + 2 / sqrt(6),
    2 / sqrt(6) - 1 / sqrt(2), -4 / sqrt(6))
  if (max(abs(q - expected)) > 1e-13)
    stop("nonzero Helmert recurrence mismatch for ", kernel)
  data.frame(kernel = kernel, demand = demand, structure = structure,
    lp = as.numeric(lp), grad_l2 = sqrt(sum(gr^2)),
    fd_max_abs = max(abs(gr - fd)),
    fd_max_rel = max(abs(gr - fd) / (1 + abs(fd))),
    q1 = q[1], q2 = q[2], q3 = q[3], q_expected_max_abs = max(abs(q - expected)),
    stringsAsFactors = FALSE)
}

out <- do.call(rbind, list(
  one("logit_supply", "logit"),
  one("logit_stackelberg", "logit", "fringe"),
  one("logit_stackelberg", "logit", "fringe", supply_model = 3L),
  one("logit_stackelberg", "logit", "stackelberg"),
  one("logit_stackelberg", "logit", "stackelberg", supply_model = 3L),
  one("ces_supply", "ces"),
  one("ces_stackelberg", "ces", "fringe"),
  one("ces_stackelberg", "ces", "fringe", supply_model = 3L),
  one("ces_stackelberg", "ces", "stackelberg"),
  one("ces_stackelberg", "ces", "stackelberg", supply_model = 3L),
  one("logit_supply", "logit", G = 1L, beta = numeric()),
  one("ces_supply", "ces", G = 1L, beta = numeric())
))
print(out, row.names = FALSE, digits = 17)

write.csv(out, Sys.getenv("GATE_CSV"), row.names=FALSE)
meta <- data.frame(package_path=find.package("antitrustBayes"), so_md5=unname(tools::md5sum(system.file("libs", "antitrustBayes.so", package="antitrustBayes"))), code_sha256=paste(vapply(mods, function(x) digest::digest(x@model_code, algo="sha256", serialize=FALSE), character(1)), collapse=";"))
write.csv(meta,paste0(Sys.getenv("GATE_CSV"),".manifest.csv"),row.names=FALSE)
