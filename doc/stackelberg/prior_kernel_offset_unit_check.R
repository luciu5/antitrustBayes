.evidence_test_fit <- function(structure="simultaneous", raw=-20) {
  dat <- c(list(N=2L,shareIn=c(.4,.6),clr_shareIn=log(c(.4,.6))-mean(log(c(.4,.6))),
    marginNormInv=c(4,5),log_marginNormInv=log(c(4,5)),margin_obs=c(1L,1L),
    market_year_idx=c(1L,1L),cluster=c(1L,1L),likelihood_mode=1L,
    structure=if(structure=="stackelberg")2L else if(structure=="fringe")1L else 0L,
    use_cutoff=as.integer(structure!="simultaneous"),use_hmt=0L,
    use_s0_market_pool=0L,use_s0_market_year_pool=0L,k_scope_code=0L,
    mky_to_partition_unit=1L),bayes_prior_manifest("baseline"))
  attr(dat,"model_spec_version") <- 6L
  methods::new("BayesLogit",stanfit=list(raw=raw),stan_data=dat,
    data=data.frame(prices=c(1,2),quantities=c(4,6),margins=c(.25,.2),
      market="m",year="1",group=c("a","b"),loo_cluster="m",
      margin_observed=TRUE,outside_price=1),
    supply_model=1L,demand_model="logit",market_side="output",price_sd=1,
    call=call("bayes_est"),metadata=list(structure=structure,
      prior_manifest="baseline",resolved_priors=bayes_prior_manifest("baseline"),
      model_prior=NA_real_))
}
.evidence_test_dims <- function(structure="simultaneous") {
  list(mu_log_a_free=1L,b_group_raw=2L,sigma_b_strat=integer(),s0_logit=1L,
    log_sigma_share=integer(),log_sigma_margin=integer(),
    cutoff_share=as.integer(structure!="simultaneous"),
    mu_b_fringe_vec=as.integer(structure=="fringe"))
}

test_that("proper constants match independently normalized prior densities", {
  dat <- .evidence_test_fit()@stan_data
  dims <- .evidence_test_dims()
  answer <- .bayes_prior_offset(dat,dims)
  expected <- dnorm(0,sd=dat$prior_alpha_sd,log=TRUE)+2*dnorm(0,log=TRUE)+
    dnorm(1,1,dat$prior_sigma_meanval_strat,log=TRUE)-
    pnorm(1/dat$prior_sigma_meanval_strat,log.p=TRUE)+
    dnorm(dat$s0_prior_mean,dat$s0_prior_mean,dat$s0_prior_sd,log=TRUE)+
    dnorm(0,sd=dat$prior_log_sigma_share_sd,log=TRUE)+
    dnorm(0,sd=dat$prior_log_sigma_margin_sd,log=TRUE)
  expect_equal(answer$offset,expected,tolerance=1e-13)
  cut <- .bayes_prior_offset(dat,.evidence_test_dims("stackelberg"))$offset
  fringe <- .bayes_prior_offset(dat,.evidence_test_dims("fringe"))$offset
  expect_equal(cut-answer$offset,-lbeta(dat$prior_cutoff_alpha,dat$prior_cutoff_beta),
    tolerance=1e-13)
  expect_equal(fringe-cut,dnorm(0,sd=dat$prior_sigma_meanval_fringe,log=TRUE),
    tolerance=1e-13)
})

test_that("correlation constant agrees with quadrature and remains stable in tails", {
  dat <- .evidence_test_fit()@stan_data
  dims <- c(.evidence_test_dims(),list(rho_raw=1L))
  for(eta in c(.01,.1,.5,1,2,50,1e4)) {
    dat$prior_lkj <- eta
    answer <- .bayes_prior_offset(dat,dims)$components[["correlation"]]
    numerical <- -log(integrate(function(z)(1-z^2)^(eta-1),-.95,.95,
      rel.tol=1e-11)$value/.95)
    expect_equal(answer,numerical,tolerance=1e-10)
  }
  for(eta in c(1e-8,1e10)) {
    dat$prior_lkj <- eta
    expect_true(is.finite(.bayes_prior_offset(dat,dims)$offset))
  }
})

test_that("improper pools and active HMT have no proper evidence route", {
  dat <- .evidence_test_fit()@stan_data
  dat$use_s0_market_pool <- dat$use_s0_market_year_pool <- 1L
  expect_error(.bayes_prior_offset(dat,.evidence_test_dims()),"without a prior")
  dat$use_s0_market_pool <- dat$use_s0_market_year_pool <- 0L
  dat$use_hmt <- 1L
  expect_error(.bayes_prior_offset(dat,.evidence_test_dims()),"HMT")
  dat$likelihood_mode <- 2L
  expect_true(is.finite(.bayes_prior_offset(dat,.evidence_test_dims())$offset))
})

test_that("bridge correction changes odds and preserves raw legacy results", {
  local_mocked_bindings(extract_log_lik=function(...)array(0,c(5,1,1)),.package="loo")
  local_mocked_bindings(bridge_sampler=function(samples,...)list(logml=samples$raw),
    .package="bridgesampling")
  local_mocked_bindings(.bayes_evidence_offset=function(fit) {
    .bayes_prior_offset(fit@stan_data,.evidence_test_dims(fit@metadata$structure))
  },.package="antitrustBayes")
  sim <- .evidence_test_fit()
  fringe <- .evidence_test_fit("fringe",raw=-25)
  stack <- .evidence_test_fit("stackelberg",raw=-25)
  legacy <- suppressWarnings(bayes_compare(sim=sim,fringe=fringe))
  proper <- suppressWarnings(bayes_compare(sim=sim,fringe=fringe,
    bridge_normalization="proper"))
  automatic <- suppressWarnings(bayes_compare(sim=sim,stack=stack))
  expect_identical(legacy$bridge_normalization,rep("legacy",2))
  expect_equal(legacy$logml,c(-20,-25))
  expect_equal(proper$logml_raw,legacy$logml)
  expect_true(legacy$post_prob[1]>.99)
  expect_true(proper$post_prob[2]>.99)
  expect_identical(automatic$bridge_normalization,rep("proper",2))
  expect_true(automatic$post_prob[2]>.99)
  expect_equal(proper$logml-proper$logml_raw,proper$prior_normalization_offset)
  expect_equal(proper$logml_sd,legacy$logml_sd)
  corrected <- attr(proper,"corrected_bridge_logml")
  expect_equal(unname(corrected[[2]]),rep(proper$logml[2],3))
})

test_that("mixture LOO units contain whole equilibria", {
  fit <- .evidence_test_fit("stackelberg")
  expect_error(.bayes_comparison_units(list(fit),"observation"),"Observation LOO")
  expect_invisible(.bayes_comparison_units(list(fit),"period"))
  fit@stan_data$k_scope_code <- 1L
  expect_error(.bayes_comparison_units(list(fit),"period"),"Period LOO")
  fit@stan_data$cluster <- c(1L,2L)
  expect_error(.bayes_comparison_units(list(fit),"cluster"),"whole Stackelberg")
  expect_invisible(.bayes_comparison_units(list(fit),"market"))
})
