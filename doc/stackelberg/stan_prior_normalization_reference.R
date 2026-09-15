# Lead's independent normalized-prior reference for the existing version-6
# graphs and the approved Stackelberg copies of those graphs. No production
# sampler target is changed. Every returned evidence offset is constant.
stan_prior_normalization_reference <- function(dat,pars,dc) {
  lp <- constants <- numeric()
  add <- function(nm,density,constant) {
    lp[nm] <<- density; constants[nm] <<- constant
  }
  normal <- function(nm,x,mu=0,sd=1,sampled_sd=FALSE,positive=FALSE) {
    if(!length(x))return(invisible(NULL))
    n <- length(x)
    stopifnot(all(is.finite(x)),all(sd>0))
    norm <- if(positive) pnorm(mu/sd,log.p=TRUE) else 0
    add(nm,sum(dnorm(x,mu,sd,log=TRUE))-n*norm,
      -n*(.5*log(2*pi)+if(sampled_sd)0 else log(sd))-n*norm)
  }
  beta <- function(nm,x,a,b) {
    if(length(x))add(nm,sum(dbeta(x,a,b,log=TRUE)),-length(x)*lbeta(a,b))
  }
  p <- function(nm) if(is.null(pars[[nm]]))numeric() else as.numeric(pars[[nm]])
  normal("global_alpha",dc$mu_log_a,dat$prior_alpha_mean,dat$prior_alpha_sd)
  normal("alpha_scale",p("sigma_log_a_vec"),0,dat$prior_sigma_alpha,positive=TRUE)
  normal("alpha_raw",p("r_market_a_raw"))
  normal("alpha_centered",p("log_a_market_pool"),dc$mu_log_a,
    if(length(p("sigma_log_a_vec")))p("sigma_log_a_vec")[1] else 1,sampled_sd=TRUE)
  if(length(p("log_kappa_market"))) normal("alpha_cournot_centered",
    dc$log_a_market_derived,dc$mu_log_a,p("sigma_log_a_vec")[1],sampled_sd=TRUE)
  normal("market_mean",p("mu_b_market_vec"))
  normal("market_scale",p("sigma_b_market_vec"),0,.5,positive=TRUE)
  normal("market_raw",p("b_market_raw"))
  normal("quality_raw",p("b_group_raw"))
  normal("quality_scale",p("sigma_b_strat"),1,dat$prior_sigma_meanval_strat,positive=TRUE)
  beta("cutoff",p("cutoff_share"),dat$prior_cutoff_alpha,dat$prior_cutoff_beta)
  normal("fringe_quality",p("mu_b_fringe_vec"),0,dat$prior_sigma_meanval_fringe)
  normal("outside_mean",p("mu_s0_pool_vec"),dat$s0_prior_mean,dat$s0_prior_sd)
  normal("outside_scale",p("sigma_s0_pool_vec"),0,dat$prior_sigma_s0_pool,positive=TRUE)
  normal("outside_market_raw",p("s0_market_raw"))
  normal("outside_period_raw",p("s0_pool_raw"))
  sd0 <- if(length(p("sigma_s0_pool_vec")))p("sigma_s0_pool_vec")[1] else 1
  mu0 <- if(length(p("mu_s0_pool_vec")))p("mu_s0_pool_vec")[1] else 0
  normal("outside_market_centered",p("eta_s0_market"),mu0,sd0,sampled_sd=TRUE)
  normal("outside_period_centered",p("eta_s0_pool"),mu0,sd0,sampled_sd=TRUE)
  normal("outside_unpooled",p("s0_logit"),dat$s0_prior_mean,dat$s0_prior_sd)
  normal("outside_regression",p("beta_s0"),0,dat$prior_sigma_beta_s0)
  beta("bargaining",p("bargaining_unit"),dat$bargaining_prior_alpha,dat$bargaining_prior_beta)
  normal("share_scale",p("log_sigma_share"),log(dat$prior_sigma_share),dat$prior_log_sigma_share_sd)
  normal("margin_scale",p("log_sigma_margin"),log(dat$prior_sigma_margin),dat$prior_log_sigma_margin_sd)
  normal("share_market_scale_raw",p("log_sigma_share_market_raw"))
  normal("margin_market_scale_raw",p("log_sigma_margin_market_raw"))
  normal("share_market_scale_hyper",p("sigma_log_sigma_share_market"),
    0,dat$prior_sigma_resid_market,positive=TRUE)
  normal("margin_market_scale_hyper",p("sigma_log_sigma_margin_market"),
    0,dat$prior_sigma_resid_market,positive=TRUE)
  nu <- p("nu_minus_2")
  if(length(nu))add("df",sum(dgamma(nu,dat$prior_nu_shape,rate=dat$prior_nu_rate,log=TRUE)),
    length(nu)*(dat$prior_nu_shape*log(dat$prior_nu_rate)-lgamma(dat$prior_nu_shape)))
  rho <- p("rho_raw")
  if(length(rho)) {
    bound <- .95;eta <- dat$prior_lkj
    integ <- integrate(function(z)(1-z^2)^(eta-1),-bound,bound,
      rel.tol=1e-12,abs.tol=1e-12)$value/bound
    kernel <- (eta-1)*log1p(-dc$rho^2)+log1p(-tanh(rho)^2)
    add("correlation",kernel-log(integ),-log(integ))
  }
  # Constraint transformations are already included by log_prob(...,TRUE).
  # This reproduces that Jacobian independently, not an extra evidence offset.
  positive_names <- c("sigma_b_market_vec","sigma_b_strat","sigma_log_a_vec",
    "sigma_s0_pool_vec","sigma_log_sigma_share_market",
    "sigma_log_sigma_margin_market","nu_minus_2")
  jac <- sum(vapply(positive_names,function(nm)sum(log(p(nm))),numeric(1)))
  for(nm in c("cutoff_share","bargaining_unit")) {
    x <- p(nm);jac <- jac+sum(log(x)+log1p(-x))
  }
  list(log_prior=sum(lp),evidence_offset=sum(constants),
    prior_components=lp,offset_components=constants,constraint_jacobian=jac)
}

if(sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  here <- dirname(sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1]))
  out <- if(length(args))args[1] else here
  source(file.path(here,"compiled_likelihood_reference.R"))
  suppressPackageStartupMessages(library(antitrustBayes))
  suppressPackageStartupMessages(library(rstan))
  ns <- asNamespace("antitrustBayes");models <- get("stanmodels",ns)
  cases <- readRDS(file.path(here,"legacy_target_snapshot.rds"))
  rows <- list();i <- 0L
  for(case in cases) {
    prep <- get(if(case$demand=="logit")"prepare_stan_data" else "prepare_ces_stan_data",ns)
    dat <- do.call(prep,case$args)
    model <- models[[paste0(case$demand,"_supply")]]
    sf <- suppressWarnings(rstan::sampling(model,data=dat,chains=1L,iter=1L,warmup=0L,
      algorithm="Fixed_param",init=0,seed=42713L,refresh=0L))
    pn <- stan_parameter_names(model)
    for(j in seq_len(nrow(case$points))) {
      u <- case$points[j,]
      dc <- rstan::constrain_pars(sf,u)
      stopifnot(all(pn%in%names(dc)),"log_lik"%in%names(dc))
      ref <- stan_prior_normalization_reference(dat,dc[pn],dc)
      target <- rstan::log_prob(sf,u,adjust_transform=TRUE)
      recovered <- ref$log_prior+sum(dc$log_lik)+ref$constraint_jacobian-target
      i <- i+1L
      rows[[i]] <- data.frame(id=case$id,point=j,offset=ref$evidence_offset,
        recovered_offset=recovered,error=abs(recovered-ref$evidence_offset))
    }
  }
  tab <- do.call(rbind,rows)
  write.csv(tab,file.path(out,"legacy_prior_full_normalization_validation.csv"),row.names=FALSE)
  print(tab,digits=6)
  stopifnot(max(tab$error)<1e-9)
  cat("PASS: 72 full normalized-prior + GQ likelihood + Jacobian comparisons recover the constant compiled-target offset.\n")
}
