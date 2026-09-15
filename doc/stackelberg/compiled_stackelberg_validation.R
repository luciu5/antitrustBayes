# Fixed-parameter validation, not posterior estimation or recovery.
# Run only after the new programs are precompiled and installed.
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1]))
out <- if(length(args))args[1] else here
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(here,"compiled_likelihood_reference.R"))
source(file.path(here,"stan_prior_normalization_reference.R"))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(numDeriv))
ns <- asNamespace("antitrustBayes")
models <- get("stanmodels",ns)
if(!all(c("logit_stackelberg","ces_stackelberg")%in%names(models))) {
  stop("PENDING: corrected precompiled Stackelberg models are not installed")
}
fixture <- function(demand,conduct,prior,rho,side,variant="regular") {
  year <- rep(c(1,2,1,2),each=4)
  market <- rep(c("m1","m1","m2","m2"),each=4)
  group <- rep(c("A","A","B","C"),4)
  p <- rep(c(1.8,2.1,2.4,2.2),4)+seq_len(16)*.015
  q <- rep(c(35,20,30,15),4)
  if(variant=="ties") {
    p[seq(4,16,by=4)] <- p[seq(3,16,by=4)]
    q[seq(4,16,by=4)] <- q[seq(3,16,by=4)]
  }
  margins <- rep(c(.22,.19,NA,.17),4)
  mask <- rep(c(TRUE,TRUE,FALSE,TRUE),4)
  if(variant=="no_margins") {mask[]<-FALSE;margins[]<-NA_real_}
  idx <- as.vector(t(matrix(seq_len(16),nrow=4,byrow=TRUE)))
  call <- list(prices=p[idx],quantities=q[idx],margins=margins[idx],
    market=market[idx],year=year[idx],group=group[idx],
    margin_observed=mask[idx],supply_model=conduct,structure="stackelberg",
    partition_prior=prior,market_side=side,priors=list(),use_hmt=FALSE,use_cutoff=TRUE,
    use_rho=!is.null(rho),outside_price=c(1.3,1.3,1.7,1.7),
    use_alpha_market_pool=variant=="pooled",
    use_s0_market_year_pool=variant=="pooled",use_s0_market_pool=FALSE,
    use_market_resid_scale=variant=="pooled",estimate_df=variant=="pooled")
  if(prior=="cutoff")call$cutoff_scale<-if(variant=="tails")1e-7 else .01
  if(prior=="ordered_k")call$k_prior<-c("0"=0,"1"=.35,"2"=.65,"3"=0)
  if(variant=="share_only")call$likelihood_mode<-2L
  if(variant=="margin_only")call$likelihood_mode<-3L
  prep <- get(if(demand=="logit")"prepare_stan_data" else "prepare_ces_stan_data",ns)
  result <- do.call(prep,call)
  attr(result,"audit_raw_quantities") <- q[idx]
  result
}
fixed <- function(model,data,pars=0) suppressWarnings(rstan::sampling(model,
  data=data,chains=1L,iter=1L,warmup=0L,algorithm="Fixed_param",
  init=if(is.list(pars))list(pars) else pars,seed=49301L,refresh=0L))
set_param <- function(pars,nm,value) {
  if(nm%in%names(pars)&&length(pars[[nm]])) {
    if(length(value)==1L)value<-rep(value,length(pars[[nm]]))
    stopifnot(length(value)==length(pars[[nm]]))
    pars[[nm]][]<-value
  };pars
}
initialize <- function(model,dat,demand,rho,variant) {
  sf0 <- fixed(model,dat)
  pars <- stan_fixed_parameters(sf0,model)
  vals <- list(log_sigma_share=log(.15),log_sigma_margin=log(.07),
    sigma_b_strat=.4,b_group_raw=seq(-.12,.12,length.out=dat$N_group),
    mu_b_market_vec=.05,sigma_b_market_vec=.3,
    b_market_raw=seq(-.2,.3,length.out=dat$N_market),
    s0_logit=qlogis((.14-1e-6)/(1-2e-6)),
    mu_s0_pool_vec=qlogis((.14-1e-6)/(1-2e-6)),sigma_s0_pool_vec=.2,
    s0_pool_raw=seq(-.2,.2,length.out=dat$N_market_year),
    s0_market_raw=seq(-.2,.2,length.out=dat$N_market),
    sigma_log_a_vec=.15,
    r_market_a_raw=seq(-.1,.1,length.out=dat$N_market),
    sigma_log_sigma_share_market=.1,sigma_log_sigma_margin_market=.1,
    log_sigma_share_market_raw=seq(-.2,.2,length.out=dat$N_market),
    log_sigma_margin_market_raw=seq(.2,-.2,length.out=dat$N_market),
    nu_minus_2=5,
    cutoff_share=if(variant=="tails") .99 else c(.27,.3),
    rho_raw=if(is.null(rho))0 else atanh(rho/.95))
  for(nm in names(vals))pars<-set_param(pars,nm,vals[[nm]])
  sf1 <- fixed(model,dat,pars);dc1<-stan_one_draw(sf1)
  if(demand=="logit") {
    pars<-set_param(pars,"mu_log_a_free",log(3*dat$price_sd))
    pars<-set_param(pars,"log_kappa_vec",log(3*dat$price_sd)+mean(log(dc1$s0_prob)))
    means<-vapply(seq_len(dat$N_market),function(mm)
      mean(log(dc1$s0_prob[dat$mky_to_market==mm])),numeric(1))
    pars<-set_param(pars,"log_kappa_market",log(3*dat$price_sd)+means)
  } else pars<-set_param(pars,"mu_log_a",log(2.5*dat$log_price_sd))
  list(sf=fixed(model,dat,pars),pars=pars)
}
max_diff <- function(a,b) {
  stopifnot(length(a)==length(b))
  if(!length(a))return(0)
  finite <- is.finite(a)&is.finite(b)
  stopifnot(all(finite | (is.infinite(a)&is.infinite(b)&sign(a)==sign(b))))
  if(any(finite))max(abs(a[finite]-b[finite])) else 0
}
rows <- list();ii<-0L
config <- expand.grid(demand=c("logit","ces"),conduct=c(1L,3L),
  prior=c("cutoff","ordered_k"),rho=c("independent","-.8","0",".8"),
  side=c("output","input"),stringsAsFactors=FALSE)
config<-config[!(config$demand=="ces"&config$side=="input"),]
config$variant<-"regular"
extra <- expand.grid(demand=c("logit","ces"),conduct=c(1L,3L),
  prior=c("cutoff","ordered_k"),variant=c("pooled","ties","no_margins",
    "share_only","margin_only","tails"),stringsAsFactors=FALSE)
extra$rho<-".8";extra$side<-"output"
config<-rbind(config,extra[,names(config)])
for(jj in seq_len(nrow(config))) {
  cfg<-config[jj,];rho<-if(cfg$rho=="independent")NULL else as.numeric(cfg$rho)
  dat<-fixture(cfg$demand,cfg$conduct,cfg$prior,rho,cfg$side,cfg$variant)
  model<-models[[paste0(cfg$demand,"_stackelberg")]]
  zz<-initialize(model,dat,cfg$demand,rho,cfg$variant)
  sf<-zz$sf;dc<-stan_one_draw(sf)
  ref<-reference_compiled_stack(dat,dc,cfg$demand)
  # These canonical names are part of the signed API/GQ specification.
  required<-c("posterior_component_prob","prior_component_prob","component_log_score",
    "unit_margin_loglik","pred_margin","pred_marginNormInv","firm_leader_prob")
  stopifnot(all(required%in%names(dc)))
  target_data<-dat;target_data$likelihood_mode<-2L
  base_sf<-fixed(model,target_data)
  base_pars<-zz$pars
  if("rho_raw"%in%names(base_pars))base_pars$rho_raw<-numeric(0)
  # Share-only removes only the active correlation parameter/prior; every
  # other constrained parameter and its transform is exactly the same.
  main_target<-rstan::log_prob(sf,rstan::unconstrain_pars(sf,zz$pars),adjust_transform=TRUE)
  base_target<-rstan::log_prob(base_sf,rstan::unconstrain_pars(base_sf,base_pars),adjust_transform=TRUE)
  rho_kernel<-if(length(zz$pars$rho_raw)) {
    lognorm <- lbeta(.5,dat$prior_lkj)+pbeta(.95^2,.5,dat$prior_lkj,log.p=TRUE)-log(.95)
    (dat$prior_lkj-1)*log1p(-dc$rho^2)+log1p(-tanh(zz$pars$rho_raw)^2)-lognorm
  } else 0
  expected_diff<-sum(ref$unit_margin_loglik)+rho_kernel
  if(dat$likelihood_mode==3L)expected_diff<-expected_diff-
    sum(reference_compiled_stack(target_data,dc,cfg$demand)$share_loglik)
  # Check the unique firm data representation independently as well.
  aggregate_error<-0
  for(tt in seq_len(dat$N_market_year)) {
    cell<-ref$cells[[tt]]
    fm<-which(dat$firm_mky_market_year==tt)
    own<-dat$firm_mky_group[fm]
    expected<-cell$size[match(own,cell$firm_order)]
    actual<-if(cfg$demand=="logit")dat$firm_mky_shareIn[fm] else dat$firm_mky_revShare[fm]
    aggregate_error<-max(aggregate_error,max_diff(actual,expected))
  }
  fp<-ref$leader_probability[cbind(dat$firm_mky_market_year,dat$firm_mky_group)]
  gradient_relative<-NA_real_
  if((cfg$variant=="regular"&&cfg$rho==".8")||cfg$variant=="pooled") {
    u<-rstan::unconstrain_pars(sf,zz$pars)
    analytic<-rstan::grad_log_prob(sf,u,adjust_transform=TRUE)
    numeric<-numDeriv::grad(function(x)rstan::log_prob(sf,x,adjust_transform=TRUE),u,
      method="Richardson",method.args=list(eps=1e-4,d=1e-4))
    gradient_relative<-max(abs(analytic-numeric)/pmax(1,abs(analytic)))
    stopifnot(is.finite(gradient_relative),gradient_relative<2e-7)
  }
  proper <- stan_prior_normalization_reference(dat,zz$pars,dc)
  density_error <- abs(proper$log_prior+sum(ref$loglik)+proper$constraint_jacobian-main_target)
  ii<-ii+1L
  rows[[ii]]<-data.frame(cfg,
    density_error=density_error,
    target_error=abs(main_target-base_target-expected_diff),
    score_error=max_diff(dc$component_log_score,ref$component_score),
    probability_error=max_diff(dc$posterior_component_prob,ref$component_probability),
    prior_error=max_diff(dc$prior_component_prob,exp(ref$component_prior)),
    unit_error=max_diff(dc$unit_margin_loglik,ref$unit_margin_loglik),
    share_error=max_diff(dc$pred_shareIn,ref$predicted_share),
    inverse_margin_error=max_diff(dc$pred_marginNormInv,ref$predicted_inverse),
    margin_error=max_diff(dc$pred_margin,ref$predicted_margin),
    role_error=max_diff(dc$firm_leader_prob,fp),aggregate_error=aggregate_error,
    likelihood_error=max_diff(dc$log_lik,ref$loglik),
    outside_error=if(cfg$demand=="ces")max_diff(dc$R0_prob,ref$outside_revenue) else 0,
    gradient_relative=gradient_relative)
  print(rows[[ii]],digits=4)
}
tab<-do.call(rbind,rows)
write.csv(tab,file.path(out,"compiled_stackelberg_validation.csv"),row.names=FALSE)
errors<-grep("_error$",names(tab),value=TRUE)
stopifnot(max(as.matrix(tab[errors]))<1e-9)
cat("PASS: 96 fixed-parameter compiled cases match independent whole-unit targets, priors, predictions and posterior roles.\n")
cat("This is a deterministic density/GQ gate; it is not posterior recovery evidence.\n")
