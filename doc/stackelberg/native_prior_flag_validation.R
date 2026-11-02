# Native specification-7 validation: compiled target equals full proper prior + likelihood + automatic Jacobian.
# Uses installed precompiled models; does not fit a posterior or change targets.
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1]))
out <- if(length(args))args[1] else here
source(file.path(here,"compiled_likelihood_reference.R"))
source(file.path(here,"stan_prior_normalization_reference.R"))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns <- asNamespace("antitrustBayes")
models <- get("stanmodels",ns)
cases <- readRDS(file.path(here,"legacy_target_snapshot.rds"))
variants <- list(
  default=list(),
  unpooled_alpha=list(use_alpha_market_pool=FALSE),
  centered_alpha=list(use_alpha_market_pool=TRUE,use_centered_alpha=TRUE),
  market_s0_raw=list(use_s0_market_pool=TRUE),
  market_s0_centered=list(use_s0_market_pool=TRUE,use_centered_s0_pool=TRUE),
  period_s0_raw=list(use_s0_market_year_pool=TRUE),
  period_s0_centered=list(use_s0_market_year_pool=TRUE,use_centered_s0_pool=TRUE),
  full_hierarchy=list(use_alpha_market_pool=TRUE,use_centered_alpha=TRUE,
    use_s0_market_year_pool=TRUE,use_centered_s0_pool=TRUE,
    use_market_resid_scale=TRUE,estimate_df=TRUE,
    X_s0=cbind(seq_len(4)/4,sin(seq_len(4)))),
  share_only=list(likelihood_mode=2L,estimate_df=TRUE),
  margin_only=list(likelihood_mode=3L,use_market_resid_scale=TRUE))
rows <- list(); ii <- 0L
for(demand in c("logit","ces"))for(conduct in c(1L,3L,5L)) {
  base <- cases[[which(vapply(cases,function(z)
    z$demand==demand && z$args$supply_model==1L &&
    isTRUE(z$args$use_cutoff) && isTRUE(z$args$use_rho) &&
    z$args$market_side=="output",logical(1)))[1]]]
  for(vn in names(variants)) {
    aa <- base$args; aa$supply_model <- conduct
    for(nm in names(variants[[vn]]))aa[[nm]] <- variants[[vn]][[nm]]
    prep <- get(if(demand=="logit")"prepare_stan_data" else "prepare_ces_stan_data",ns)
    dat <- do.call(prep,aa)
    # Non-default fixed hyperparameters exercise constants rather than only
    # their dimensions. Bargaining bounds act through a unit-Beta parameter.
    dat$prior_alpha_sd <- .7; dat$prior_sigma_meanval_strat <- .6
    dat$prior_sigma_meanval_fringe <- .8; dat$prior_sigma_resid_market <- .4
    dat$prior_lkj <- 1.4; dat$prior_nu_shape <- 2.3; dat$prior_nu_rate <- .7
    dat$bargaining_prior_alpha <- 2.2; dat$bargaining_prior_beta <- 3.1
    model <- models[[paste0(demand,"_supply")]]
    sf <- suppressWarnings(rstan::sampling(model,data=dat,chains=1L,iter=1L,
      warmup=0L,algorithm="Fixed_param",init=0,seed=78102L,refresh=0L))
    pars <- stan_fixed_parameters(sf,model)
    u0 <- rstan::unconstrain_pars(sf,pars)
    set.seed(971L+ii)
    points <- rbind(u0,u0+rnorm(length(u0),sd=.09),u0+rnorm(length(u0),sd=.2))
    for(j in seq_len(nrow(points))) {
      dc <- rstan::constrain_pars(sf,points[j,])
      ref <- stan_prior_normalization_reference(dat,dc[names(pars)],dc)
      target <- rstan::log_prob(sf,points[j,],adjust_transform=TRUE)
      recovered <- ref$log_prior+sum(dc$log_lik)+ref$constraint_jacobian-target
      ii <- ii+1L
      rows[[ii]] <- data.frame(demand=demand,conduct=conduct,variant=vn,
        point=j,offset=ref$evidence_offset,recovered_offset=recovered,
        error=abs(recovered))
    }
  }
}
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,"native_prior_flag_normalization_validation.csv"),row.names=FALSE)
print(aggregate(error~demand+conduct+variant,tab,max),digits=6)
stopifnot(nrow(tab)==180L,max(tab$error)<1e-9)
cat("PASS: 180 native normalized compiled prior/likelihood/Jacobian checks across pooling, centering, regression, residual-scale, DF, bargaining and likelihood-mode flags.\n")
