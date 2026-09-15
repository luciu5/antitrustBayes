# Native one-firm/one-product array and generated-quantity dimension gate.
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns<-asNamespace('antitrustBayes');mods<-get('stanmodels',ns);rows<-list();i<-0L
for(demand in c('logit','ces'))for(conduct in c(1L,3L))
for(side in if(demand=='logit')c('output','input') else 'output')
for(prior in c('cutoff','ordered_k')) {
 aa<-list(prices=3,quantities=10,margins=.2,market='m',year='1',group='a',
  supply_model=conduct,market_side=side,use_cutoff=TRUE,use_rho=FALSE,use_hmt=FALSE,
  use_alpha_market_pool=FALSE,priors=list(),structure='stackelberg',partition_prior=prior)
 prep<-get(if(demand=='logit')'prepare_stan_data' else 'prepare_ces_stan_data',ns)
 dat<-do.call(prep,aa)
 sf<-suppressWarnings(rstan::sampling(mods[[paste0(demand,'_stackelberg')]],data=dat,
  chains=1L,iter=2L,warmup=0L,algorithm='Fixed_param',init=0,refresh=0,seed=3157L))
 stopifnot(sf@mode==0L)
 z<-rstan::extract(sf)
 stopifnot(length(z$posterior_component_prob)==2L,
  max(abs(z$posterior_component_prob-1))<1e-12,
  max(abs(z$firm_leader_prob-1))<1e-12,
  all(is.finite(z$pred_marginNormInv)),all(z$pred_marginNormInv>0),
  all(is.finite(z$log_lik)),all(is.finite(z$unit_margin_loglik)))
 i<-i+1L;rows[[i]]<-data.frame(demand,conduct,side,prior,N_firm=dat$N_firm_mky,
  N_component=dat$N_partition_component,finite_density=TRUE,role_error=0)
}
z<-do.call(rbind,rows)
filearg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1L])
write.csv(z,file.path(dirname(filearg),'compiled_singleton_validation.csv'),row.names=FALSE)
print(z);stopifnot(nrow(z)==12L)
cat('PASS: 12 native one-product/one-firm cells with two stored draws.\n')
