suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
here<-dirname(sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]))
source(file.path(here,'compiled_likelihood_reference.R'))
cases <- readRDS(file.path(here,'legacy_target_snapshot.rds'))
ns<-asNamespace('antitrustBayes');mods<-get('stanmodels',ns)
for(i in seq_along(cases)) {
 z<-cases[[i]];prep<-get(if(z$demand=='logit')'prepare_stan_data' else 'prepare_ces_stan_data',ns)
 dat<-do.call(prep,z$args);mod<-mods[[paste0(z$demand,'_supply')]]
 sf<-suppressWarnings(rstan::sampling(mod,data=dat,chains=1L,iter=1L,warmup=0L,algorithm='Fixed_param',init=0,refresh=0L,seed=49301L))
 z$likelihood<-apply(z$points,1,function(u)sum(rstan::constrain_pars(sf,u)$log_lik))
 cases[[i]]<-z
}
saveRDS(cases,file.path(here,'legacy_target_density_snapshot.rds'))
cat('PASS: 72 archived compiled baseline likelihood evaluations.\n')
