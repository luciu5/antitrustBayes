# Specification-7 density regression against archived version-6 targets.
# Constants are now native; the CES outside-revenue correction is a separately
# measured likelihood change, not falsely called an additive prior constant.
args<-commandArgs(TRUE)
here<-dirname(sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]))
out<-if(length(args))args[1] else here
source(file.path(here,'compiled_likelihood_reference.R'))
source(file.path(here,'stan_prior_normalization_reference.R'))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns<-asNamespace('antitrustBayes');mods<-get('stanmodels',ns)
cases<-readRDS(file.path(here,'legacy_target_density_snapshot.rds'))
rows<-list();ii<-0L
for(z in cases) {
 prep<-get(if(z$demand=='logit')'prepare_stan_data' else 'prepare_ces_stan_data',ns)
 dat<-do.call(prep,z$args);mod<-mods[[paste0(z$demand,'_supply')]]
 sf<-suppressWarnings(rstan::sampling(mod,data=dat,chains=1L,iter=1L,warmup=0L,algorithm='Fixed_param',init=0,refresh=0L,seed=49301L))
 pn<-stan_parameter_names(mod)
 for(j in seq_len(nrow(z$points))) {
  u<-z$points[j,];dc<-rstan::constrain_pars(sf,u)
  proper<-stan_prior_normalization_reference(dat,dc[pn],dc)
  ll<-sum(dc$log_lik);target<-rstan::log_prob(sf,u,adjust_transform=TRUE)
  expected<-z$target[j]+proper$evidence_offset+ll-z$likelihood[j]
  g<-rstan::grad_log_prob(sf,u,adjust_transform=TRUE)
  ii<-ii+1L
  rows[[ii]]<-data.frame(id=z$id,demand=z$demand,point=j,
   prior_constant=proper$evidence_offset,likelihood_change=ll-z$likelihood[j],
   target_shift_error=abs(target-expected),
   proper_density_error=abs(proper$log_prior+ll+proper$constraint_jacobian-target),
   unchanged_gradient_error=if(z$demand=='logit')max(abs(g-z$gradient[j,])) else NA_real_)
 }
}
tab<-do.call(rbind,rows)
write.csv(tab,file.path(out,'native_compiled_target_comparison.csv'),row.names=FALSE)
print(tab,digits=6)
stopifnot(nrow(tab)==72L,max(tab$target_shift_error)<1e-9,max(tab$proper_density_error)<1e-9,
 max(tab$unchanged_gradient_error,na.rm=TRUE)<1e-9,
 max(abs(tab$likelihood_change[tab$demand=='logit']))<1e-9)
cat('PASS: 72 native proper-density/target-shift cases; all Logit gradients unchanged. CES physical outside conversion produces the measured likelihood correction.\n')
