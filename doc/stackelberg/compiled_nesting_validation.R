# Independent fixed-parameter endpoint and physical-share FOC gate.
args<-commandArgs(TRUE)
here<-dirname(sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]))
out<-if(length(args))args[1] else here
source(file.path(here,'compiled_likelihood_reference.R'))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns<-asNamespace('antitrustBayes');mods<-get('stanmodels',ns)
fixed<-function(mod,dat,pars=0) {
 sf<-suppressWarnings(rstan::sampling(mod,data=dat,chains=1L,iter=1L,warmup=0L,algorithm='Fixed_param',init=if(is.list(pars))list(pars) else pars,refresh=0L,seed=14217L))
 stopifnot(methods::is(sf,'stanfit'),sf@mode==0L)
 dc<-stan_one_draw(sf)
 stopifnot(all(c('pred_marginNormInv','log_lik')%in%names(dc)))
 sf
}
rows<-list();ii<-0L
for(demand in c('logit','ces'))for(conduct in c(1L,3L))
for(side in if(demand=='logit')c('output','input') else 'output')
for(multi in c(FALSE,TRUE))for(boundary in c(FALSE,TRUE)) {
 owner<-if(multi)c('A','A','B','C') else c('A','B','C')
 p<-if(multi)c(1.8,2.1,2.4,2.2) else c(1.8,2.4,2.2)
 q<-if(multi)c(35,20,30,15) else c(55,30,15)
 if(boundary) {
  p<-if(multi)c(1000,950,3,4) else c(1000,3,4)
  q<-if(multi)c(.6,.4,1e-11,2e-11) else c(1,1e-11,2e-11)
 }
 prep<-get(if(demand=='logit')'prepare_stan_data' else 'prepare_ces_stan_data',ns)
 aa<-list(prices=p,quantities=q,margins=rep(.2,length(p)),market=rep('m',length(p)),year=rep(1,length(p)),group=owner,supply_model=conduct,market_side=side,use_rho=FALSE,use_hmt=FALSE,use_alpha_market_pool=FALSE,priors=list(),outside_price=1)
 dsim<-do.call(prep,c(aa,list(structure='simultaneous',use_cutoff=FALSE)))
 dstk<-do.call(prep,c(aa,list(structure='stackelberg',use_cutoff=TRUE,partition_prior='ordered_k',k_prior=setNames(1,length(unique(owner))))))
 m1<-mods[[paste0(demand,'_supply')]];m2<-mods[[paste0(demand,'_stackelberg')]]
 init<-fixed(m1,dsim);pars<-stan_fixed_parameters(init,m1)
 # Use exactly the inverse bounded-logistic normalization.
 for(nm in c('s0_logit','mu_s0_pool_vec'))if(nm%in%names(pars)&&length(pars[[nm]]))pars[[nm]][]<-qlogis(((if(boundary)2e-6 else .14)-1e-6)/(1-2e-6))
 a<-fixed(m1,dsim,pars);b<-fixed(m2,dstk,pars)
 dc1<-stan_one_draw(a);dc2<-stan_one_draw(b)
 R0<-if(demand=='ces')dc1$R0_prob[1] else dc1$s0_prob[1]
 firm_cond<-tapply(if(demand=='ces')p*q/sum(p*q) else q/sum(q),owner,sum)
 S<-as.numeric(firm_cond[owner])*(1-R0)
 inv<-if(demand=='logit') {
  h<-if(conduct==1L)1/(1-S) else 1+S/R0
  p*dc1$alpha_structural_market[1]/h
 } else {
  g<-dc1$gamma_structural_market[1]
  if(conduct==1L)g-(g-1)*S else 1/(1/g+(g-1)*S/(g*(1+(g-1)*R0)))
 }
 ii<-ii+1L
 rows[[ii]]<-data.frame(demand,conduct,side,multi,boundary,outside=R0,
  inverse_error=max(abs(dc1$pred_marginNormInv-inv)),
  nesting_error=max(abs(dc1$pred_marginNormInv-dc2$pred_marginNormInv)),
  likelihood_error=max(abs(dc1$log_lik-dc2$log_lik)),
  role_error=max(abs(dc2$firm_leader_prob-1)),
  normalized_target_error=abs(
   rstan::log_prob(a,rstan::unconstrain_pars(a,stan_fixed_parameters(a,m1)),adjust_transform=TRUE)-
   rstan::log_prob(b,rstan::unconstrain_pars(b,stan_fixed_parameters(b,m2)),adjust_transform=TRUE)))
}
tab<-do.call(rbind,rows);write.csv(tab,file.path(out,'compiled_nesting_validation.csv'),row.names=FALSE)
print(tab,digits=5)
stopifnot(nrow(tab)==24L,max(tab$inverse_error)<1e-9,max(tab$nesting_error)<1e-9,
 max(tab$likelihood_error)<1e-9,max(tab$role_error)<1e-12,
 max(tab$normalized_target_error)<1e-9)
cat('PASS: 24 all-leader ordinary-conduct comparisons, including dominant firms and very small CES outside revenue shares, with independently computed exact physical-share FOCs and matching proper targets.\n')
