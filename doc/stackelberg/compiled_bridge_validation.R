# Compiled fixed-draw posterior plumbing, not posterior recovery.
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns<-asNamespace('antitrustBayes');mods<-get('stanmodels',ns);rows<-list();i<-0L
for(demand in c('logit','ces'))for(conduct in c(1L,3L))
for(side in if(demand=='logit')c('output','input') else 'output')
for(prior in c('cutoff','ordered_k'))for(layout in c('multi_leader','multi_follower')) {
 owner<-if(layout=='multi_leader')c('A','A','B','C','D','E') else c('A','B','B','C','D','E')
 q<-if(layout=='multi_leader')c(30,20,20,10,12,8) else c(50,10,10,10,12,8)
 p<-seq(2,3,length.out=6);df<-data.frame(prices=rep(p,2)+rep(c(0,.03),each=6),
  quantities=rep(q,2),margins=.2,market='m',year=rep(c('1','2'),each=6),
  group=rep(owner,2),margin_observed=TRUE,outside_price=1,loo_cluster='m')
 aa<-c(as.list(df[c('prices','quantities','margins','market','year','group')]),
  list(supply_model=conduct,market_side=side,use_cutoff=TRUE,use_rho=FALSE,
   use_hmt=FALSE,use_alpha_market_pool=FALSE,priors=list(),structure='stackelberg',
   partition_prior=prior))
 prep<-get(if(demand=='logit')'prepare_stan_data' else 'prepare_ces_stan_data',ns)
 dat<-do.call(prep,aa)
 sf<-suppressWarnings(rstan::sampling(mods[[paste0(demand,'_stackelberg')]],data=dat,
  chains=1L,iter=2L,warmup=0L,algorithm='Fixed_param',init=0,refresh=0,seed=7153L))
 stopifnot(sf@mode==0L);z<-rstan::extract(sf)
 fit<-methods::new('BayesLogit',stanfit=sf,stan_data=dat,data=df,
  supply_model=conduct,demand_model=demand,market_side=side,
  price_sd=if(demand=='logit')dat$price_sd else dat$log_price_sd,
  call=call('bayes_est'),metadata=list(structure='stackelberg',partition_prior=prior,
   k_scope='market',cutoff_scale=.01,prior_normalization=list(
    mode='native_stan_normalized',model_spec_version=7L,partition_support='1:N')))
 sumry<-bayes_cutoff_summary(fit)$by_market_year
 stopifnot(nrow(sumry)==2L)
 summary_error<-0
 for(t in 1:2) {
  ids<-which(dat$firm_mky_market_year==t);u<-dat$mky_to_partition_unit[t]
  cc<-seq.int(dat$partition_unit_start[u],length.out=dat$partition_unit_count[u])
  prob<-colMeans(z$posterior_component_prob[,cc,drop=FALSE])
  sizes<-if(demand=='logit')dat$firm_mky_shareIn[ids] else dat$firm_mky_revShare[ids]
  core<-cumsum(as.numeric(sizes))[dat$partition_component_k[cc]]
  r<-which(sumry$year==as.character(t))
  summary_error<-max(summary_error,max(abs(sumry$realized_k_posterior[[r]]-prob)),
   abs(sumry$core_share_conditional[r]-sum(core*prob)),
   max(abs(sumry$firm_leader_probability[[r]]-colMeans(z$firm_leader_prob[,ids,drop=FALSE]))))
 }
 native<-bayes_stackelberg_model(fit,draw=2L,market='m',year='1',fixed_k=1L,
  insideSize=100,labels=paste0('sku',1:6))
 stopifnot(methods::is(native,if(demand=='logit')'StackelbergLogit' else 'StackelbergCES'),
  identical(native@leadersPre,'A'),is.list(native@diagnostics$bayes_reconstruction$normalization),
  native@diagnostics$bayes_reconstruction$normalization$model_spec_version==7L)
 inside_error<-abs(native@insideSize-if(demand=='ces')sum(p*q) else 100)
 # Posterior exact game selection is shared across both years under direct k.
 post_owner<-owner;post_owner[post_owner=='B']<-'A'
 costs<-setNames(.03,'sku6')
 sc<-antitrust::counterfactual(ownership=post_owner,costs=costs)
 set.seed(1709L)
 sim<-bayes_sim(fit,sc,n_draws=2L,market='m',year=c('1','2'),mc.cores=1L,
  insideSize=100,labels=paste0('sku',1:6))
 stopifnot(all(sim@n_valid==2L),sum(vapply(sim@failures,function(x)sum(!is.na(x)),integer(1)))==0L,
  identical(sim@metadata$simulation_route,'coordination::stackelberg_simulate'))
 if(prior=='ordered_k')stopifnot(all(vapply(sim@metadata$k_cache,function(k)length(unique(k))==1L,logical(1))))
 baseline_error<-max(vapply(names(sim@price_pre),function(key) {
  tt<-sub('^m_','',key);pp<-df$prices[df$year==tt]
  max(abs(sim@price_pre[[key]]-matrix(rep(pp,each=2),nrow=2)))
 },numeric(1)))
 normalized_error<-max(abs(unlist(lapply(sim@sims,function(uu)lapply(uu,function(obj)
  obj@diagnostics$residualsPost$normalizedAll)))))
 stopifnot(summary_error<1e-10,inside_error<1e-10,baseline_error<2e-7,
  normalized_error<2e-10)
 i<-i+1L;rows[[i]]<-data.frame(demand,conduct,side,prior,layout,summary_error,
  inside_error,baseline_error,normalized_error,n_valid=sum(sim@n_valid))
}
z<-do.call(rbind,rows)
filearg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1L])
write.csv(z,file.path(dirname(filearg),'compiled_bridge_validation.csv'),row.names=FALSE)
print(z,digits=5);stopifnot(nrow(z)==24L)
cat('PASS: 24 native fixed-draw summary/reconstruction/combined-shock panels; direct k shared across years.\n')
