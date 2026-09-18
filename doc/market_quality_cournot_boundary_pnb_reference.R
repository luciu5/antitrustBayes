# Independent PNB observed-anchor oracle; reads saved fits/plans only.
# No production solver, native target evaluation, HMC or bridge call.
suppressPackageStartupMessages(library(antitrustBayes))
args <- commandArgs(trailingOnly=TRUE)
mq <- normalizePath(if (length(args)) args[[1L]] else
  'doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality',
  mustWork=TRUE)
out <- if (length(args)>=2L) args[[2L]] else file.path(mq,
  paste0('cournot-oracle-replay-', format(Sys.time(), '%Y%m%dT%H%M%S')))
stopifnot(!dir.exists(out)); dir.create(out)
logsum <- function(x) {v<-max(x); v+log(sum(exp(x-v)))}
wlog <- function(l) {
  lo<-min(-50,l-1); hi<-max(1,log(max(1,l))+1)
  exp(uniroot(function(t) t+exp(t)-l,c(lo,hi),tol=1e-13)$root)
}
all <- list()
for (case in c('CF','CS')) {
  fit<-readRDS(file.path(mq,'pnb-archived-retry-01',paste0(case,'-fit.rds')))
  pol<-file.path(mq,'pnb-archived-retry-01/postprocess-run-01/analysis/policy')
  scenario<-readRDS(file.path(pol,paste0(case,'-policy-status.rds')))
  stopifnot(identical(as.character(scenario$acquirer),'298816'),
    identical(as.character(scenario$target),'889018'))
  plan<-readRDS(file.path(pol,paste0(case,'-draw-plan.rds')))
  ar<-readRDS(file.path(pol,paste0(case,'-observed-draws.rds')))
  stopifnot(identical(plan$draw_indices,ar$draw_indices),identical(fit@market_side,'input'))
  obs<-fit@data[match(ar$product_labels,as.character(fit@data$group)),,drop=FALSE]
  stopifnot(identical(as.character(obs$group),ar$product_labels),!anyDuplicated(obs$group))
  p<-obs$prices; r<-obs$quantities/sum(obs$quantities)
  owner<-as.character(obs$group); postowner<-owner
  postowner[postowner==as.character(scenario$target)]<-as.character(scenario$acquirer)
  ranked<-owner[order(-r,owner)]
  stopifnot(identical(ranked[1:3],c('354712','298816','889018')))
  draws<-rstan::extract(fit@stanfit,pars=c('a_market','s0_prob'))
  aa<-matrix(draws$a_market,ncol=1); ss<-matrix(draws$s0_prob,ncol=1)
  result<-lapply(seq_along(plan$draw_indices),function(pos) {
    index<-plan$draw_indices[pos]; a<-aa[index,1]/fit@price_sd; s0<-ss[index,1]
    k<-ar$k[pos]; roles<-ranked[seq_len(k)]
    ratios<-(1-s0)*r/s0
    hpre<-if(case=='CF') ifelse(owner %in% roles,1+ratios,1) else
      ifelse(owner %in% roles,1+ratios/(1+sum(ratios[!owner %in% roles])),1+ratios)
    values<-p+hpre/a; logaj<-log(ratios)+hpre
    firms<-unique(postowner)
    logaf<-setNames(vapply(firms,function(f)logsum(logaj[postowner==f]),numeric(1)),firms)
    postroles<-unique(postowner[owner %in% roles])
    xx<-setNames(numeric(length(firms)),firms)
    if(case=='CF') {
      xx[]<-vapply(firms,function(f)if(f %in% postroles)wlog(logaf[[f]]-1)else exp(logaf[[f]]-1),numeric(1))
      hh<-ifelse(firms %in% postroles,1+xx,1)
    } else {
      follows<-setdiff(firms,postroles)
      xx[follows]<-vapply(follows,function(f)wlog(logaf[[f]]-1),numeric(1))
      reserve<-1+sum(xx[follows])
      xx[postroles]<-vapply(postroles,function(f)reserve*wlog(logaf[[f]]-1-log(reserve)),numeric(1))
      hh<-ifelse(firms %in% postroles,1+xx/reserve,1+xx)
    }
    names(hh)<-firms; hpost<-hh[postowner]; prices<-values-hpost/a
    u<-log(ratios)+a*(prices-p); center<-max(c(0,u)); eu<-exp(u-center); e0<-exp(-center)
    shares<-eu/(e0+sum(eu)); outside<-e0/(e0+sum(eu))
    sx<-setNames(vapply(firms,function(f)sum(shares[postowner==f])/outside,numeric(1)),firms)
    independent_h<-if(case=='CF') ifelse(firms %in% postroles,1+sx,1) else
      ifelse(firms %in% postroles,1+sx/(1+sum(sx[!firms %in% postroles])),1+sx)
    names(independent_h)<-firms
    residual<-max(abs(prices-(values-independent_h[postowner]/a))/pmax(abs(prices),1))
    error<-if(is.na(ar$failures[pos]))max(abs((prices/p-1)-ar$price_change[pos,]))else NA_real_
    data.frame(case=case,anchor='observed',position=pos,draw_index=index,k=k,
      saved_failed=!is.na(ar$failures[pos]),min_price=min(prices),min_price_product=owner[which.min(prices)],
      target_price=prices[owner=='889018'],nonpositive=any(prices<=0),
      max_foc_residual=residual,max_success_price_change_difference=error,
      alpha=a,s0=s0)
  })
  tab<-do.call(rbind,result)
  write.csv(tab,file.path(out,paste0(case,'-oracle.csv')),row.names=FALSE)
  stopifnot(max(tab$max_foc_residual)<2e-10,
    max(tab$max_success_price_change_difference,na.rm=TRUE)<2e-9)
  print(table(saved_failed=tab$saved_failed,nonpositive=tab$nonpositive))
  print(tab[tab$saved_failed,c('case','k','min_price','min_price_product')])
  all[[case]]<-tab
}
joined<-do.call(rbind,all)
write.csv(joined,file.path(out,'oracle-all.csv'),row.names=FALSE)
print(aggregate(cbind(max_foc_residual,max_success_price_change_difference)~case,joined,max,na.rm=TRUE))
