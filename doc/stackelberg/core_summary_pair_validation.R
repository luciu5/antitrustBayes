# Independent explicit draw/component records; no flattened-index oracle.
suppressPackageStartupMessages(library(antitrustBayes))
m<-bayes_market(c(2,2,2),c(8,2,1),rep(.3,3),c('a','b','c'),outside_price=1,market_side='input')
rows<-list()
for(d in c('logit','ces')) {
 p<-if(d=='logit')bayes_params(alpha=c(.7,.7),s0=c(.1,.9),cutoff_share=c(.2,.8),cutoff_scale=.1,supply_model=1L,structure='stackelberg') else bayes_params(gamma=c(1.5,1.5),s0=c(.1,.9),cutoff_share=c(.2,.8),cutoff_scale=.1,supply_model=1L,structure='stackelberg')
 f<-bayes_specify(p,m);z<-bayes_cutoff_summary(f)$by_market_year;s<-c(8,2,1)/11;rec<-list();i<-0
 for(dd in 1:2) {
  cdf<-plogis((s-p$cutoff_share[dd])/.1)
  pr<-c(cdf[1]-cdf[2],cdf[2]-cdf[3],cdf[3])/cdf[1]
  out<-if(d=='logit')p$s0[dd] else p$s0[dd]/(p$s0[dd]+2*(1-p$s0[dd]))
  for(k in 1:3){i<-i+1;rec[[i]]<-data.frame(draw=dd,k=k,value=sum(s[1:k])*(1-out),weight=pr[k]/2)}
 }
 pairs<-do.call(rbind,rec);expected<-sum(pairs$value*pairs$weight)
 pairs<-pairs[order(pairs$value),]
 med<-pairs$value[which(cumsum(pairs$weight)>=.5-1e-12)[1]]
 rows[[d]]<-data.frame(demand=d,explicit_loop_mean=expected,
  production_mean=z$core_share_unconditional,explicit_sorted_median=med,
  production_median=z$core_share_unconditional_median)
 stopifnot(abs(expected-z$core_share_unconditional)<1e-12,
  abs(med-z$core_share_unconditional_median)<1e-12)
}
z<-do.call(rbind,rows);print(z,digits=13)
filearg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1L])
write.csv(z,file.path(dirname(filearg),'core_summary_pair_validation.csv'),row.names=FALSE)
cat('PASS: explicit independent pairs confirm joint outside/core means and weighted medians.\n')
