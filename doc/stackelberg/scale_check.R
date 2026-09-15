# Prespecified scale sensitivity: probability geometry, not an estimated fit.
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "doc/stackelberg"
load("data/pnb.rda")
sizes <- sort(tapply(pnb$quantities,pnb$group,sum)/sum(pnb$quantities),decreasing=TRUE)
scales <- c(.0025,.005,.01,.02,.2)
cutoff <- .03
rows <- lapply(scales,function(tau) {
  vv<-plogis((sizes-cutoff)/tau)
  pp<-c(1-vv[1],head(vv,-1)-tail(vv,-1),tail(vv,1))
  stopifnot(all(pp>=0),abs(sum(pp)-1)<1e-14)
  kk<-0:length(sizes); cum<-cumsum(pp)
  data.frame(scale=tau,width_10_90=2*log(9)*tau,k_hard=sum(sizes>cutoff),k_mixture_mean=sum(pp*kk),k_mixture_lower=kk[which(cum>=.025)[1]],k_mixture_upper=kk[which(cum>=.975)[1]],largest_follower_prob=1-vv[1],smallest_follower_prob=1-tail(vv,1),near_cutoff_follower_prob=1-vv[which.min(abs(sizes-cutoff))])
})
tab<-do.call(rbind,rows); write.csv(tab,file.path(out,"scale_geometry.csv"),row.names=FALSE)
print(tab,digits=5)
# Reproducible prior sensitivity with the unchanged application-selectable
# default beta(3,100) cutoff prior. No historical 1--5 restriction is imposed.
set.seed(230913); cc<-rbeta(20000,3,100); eps<-qlogis(runif(length(cc)))
pr<-lapply(scales,function(tau) {
  kh<-vapply(cc,function(cut)sum(sizes>cut),integer(1))
  kr<-vapply(cc+tau*eps,function(cut)sum(sizes>cut),integer(1))
  data.frame(scale=tau,k_hard_mean=mean(kh),k_realized_mean=mean(kr),P_hard_k_1_5=mean(kh>=1&kh<=5),P_realized_k_1_5=mean(kr>=1&kr<=5))
})
write.csv(do.call(rbind,pr),file.path(out,"scale_prior_sensitivity.csv"),row.names=FALSE)
# Single-firm demonstration that a weighted-sector leader plug-in changes
# a hard-role-invariant monopoly markup.
S<-.35; w<-.5; B<-w*S^2/(1-S+S^2)
demo<-data.frame(hard_leader=1/(1-S),hard_follower=1/(1-S),coherent_mixture=1/(1-S),weighted_sector_leader=1/(1-S/(1-B)))
write.csv(demo,file.path(out,"weighted_sector_counterexample.csv"),row.names=FALSE)
stopifnot(abs(demo$hard_leader-demo$weighted_sector_leader)>1e-3)
cat("PASS: stable ordered-partition probability and prespecified scale checks.\n")
