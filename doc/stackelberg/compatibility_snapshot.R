# Executable comparison of legacy structural realizations, no Stan fitting.
# Usage: Rscript this-file package-directory output-rds
args<-commandArgs(TRUE); stopifnot(length(args)==2)
pkgload::load_all(args[1],quiet=TRUE)
o<-c("A","A","B","C","C"); op<-o; op[op=="C"]<-"B"
p<-c(2.1,2.4,2.2,2.5,2.3); q<-c(25,15,30,20,10)
ans<-list()
for(demand in c("logit","ces")) for(side in if(demand=="logit") c("output","input") else "output") for(supply in 1:5) for(cutoff in c(FALSE,TRUE)) {
  if(supply==2&&cutoff) next
  pars<-if(demand=="logit") bayes_params(alpha=2.1,s0=.25,supply_model=supply,use_cutoff=cutoff,cutoff_share=.15) else bayes_params(gamma=3.5,s0=.25,supply_model=supply,use_cutoff=cutoff,cutoff_share=.15)
  mk<-bayes_market(prices=p,quantities=q,margins=rep(.15,5),group=o,market_side=side,insideSize=100,outside_price=1.7)
  ob<-build_conduct_object(pars,1,mk,o,o)
  ff<-suppressWarnings(as_antitrust_fit(ob))
  base<-ff@model
  sims<-lapply(list(merger=antitrust::counterfactual(ownership=op),cost=antitrust::counterfactual(costs=c(0,0,.1,0,0)),combined=antitrust::counterfactual(ownership=op,costs=c(0,0,.1,0,0))),function(scenario) {
    tryCatch({ x<-suppressWarnings(antitrust::simulate(ff,scenario)); list(pricePre=x@pricePre,pricePost=x@pricePost,sharePre=antitrust::calcShares(x,TRUE),sharePost=antitrust::calcShares(x,FALSE),mcPre=antitrust::calcMC(x,TRUE),mcPost=antitrust::calcMC(x,FALSE),marginPre=antitrust::calcMargins(x,TRUE),marginPost=antitrust::calcMargins(x,FALSE)) },error=function(e) list(error=conditionMessage(e)))
  })
  key<-paste(demand,side,supply,cutoff,sep="/")
  ans[[key]]<-list(class=class(base),baseline=list(price=base@pricePre,shares=antitrust::calcShares(base,TRUE),margins=antitrust::calcMargins(base,TRUE),costs=antitrust::calcMC(base,TRUE)),simulations=sims)
}
saveRDS(ans,args[2]); cat("Saved",length(ans),"legacy fixtures; errors retained for comparison.\n")
