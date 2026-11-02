# Published bayes_params already supported arbitrary positive cutoff_scale.
# Preserve these specified-fit paths independently of the newly appended
# estimation tuning argument. Run with baseline/new libraries in separate R.
# Rscript this.R baseline|compare baseline.rds output-directory
args <- commandArgs(TRUE);stopifnot(length(args)==3L)
mode <- args[1];snapshot <- args[2];out <- args[3]
suppressPackageStartupMessages(library(antitrustBayes))
ns <- asNamespace("antitrustBayes")
build <- get("build_conduct_object",ns)
owner <- c("A","A","B","C","C")
post <- owner;post[post=="C"] <- "B"
prices <- c(2.1,2.4,2.2,2.5,2.3)
quantities <- c(25,15,30,20,10)
results <- list()
for(demand in c("logit","ces"))
  for(side in if(demand=="logit")c("output","input") else "output")
    for(conduct in c(1L,3L))for(cutoff in c(FALSE,TRUE))
      for(scale in c(.05,.15,.2)) {
        slope <- if(demand=="logit")list(alpha=2.1) else list(gamma=3.5)
        pars <- do.call(bayes_params,c(slope,list(s0=.25,supply_model=conduct,
          use_cutoff=cutoff,cutoff_share=.15,cutoff_scale=scale)))
        market <- bayes_market(prices,quantities,rep(.15,5),owner,
          outside_price=1.7,market_side=side,insideSize=100)
        fit <- suppressWarnings(as_antitrust_fit(build(pars,1L,market,owner,owner)))
        model <- fit@model
        scenarios <- list(
          merger=antitrust::counterfactual(ownership=post),
          cost=antitrust::counterfactual(costs=c(0,0,.1,0,0)),
          combined=antitrust::counterfactual(ownership=post,costs=c(0,0,.1,0,0)))
        simulations <- lapply(scenarios,function(scenario)tryCatch({
          x <- suppressWarnings(antitrust::simulate(fit,scenario))
          list(price=x@pricePost,shares=antitrust::calcShares(x,FALSE),
            margins=antitrust::calcMargins(x,FALSE),costs=antitrust::calcMC(x,FALSE))
        },error=function(e)list(error=conditionMessage(e))))
        key <- paste(demand,side,conduct,cutoff,scale,sep="/")
        results[[key]] <- list(scale=pars$cutoff_scale,class=class(model),
          baseline=list(price=model@pricePre,shares=antitrust::calcShares(model,TRUE),
            margins=antitrust::calcMargins(model,TRUE),costs=antitrust::calcMC(model,TRUE)),
          simulations=simulations)
      }
if(mode=="baseline") {
  saveRDS(results,snapshot)
  cat("PASS: saved",length(results),"published specified-parameter scale fixtures.\n")
} else if(mode=="compare") {
  old <- readRDS(snapshot);stopifnot(identical(names(old),names(results)))
  differences <- vapply(names(results),function(k)
    isTRUE(all.equal(results[[k]],old[[k]],tolerance=1e-10,check.attributes=TRUE)),logical(1))
  write.csv(data.frame(id=names(results),unchanged=differences),
    file.path(out,"legacy_parameter_scale_comparison.csv"),row.names=FALSE)
  stopifnot(all(differences))
  cat("PASS: all",length(results),"legacy specified-scale fixtures unchanged.\n")
} else stop("mode must be baseline or compare")
