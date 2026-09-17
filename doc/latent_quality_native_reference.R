## Independent root check: existing native solvers accept fixed-cost latent demand.
## No coordination source changes and no Bayesian fitting.
suppressPackageStartupMessages(library(antitrust))
suppressPackageStartupMessages(library(coordination))
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else 'doc'
source(file.path(script_dir,'latent_quality_policy_reference.R'))
## Separate interior fixture: alpha=6 satisfies p_j > Cournot/core dollar
## markups here. The earlier alpha=2 reference is retained as noninterior failure.
alpha <- -6
lg <- latent_logit(x,s0,alpha,p,p0)
owner <- c('a','b','c')
post_owner <- c('a','a','c')
rows <- list()
for (demand in c('logit','ces')) for (conduct in c('bertrand','cournot')) {
  g <- 10
  latent <- if (demand == 'logit') lg else latent_ces_output(cx,s0,g,p,p0)
  os <- observed_s
  if (demand == 'ces') os <- p * observed_s / (p0*s0 + sum(p*observed_s))
  for (structure in c('stackelberg','fringe')) {
    common <- list(prices=p,shares=os,ownerPre=owner,demand=demand,
                   conduct=conduct,normIndex=NA,priceOutside=p0,
                   insideSize=1,labels=owner)
    if(demand=='logit') common$alpha <- abs(alpha) else common$gamma <- g
    if(structure=='stackelberg') {
      common$leadersPre <- 'a'
      obs <- do.call(coordination::stackelberg, common)
      residual <- coordination::stackelberg_residuals
      simulate <- coordination::stackelberg_simulate
    } else {
      common$corePre <- 'a'
      obs <- do.call(coordination::core_fringe, common)
      residual <- coordination::core_fringe_residuals
      simulate <- coordination::core_fringe_simulate
    }
    old_cost <- obs@mcPre
    old_size <- obs@mktSize
    obj <- obs
    obj@slopes$meanval <- latent$meanval
    predicted <- if(demand=='ces') {
      antitrust::calcShares(obj,preMerger=TRUE,revenue=TRUE)
    } else antitrust::calcShares(obj,preMerger=TRUE)
    expected <- if(demand=='logit') latent$s else latent$revenue
    prediction_error <- max(abs(predicted-expected))
    obj@pricePre <- antitrust::calcPrices(obj,preMerger=TRUE)
    pre_res <- residual(obj,TRUE)
    pre_error <- max(abs(pre_res$all))
    for (scenario in c('merger','cost','combined')) {
      args <- list(object=obj,ownerPost=if(scenario=='cost')owner else post_owner,
                   mcDelta=rep(if(scenario=='merger')0 else .1,length(p)))
      if(structure=='stackelberg') args$leadersPost <- 'a' else args$corePost <- 'a'
      out <- do.call(simulate,args)
      post_res <- residual(out,FALSE)
      expected_cost <- old_cost*(1+args$mcDelta)
      checks <- c(prediction_error,pre_error,max(abs(post_res$all)),
                  max(abs(obj@mcPre-old_cost)),max(abs(out@mcPre-old_cost)),
                  max(abs(out@mcPost-expected_cost)),abs(out@mktSize-old_size),
                  max(abs(out@slopes$meanval-obj@slopes$meanval)))
      rows[[length(rows)+1L]] <- data.frame(demand,conduct,structure,scenario,
        prediction_error,pre_error,post_error=max(abs(post_res$all)),
        max_invariant_error=max(checks[4:8]),
        baseline_price_discrepancy=max(abs(obj@pricePre-p)),
        passed=all(is.finite(checks))&&max(checks)<1e-7)
    }
  }
}
x<-do.call(rbind,rows)
write.csv(x,'/tmp/root_latent_quality_native_validation.csv',row.names=FALSE)
print(x,row.names=FALSE)
stopifnot(all(x$passed))
