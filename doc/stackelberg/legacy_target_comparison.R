# Compare full compiled legacy posterior targets at identical parameter values.
# Run baseline in a separate R process with its installed library first.
# Usage: R_LIBS=baseline:deps Rscript this.R baseline snapshot.rds output-dir
#        R_LIBS=new:deps Rscript this.R compare snapshot.rds output-dir
args <- commandArgs(TRUE)
mode <- if(length(args)) args[1] else "baseline"
snapshot <- if(length(args)>1) args[2] else "legacy_target_snapshot.rds"
out <- if(length(args)>2) args[3] else "doc/stackelberg"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
ns <- asNamespace("antitrustBayes")
fixture <- function(demand,conduct,cutoff,rho,side) {
  # Interleave two periods and two markets, with one multiproduct firm.
  year <- rep(c(1,2,1,2),each=4)
  market <- rep(c("m1","m1","m2","m2"),each=4)
  group <- rep(c("A","A","B","C"),4)
  idx <- as.vector(t(matrix(seq_len(16),nrow=4,byrow=TRUE)))
  x <- list(prices=(rep(c(1.8,2.1,2.4,2.2),4)+seq_len(16)*.015)[idx],
    quantities=rep(c(35,20,30,15),4)[idx],
    margins=rep(c(.22,.19,NA,.17),4)[idx],market=market[idx],
    year=year[idx],group=group[idx],supply_model=conduct,
    margin_observed=rep(c(TRUE,TRUE,FALSE,TRUE),4)[idx],
    use_cutoff=cutoff,use_rho=rho,use_hmt=FALSE,market_side=side,priors=list(),
    outside_price=c(1.3,1.3,1.7,1.7))
  prep <- get(if(demand=="logit") "prepare_stan_data" else
    "prepare_ces_stan_data",ns)
  # The preparation API accepts one outside price per sorted market-year.
  list(args=x,data=do.call(prep,x),demand=demand)
}
instantiate <- function(fx) {
  model <- get("stanmodels",ns)[[paste0(fx$demand,"_supply")]]
  suppressWarnings(rstan::sampling(model,data=fx$data,chains=1L,iter=1L,
    warmup=0L,algorithm="Fixed_param",init=0,seed=48391L,refresh=0L))
}
if(mode=="baseline") {
  cases <- list(); i <- 0L
  for(demand in c("logit","ces")) for(conduct in c(1L,3L))
    for(cutoff in c(FALSE,TRUE)) for(rho in c(FALSE,TRUE))
      for(side in if(demand=="logit") c("output","input") else "output") {
        fx <- fixture(demand,conduct,cutoff,rho,side)
        sf <- instantiate(fx)
        dc <- rstan::extract(sf,permuted=TRUE)
        # Use the model's actual initialized constrained parameter names.
        pnames <- setdiff(names(sf@par_dims),"lp__")
        pars <- lapply(pnames,function(nm) {
          a <- dc[[nm]]; dims <- sf@par_dims[[nm]]
          if(!length(dims)) as.numeric(a)[1] else
            array(as.numeric(a)[seq_len(prod(dims))],dim=dims)
        }); names(pars) <- pnames
        u <- rstan::unconstrain_pars(sf,pars)
        set.seed(1820+i)
        points <- rbind(u,u+rnorm(length(u),sd=.12),u+rnorm(length(u),sd=.3))
        target <- apply(points,1,function(z) rstan::log_prob(sf,z,adjust_transform=TRUE))
        gradient <- t(apply(points,1,function(z) rstan::grad_log_prob(sf,z,adjust_transform=TRUE)))
        stopifnot(all(is.finite(target)),all(is.finite(gradient)))
        i <- i+1L
        cases[[i]] <- list(id=paste(demand,conduct,cutoff,rho,side,sep="_"),
          args=fx$args,demand=demand,points=points,target=target,
          gradient=gradient,parameter_dims=sf@par_dims)
      }
  saveRDS(cases,snapshot)
  cat("PASS: saved",length(cases),"compiled legacy fixtures at three parameter points each.\n")
} else if(mode=="compare") {
  cases <- readRDS(snapshot); rows <- vector("list",length(cases))
  for(i in seq_along(cases)) {
    old <- cases[[i]]
    prep <- get(if(old$demand=="logit") "prepare_stan_data" else "prepare_ces_stan_data",ns)
    sf <- instantiate(list(demand=old$demand,data=do.call(prep,old$args)))
    # New generated quantities may be added; every existing field retains
    # its dimensions and the unconstrained vector must still be accepted.
    stopifnot(identical(sf@par_dims[names(old$parameter_dims)],old$parameter_dims))
    target <- apply(old$points,1,function(z) rstan::log_prob(sf,z,adjust_transform=TRUE))
    gradient <- t(apply(old$points,1,function(z) rstan::grad_log_prob(sf,z,adjust_transform=TRUE)))
    rows[[i]] <- data.frame(id=old$id,target=max(abs(target-old$target)),
      gradient=max(abs(gradient-old$gradient)))
  }
  tab <- do.call(rbind,rows)
  write.csv(tab,file.path(out,"legacy_compiled_target_comparison.csv"),row.names=FALSE)
  print(tab,digits=5)
  stopifnot(max(tab$target)<1e-9,max(tab$gradient)<1e-9)
  cat("PASS: full legacy posterior targets and gradients are unchanged.\n")
} else stop("mode must be baseline or compare")
