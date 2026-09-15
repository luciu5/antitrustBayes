# Test the independent R oracle itself before comparing it with new Stan code.
argv <- commandArgs(TRUE)
here <- dirname(sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1]))
out <- if(length(argv))argv[1] else here
source(file.path(here,"compiled_likelihood_reference.R"))
source(file.path(here,"synthetic_dgp.R"))
mock_data <- function(z,demand,conduct,prior,rho) {
  obs <- z$data; gs <- sort(unique(obs$group)); T <- length(z$truth)
  group <- match(obs$group,gs); tt <- as.integer(obs$year)
  p <- obs$prices; qc <- ave(obs$quantities,tt,FUN=function(x)x/sum(x))
  s0 <- z$truth[[1]]$outside_quantity
  sd_p <- sd(p); sd_lp <- sd(log(p)); pstd <- (p-mean(p))/sd_p
  lpstd <- (log(p)-mean(log(p)))/sd_lp
  bq <- vapply(seq_along(gs),function(gg) {
    ii <- which(group==gg)[1]
    if(demand=="logit") log(qc[ii])-qlogis(s0)+3*(p[ii]-mean(p)) else
      log(qc[ii])-qlogis(s0)+3.5*(log(p[ii])-mean(log(p))-log(z$parameters$outside_price))
  },numeric(1))
  bq <- bq+(seq_along(gs)-mean(seq_along(gs)))*.013
  obs$margins[2] <- NA_real_
  n <- length(gs); units <- if(prior=="cutoff")T else 1L
  dat <- list(N=nrow(obs),N_market_year=T,N_market=1L,N_group=n,
    market=rep(1L,nrow(obs)),group=group,market_year_idx=tt,
    price=p,shareIn=qc,log_shareIn=log(qc),price_std=pstd,price_sd=sd_p,
    log_price_std=lpstd,log_price_sd=sd_lp,market_sign=-1L,
    mky_to_market=rep(1L,T),outside_price=rep(z$parameters$outside_price,T),
    margin_obs=as.integer(is.finite(obs$margins)),
    log_marginNormInv=ifelse(is.finite(obs$margins),log(1/obs$margins),0),
    structure=2L,supply_model=if(conduct=="bertrand")1L else 3L,
    likelihood_mode=1L,use_rho=as.integer(!is.null(rho)),
    obs_order=order(tt,seq_len(nrow(obs))),cutoff_scale=.005,
    partition_prior_code=as.integer(prior=="ordered_k"),
    N_partition_unit=units,N_partition_component=units*n,
    mky_to_partition_unit=if(prior=="cutoff")seq_len(T) else rep(1L,T),
    partition_unit_start=seq(1L,by=n,length.out=units),
    partition_unit_count=rep(n,units),
    partition_component_k=rep(seq_len(n),units),
    partition_prior_logprob=rep(-log(n),units*n))
  dc <- list(s0_prob=rep(s0,T),b_quality_all=bq,b_market=0,
    a_market=if(demand=="logit")3*sd_p else 2.5*sd_lp,
    alpha_structural_market=3,gamma_structural_market=3.5,
    sigma_share_market=.15,sigma_margin_market=.07,nu=5,
    rho=if(is.null(rho))0 else rho,cutoff_share=.2)
  list(dat=dat,dc=dc,obs=obs)
}
rows <- list();i <- 0L
for(d in c("logit","ces")) for(cn in c("bertrand","cournot"))
  for(layout in c("single","multi_leader","multi_follower"))
    for(prior in c("cutoff","ordered_k")) for(rho in list(NULL,-.8,0,.8)) {
      z <- synthetic_stackelberg_panel(d,cn,1103,layout,normalize_bayes=TRUE)
      a <- mock_data(z,d,cn,prior,rho);dat <- a$dat;dc <- a$dc
      r <- reference_compiled_stack(dat,dc,d)
      cells <- vector("list",dat$N_market_year)
      for(tt in seq_len(dat$N_market_year)) {
        cell <- r$cells[[tt]];ix <- cell$ix
        cells[[tt]] <- boundary_component_reference(cell$size,.2,.005,d,cn,
          dat$price[ix],as.character(dat$group[ix]),cell$firm_share,
          if(d=="logit")dc$s0_prob[tt] else r$outside_revenue[tt],a$obs$margins[ix],
          alpha=3,gamma=3.5,sigma=.07*sqrt(3/5),df=5,rho=rho,
          share_error=dat$log_shareIn[ix]-log(r$predicted_share[ix]),sigma_share=.15,
          nonempty=TRUE)
      }
      if(prior=="cutoff") {
        expected <- sum(vapply(cells,function(x)x$target,numeric(1)))
        post <- unlist(lapply(cells,function(x)x$probabilities[-1]))
      } else {
        n <- dat$N_group
        joint <- ordered_k_component_reference(cells,
          setNames(rep(1/n,n),seq_len(n)))
        expected <- joint$target;post <- joint$probabilities
      }
      i <- i+1L
      rows[[i]] <- data.frame(demand=d,conduct=cn,layout=layout,prior=prior,
        rho=if(is.null(rho))"independent" else as.character(rho),
        target_error=abs(sum(r$unit_margin_loglik)-expected),
        probability_error=max(abs(r$component_probability-post)),
        accounting_error=abs(sum(r$margin_loglik)-expected))
    }
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,"compiled_likelihood_reference_unit_validation.csv"),row.names=FALSE)
stopifnot(max(tab$target_error)<1e-10,max(tab$probability_error)<1e-10,
          max(tab$accounting_error)<1e-10)
print(tab,digits=4)
cat("PASS: independent whole-unit oracle agrees with the separately coded cell/ordered-k reference in 96 panels.\n")
cat("This validates the reference; new compiled Stan agreement is a separate pending gate.\n")
