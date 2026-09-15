# Independent compiled-Stan validation helpers. This reference reconstructs
# firm aggregates from product rows, rather than trusting production aggregates.
script_root <- function() dirname(sub("^--file=", "", commandArgs()[
  grepl("^--file=",commandArgs())][1]))
source(file.path(script_root(),"boundary_likelihood_reference.R"))

stan_parameter_names <- function(model) {
  code <- model@model_code
  block <- regmatches(code,regexec("(?s)\\nparameters\\s*\\{(.*?)\\n\\}",
                                  code,perl=TRUE))[[1]][2]
  if (is.na(block)) stop("Cannot locate the actual Stan parameter block")
  block <- gsub("//[^\\n]*", "",block)
  declarations <- strsplit(block,";",fixed=TRUE)[[1]]
  declarations <- trimws(declarations[nzchar(trimws(declarations))])
  vapply(declarations,function(z) {
    x <- regmatches(z,regexec("([A-Za-z_][A-Za-z_0-9]*)\\s*$",z,perl=TRUE))[[1]]
    if(length(x)!=2) stop("Unrecognized Stan parameter declaration: ",z)
    x[2]
  },character(1),USE.NAMES=FALSE)
}

stan_fixed_parameters <- function(sf,model) {
  draws <- rstan::extract(sf,permuted=TRUE)
  pn <- stan_parameter_names(model)
  ans <- lapply(pn,function(nm) {
    dims <- sf@par_dims[[nm]]
    a <- as.numeric(draws[[nm]])
    if(!length(dims)) a[1] else array(a[seq_len(prod(dims))],dim=dims)
  });names(ans) <- pn;ans
}

stan_one_draw <- function(sf) {
  draws <- rstan::extract(sf,permuted=TRUE)
  ans <- lapply(names(draws),function(nm) {
    dims <- sf@par_dims[[nm]]
    if(is.null(dims)||!length(dims)) as.numeric(draws[[nm]])[1] else
      array(as.numeric(draws[[nm]])[seq_len(prod(dims))],dim=dims)
  });names(ans)<-names(draws);ans
}

reference_compiled_stack <- function(dat,dc,demand) {
  stopifnot(dat$structure==2L,dat$supply_model%in%c(1L,3L))
  conduct <- if(dat$supply_model==1L) "bertrand" else "cournot"
  N <- dat$N; T <- dat$N_market_year
  pred_log_share <- numeric(N); share_ell <- numeric(N)
  outside_revenue <- numeric(T)
  cells <- vector("list",T)
  for(tt in seq_len(T)) {
    ix <- which(dat$market_year_idx==tt)
    mm <- dat$mky_to_market[tt]
    s0 <- dc$s0_prob[tt]
    conditional_q <- dat$shareIn[ix]
    avgp <- sum(dat$price[ix]*conditional_q)
    r0 <- plogis(qlogis(s0)+log(dat$outside_price[tt])-log(avgp))
    outside_revenue[tt] <- r0
    # Prefer raw quantity observations when the audit fixture supplies them.
    # Reconstructing revenue through already normalized quantity shares takes
    # a different floating-point path and magnifies one-ulp size differences
    # at tau=1e-7. Firm aggregation remains independent of production firm rows.
    raw_q <- attr(dat,"audit_raw_quantities",exact=TRUE)
    revenue_cond <- if (is.null(raw_q)) dat$price[ix]*conditional_q/avgp else {
      raw_rev <- dat$price[ix]*raw_q[ix]
      raw_rev/sum(raw_rev)
    }
    # Reconstruct unique firms from observations. No repeated firm row is summed.
    gs <- sort(unique(dat$group[ix]))
    size <- vapply(gs,function(gg) sum(if(demand=="logit")
      conditional_q[dat$group[ix]==gg] else revenue_cond[dat$group[ix]==gg]),numeric(1))
    fs <- size*(1-if(demand=="logit")s0 else r0)
    ord <- order(-size,gs); gs <- gs[ord]; fs <- fs[ord];size <- size[ord]
    names(fs)<-as.character(gs); names(size)<-as.character(gs)
    quality <- dc$b_quality_all[dat$group[ix]]
    if(demand=="logit") {
      delta <- dc$b_market[mm]+quality+
        dat$market_sign*dc$a_market[mm]*dat$price_std[ix]
    } else delta <- dc$b_market[mm]+quality+
      dat$market_sign*(dc$a_market[mm]+dat$log_price_sd)*dat$log_price_std[ix]+
      dc$gamma_structural_market[mm]*log(dat$outside_price[tt])
    pred_log_share[ix] <- qlogis(s0)+delta
    if(dat$likelihood_mode!=3L) {
      sc <- dc$sigma_share_market[mm]*sqrt((dc$nu-2)/dc$nu)
      share_ell[ix] <- dt((dat$log_shareIn[ix]-pred_log_share[ix])/sc,
                         df=dc$nu,log=TRUE)-log(sc)
    }
    predictions <- matrix(NA_real_,length(ix),length(gs)+1L)
    ell <- matrix(0,length(ix),length(gs)+1L)
    for(k in 0:length(gs)) {
      leaders <- if(k) as.character(gs[seq_len(k)]) else character()
      predictions[,k+1] <- boundary_margin_predictions(demand,conduct,
        dat$price[ix],as.character(dat$group[ix]),fs,
        if(demand=="logit")s0 else r0,leaders,
        alpha=if(demand=="logit")dc$alpha_structural_market[mm] else NULL,
        gamma=if(demand=="ces")dc$gamma_structural_market[mm] else NULL)
      usable <- dat$margin_obs[ix]==1L & dat$likelihood_mode!=2L
      active_rho <- dat$use_rho==1L && dat$likelihood_mode==1L && any(dat$margin_obs==1L)
      if(any(usable)) ell[usable,k+1] <- boundary_margin_density(
        dat$log_marginNormInv[ix[usable]],log(1/predictions[usable,k+1]+1e-9),
        dc$sigma_margin_market[mm],df=dc$nu,
        rho=if(active_rho)dc$rho else NULL,
        share_error=if(active_rho)dat$log_shareIn[ix[usable]]-pred_log_share[ix[usable]] else NULL,
        sigma_share=dc$sigma_share_market[mm])
    }
    cells[[tt]] <- list(ix=ix,size=size,firm_share=fs,firm_order=gs,
                         predictions=predictions,ell=ell)
  }
  C <- dat$N_partition_component; U <- dat$N_partition_unit
  lp <- score <- numeric(C); post <- numeric(C); unit_ell <- numeric(U)
  expected_margin <- expected_inverse <- numeric(N)
  leader_probability <- matrix(NA_real_,T,dat$N_group)
  sequential <- numeric(N)
  for(u in seq_len(U)) {
    cc <- seq.int(dat$partition_unit_start[u],length.out=dat$partition_unit_count[u])
    kk <- dat$partition_component_k[cc]
    stopifnot(all(kk>=1L))
    periods <- which(dat$mky_to_partition_unit==u)
    if(dat$partition_prior_code==0L) {
      stopifnot(length(periods)==1L)
      tt <- periods[1]
      lp[cc] <- boundary_nonempty_log_probabilities(cells[[tt]]$size,
        dc$cutoff_share[dat$mky_to_market[tt]],dat$cutoff_scale)[kk+1L]
    } else lp[cc] <- dat$partition_prior_logprob[cc]
    component_ell <- numeric(length(cc))
    for(tt in periods) component_ell <- component_ell+
      colSums(cells[[tt]]$ell[,kk+1L,drop=FALSE])
    score[cc] <- lp[cc]+component_ell
    unit_ell[u] <- reference_log_sum_exp(score[cc])
    post[cc] <- exp(score[cc]-unit_ell[u])
    for(tt in periods) {
      cell <- cells[[tt]]
      pp <- cell$predictions[,kk+1L,drop=FALSE]
      expected_margin[cell$ix] <- as.numeric(pp%*%post[cc])
      expected_inverse[cell$ix] <- as.numeric((1/pp)%*%post[cc])
      for(rr in seq_along(cell$firm_order)) {
        leader_probability[tt,cell$firm_order[rr]] <- sum(post[cc][kk>=rr])
      }
    }
    # Prefix increments merely account for the full unit likelihood.
    ix <- dat$obs_order[dat$market_year_idx[dat$obs_order]%in%periods]
    prefix <- lp[cc]; prev <- reference_log_sum_exp(prefix)
    for(j in ix) {
      tt <- dat$market_year_idx[j]
      cell <- cells[[tt]]
      prefix <- prefix+cell$ell[match(j,cell$ix),kk+1L]
      now <- reference_log_sum_exp(prefix)
      sequential[j] <- now-prev;prev <- now
    }
    stopifnot(abs(sum(sequential[ix])-unit_ell[u])<1e-10)
  }
  list(component_prior=lp,component_score=score,component_probability=post,
       unit_margin_loglik=unit_ell,margin_loglik=sequential,share_loglik=share_ell,
       loglik=sequential+share_ell,loglik_target=sum(unit_ell)+sum(share_ell),
       predicted_margin=expected_margin,predicted_inverse=expected_inverse,
       predicted_share=exp(pred_log_share),outside_revenue=outside_revenue,
       leader_probability=leader_probability,cells=cells)
}

if(sys.nframe()==0L) {
  suppressPackageStartupMessages(library(rstan))
  suppressPackageStartupMessages(library(antitrustBayes))
  ns <- asNamespace("antitrustBayes")
  dat <- get("prepare_stan_data",ns)(prices=c(1.8,2.1,2.4,2.2),
    quantities=c(35,20,30,15),margins=c(.22,.19,.18,.17),
    market=rep("m",4),year=rep(1,4),group=c("A","A","B","C"),
    supply_model=1L,priors=list(),use_hmt=FALSE,use_alpha_market_pool=FALSE,
    use_s0_market_year_pool=FALSE,use_s0_market_pool=FALSE)
  model <- get("stanmodels",ns)$logit_supply
  sf <- suppressWarnings(rstan::sampling(model,data=dat,chains=1L,iter=1L,
    warmup=0L,algorithm="Fixed_param",init=0,seed=1004L,refresh=0L))
  pars <- stan_fixed_parameters(sf,model)
  stopifnot(setequal(names(pars),stan_parameter_names(model)))
  u <- rstan::unconstrain_pars(sf,pars)
  stopifnot(all(is.finite(u)),is.finite(rstan::log_prob(sf,u,adjust_transform=TRUE)))
  dc <- stan_one_draw(sf)
  stopifnot(length(dc$pred_shareIn)==4L,length(dc$s0_prob)==1L)
  cat("PASS: independently parse actual parameter block and round-trip a compiled initialized target.\n")
}
