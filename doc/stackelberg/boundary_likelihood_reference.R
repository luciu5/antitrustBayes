# Independent R reference for the approved marginalized timing likelihood.
# No equilibrium solver; these are audited observed-equilibrium margin moments.
boundary_log_probabilities <- function(size,cutoff,scale) {
  stopifnot(length(size)>0,all(diff(size)<=0),scale>0)
  n<-length(size); x<-(size-cutoff)/scale
  lp<-rep(-Inf,n+1)
  lp[1]<-plogis(-x[1],log.p=TRUE)
  lp[n+1]<-plogis(x[n],log.p=TRUE)
  if(n>1) for(k in 1:(n-1)) {
    gap<-(size[k+1]-size[k])/scale
    if(gap<0) lp[k+1]<-plogis(x[k],log.p=TRUE)+plogis(-x[k+1],log.p=TRUE)+log(-expm1(gap))
  }
  lp
}
boundary_nonempty_log_probabilities <- function(size,cutoff,scale) {
  lp <- boundary_log_probabilities(size,cutoff,scale)
  lp[1] <- -Inf
  lp[-1] <- lp[-1]-plogis((size[1]-cutoff)/scale,log.p=TRUE)
  lp
}
reference_log_sum_exp <- function(x) {
  mx<-max(x); if(!is.finite(mx)) return(mx)
  mx+log(sum(exp(x-mx)))
}
# sigma_margin is a marginal standard deviation, as in the package. A NULL
# rho selects independent univariate t errors; rho=0 still selects a joint t.
boundary_margin_density <- function(y,pred,sigma_margin,df=5,rho=NULL,
                                    share_error=NULL,sigma_share=NULL) {
  stopifnot(df>2,all(sigma_margin>0))
  if(is.null(rho)) {
    scale <- sigma_margin*sqrt((df-2)/df)
    return(dt((y-pred)/scale,df=df,log=TRUE)-log(scale))
  }
  stopifnot(abs(rho)<1,length(share_error)==length(y),all(sigma_share>0))
  mu <- pred+rho*sigma_margin/sigma_share*share_error
  scale <- sigma_margin*sqrt(1-rho^2)*
    sqrt((df-2+(share_error/sigma_share)^2)/(df+1))
  dt((y-mu)/scale,df=df+1,log=TRUE)-log(scale)
}
boundary_margin_predictions <- function(demand,conduct,prices,owner,
                                        firm_share,outside,leaders,alpha=NULL,gamma=NULL) {
  firms<-names(firm_share); R<-as.numeric(firm_share)
  lead<-firms%in%leaders; ff<-!lead
  if(demand=="logit") {
    if(conduct=="bertrand") {
      phi<-1/(1-sum(R[ff]^2/(1-R[ff]+R[ff]^2)))
      hf<-ifelse(lead,1/(1-phi*R),1/(1-R))
    } else hf<-ifelse(lead,1+R/(outside+sum(R[ff])),1+R/outside)
    hf[match(owner,firms)]/(alpha*prices)
  } else {
    h<-gamma-1
    if(conduct=="bertrand") {
      af<-h*R/(1-R+h*(1-R+R^2))
      phi<-1/(1-sum(R[ff]*af[ff]))
      mf<-ifelse(lead,1/(gamma-h*phi*R),1/(1+h*(1-R)))
    } else {
      H<-1+h*outside; b<-h/(gamma*H)
      af<-h*(H-2*R+gamma*h*outside*R/H)/(H*(H+(gamma-2)*R))
      cf<-sum(R[ff]*af[ff])
      mf<-ifelse(lead,1/gamma+b*R/(1+cf),1/gamma+b*R)
    }
    mf[match(owner,firms)]
  }
}
boundary_component_reference <- function(size,cutoff,scale,demand,conduct,
                                         prices,owner,firm_share,outside,margins,
                                         alpha=NULL,gamma=NULL,sigma=.04*sqrt(3/5),df=5,
                                         rho=NULL,share_error=NULL,sigma_share=NULL,
                                         nonempty=FALSE) {
  stopifnot(!is.null(names(size)),setequal(names(size),names(firm_share)))
  oo<-order(-size,names(size)); size<-size[oo]
  firms<-names(size); n<-length(firms)
  lp<-if(nonempty)boundary_nonempty_log_probabilities(size,cutoff,scale) else
    boundary_log_probabilities(size,cutoff,scale)
  pred<-matrix(NA_real_,length(prices),n+1)
  ell<-matrix(0,length(prices),n+1)
  usable<-is.finite(margins)&margins>0
  for(k in 0:n) {
    lead<-if(k) firms[seq_len(k)] else character()
    pred[,k+1]<-boundary_margin_predictions(demand,conduct,prices,owner,
       firm_share,outside,lead,alpha,gamma)
    ell[usable,k+1]<-boundary_margin_density(log(1/margins[usable]),
       log(1/pred[usable,k+1]+1e-9),sigma/sqrt((df-2)/df),df,rho,
       if(is.null(share_error))NULL else share_error[usable],sigma_share)
  }
  score<-lp+colSums(ell); target<-reference_log_sum_exp(score)
  probs<-exp(score-target)
  # Sequential accounting factors have an invariant sum, not individual LOO meaning.
  row_loglik<-numeric(length(prices)); prefix<-lp; prev<-reference_log_sum_exp(prefix)
  for(j in seq_along(prices)) {
    prefix<-prefix+ell[j,]; now<-reference_log_sum_exp(prefix)
    row_loglik[j]<-now-prev; prev<-now
  }
  list(target=target,log_prior=lp,probabilities=probs,predicted_margins=pred,
       component_log_density=ell,row_loglik=row_loglik,firm_order=firms)
}

ordered_k_component_reference <- function(cells,prior=NULL) {
  n<-min(vapply(cells,function(x)ncol(x$component_log_density)-1L,integer(1)))
  if(is.null(prior))prior<-setNames(rep(1/(n+1),n+1),0:n)
  kk<-as.integer(names(prior))
  stopifnot(length(kk)==length(prior),all(kk>=0&kk<=n),all(is.finite(prior)),
            all(prior>=0),abs(sum(prior)-1)<1e-12)
  ell<-do.call(rbind,lapply(cells,function(x)x$component_log_density[,kk+1,drop=FALSE]))
  lp<-log(prior); score<-lp+colSums(ell); target<-reference_log_sum_exp(score)
  probs<-exp(score-target)
  sequential<-numeric(nrow(ell)); prefix<-lp; prev<-reference_log_sum_exp(lp)
  for(j in seq_len(nrow(ell))) {
    prefix<-prefix+ell[j,]; now<-reference_log_sum_exp(prefix)
    sequential[j]<-now-prev;prev<-now
  }
  list(target=target,probabilities=probs,k=kk,row_loglik=sequential,
       component_log_density=ell,component_score=score)
}

if(sys.nframe()==0L) {
  argv<-commandArgs(TRUE)
  script<-sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
  source(file.path(dirname(script),"synthetic_dgp.R"))
  max_error<-0
  for(d in c("logit","ces")) for(cn in c("bertrand","cournot"))
    for(layout in c("single","multi_leader","multi_follower")) {
      zz<-synthetic_stackelberg_panel(d,cn,1103,layout)
      for(tt in seq_along(zz$truth)) {
        t<-zz$truth[[tt]]; obs<-zz$data[zz$data$year==as.character(tt),]
        fs<-tapply(if(d=="logit")t$quantity_share else t$revenue_share,obs$group,sum)
        sz<-setNames(t$size,t$firm_ids)
        fn<-function(order=seq_len(nrow(obs))) boundary_component_reference(sz,.2,.005,d,cn,
           obs$prices[order],obs$group[order],fs,if(d=="logit")t$outside_quantity else t$outside_revenue,
           obs$margins[order],alpha=3,gamma=3.5)
        r<-fn(); p<-fn(rev(seq_len(nrow(obs))))
        max_error<-max(max_error,abs(sum(exp(r$log_prior))-1),abs(sum(r$probabilities)-1),
                       abs(sum(r$row_loglik)-r$target),abs(r$target-p$target))
      }
    }
  tie<-boundary_log_probabilities(c(.4,.2,.2,.1),.2,.005)
  stopifnot(is.infinite(tie[3])&&tie[3]<0)
  for(cut in c(-1,0,.2,1,10)) for(sc in c(1e-7,.005,.2)) {
    lp<-boundary_log_probabilities(c(.4,.2,.2,.1),cut,sc)
    max_error<-max(max_error,abs(sum(exp(lp))-1))
  }
  # Single-firm timing cannot change monopoly incentives.
  for(d in c("logit","ces")) for(cn in c("bertrand","cournot")) {
    rr<-boundary_component_reference(c(A=.7),.4,.2,d,cn,2,"A",c(A=.6),.4,.3,alpha=3,gamma=3.5)
    max_error<-max(max_error,abs(diff(rr$predicted_margins)),abs(sum(rr$row_loglik)-rr$target))
  }
  # Independently compare conditional densities with a direct bivariate t
  # density divided by its share marginal, including rho=0 dependence.
  df<-5; ss<-.12; sm<-.07
  for(rho in c(-.85,0,.8)) for(e1 in c(-.2,0,.15)) for(e2 in c(-.1,0,.2)) {
    cov<-matrix(c(ss^2,rho*ss*sm,rho*ss*sm,sm^2),2,2)
    scale<-cov*(df-2)/df; e<-c(e1,e2)
    joint<-lgamma((df+2)/2)-lgamma(df/2)-log(df*pi)-
      .5*as.numeric(determinant(scale,logarithm=TRUE)$modulus)-
      (df+2)/2*log1p(as.numeric(crossprod(e,solve(scale,e)))/df)
    share_scale<-ss*sqrt((df-2)/df)
    marginal<-dt(e1/share_scale,df=df,log=TRUE)-log(share_scale)
    conditional<-boundary_margin_density(e2,0,sm,df,rho,e1,ss)
    max_error<-max(max_error,abs(conditional-(joint-marginal)))
  }
  # Correlated masked-row mixtures preserve joint scores under product order.
  zz<-synthetic_stackelberg_panel("ces","cournot",2207,"multi_follower")
  t<-zz$truth[[1]]; obs<-zz$data[zz$data$year=="1",]
  fs<-tapply(t$revenue_share,obs$group,sum); sz<-setNames(t$size,t$firm_ids)
  measured<-obs$margins; measured[2]<-NA_real_; err<-seq_len(nrow(obs))*.007-.02
  correlated<-function(order)boundary_component_reference(sz,.2,.005,"ces","cournot",
    obs$prices[order],obs$group[order],fs,t$outside_revenue,measured[order],gamma=3.5,
    rho=.6,share_error=err[order],sigma_share=.12)
  r<-correlated(seq_len(nrow(obs))); rp<-correlated(rev(seq_len(nrow(obs))))
  max_error<-max(max_error,abs(r$target-rp$target),abs(sum(r$row_loglik)-r$target),
                 abs(r$row_loglik[2]))
  no_margins<-boundary_component_reference(sz,.2,.005,"ces","cournot",obs$prices,
    obs$group,fs,t$outside_revenue,rep(NA_real_,nrow(obs)),gamma=3.5)
  max_error<-max(max_error,abs(no_margins$target),
                 max(abs(no_margins$probabilities-exp(no_margins$log_prior))))
  # One shared k across periods: logsumexp after adding all cell densities.
  cells<-lapply(1:2,function(tt) {
    t<-zz$truth[[tt]];obs<-zz$data[zz$data$year==as.character(tt),]
    fs<-tapply(t$revenue_share,obs$group,sum);sz<-setNames(t$size,t$firm_ids)
    boundary_component_reference(sz,.2,.005,"ces","cournot",obs$prices,obs$group,
      fs,t$outside_revenue,obs$margins,gamma=3.5)
  })
  prior<-c(`1`=.25,`2`=.5,`3`=.25)
  joint<-ordered_k_component_reference(cells,prior)
  separate<-sum(vapply(cells,function(x)ordered_k_component_reference(list(x),prior)$target,numeric(1)))
  manual<-reference_log_sum_exp(log(prior)+vapply(as.integer(names(prior)),function(k)
    sum(vapply(cells,function(x)sum(x$component_log_density[,k+1]),numeric(1))),numeric(1)))
  reversed<-ordered_k_component_reference(rev(cells),prior)
  max_error<-max(max_error,abs(joint$target-manual),abs(joint$target-reversed$target),
                 abs(sum(joint$row_loglik)-joint$target),abs(sum(joint$probabilities)-1))
  stopifnot(abs(joint$target-separate)>1e-6)
  stopifnot(max_error<2e-12)
  cat("PASS: marginalized boundary reference; probability, joint-target, ties, tail and permutation checks. Max error",max_error,"\n")
}
