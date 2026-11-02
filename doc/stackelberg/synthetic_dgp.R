# Independent structural DGP, no package/equilibrium-solver dependency.
# Every baseline price is inverted from fixed demand quality, then audited
# hard-role FOCs recover constant costs. Roles are firm-level and publicly fixed.
synthetic_stackelberg_panel <- function(demand=c("logit","ces"),
                                       conduct=c("bertrand","cournot"),
                                       seed=1103L, layout=c("single","multi_leader","multi_follower"),
                                       normalize_bayes=FALSE,
                                       partition=c("cutoff","ordered_k"),k=2L,
                                       structure=c("stackelberg","fringe")) {
  demand<-match.arg(demand); conduct<-match.arg(conduct); layout<-match.arg(layout)
  partition<-match.arg(partition); structure<-match.arg(structure)
  set.seed(seed)
  firms<-LETTERS[1:5]
  owner<-switch(layout,single=firms,multi_leader=c("A","A",firms[-1]),multi_follower=c("A","B","B","C","D","E"))
  split<-rep(1,length(owner)); split[which(duplicated(owner)|duplicated(owner,fromLast=TRUE))]<-c(.55,.45)
  fi<-match(owner,firms); n<-length(owner)
  alpha<-3; gamma<-3.5; s0<-.3; p0<-1.3; cutoff<-.2; tau<-.005; noise_sd<-.04
  firm0<-c(.34,.22,.18,.14,.12); base_q<-firm0[fi]*split
  firm_ref<-2+seq_along(firms)*.1
  base_p<-if(demand=="logit") firm_ref[fi]-log(split)/alpha else firm_ref[fi]*split^(-1/gamma)
  qpaths<-t(vapply(1:8,function(tt) {
    fs<-c(.34,.22+.07*sin(2*pi*tt/8),.18+.025*cos(2*pi*tt/8),.14,NA_real_)
    fs[5]<-1-sum(fs[1:4]);fs[fi]*split
  },numeric(n)))
  pricepaths<-if(demand=="logit") {
    matrix(base_p,8,n,byrow=TRUE)-log(qpaths/matrix(base_q,8,n,byrow=TRUE))/alpha
  } else matrix(base_p,8,n,byrow=TRUE)*(matrix(base_q,8,n,byrow=TRUE)/qpaths)^(1/gamma)
  if(normalize_bayes) {
    # The old single-market Stan model centers firm quality and sets its
    # market intercept to zero. Derive compatible outside truth before fitting;
    # observed prices and conditional shares do not change with this choice.
    first<-!duplicated(owner)
    if(demand=="logit") {
      qnorm<-log(base_q)+alpha*(base_p-p0)
      eta<-mean(qnorm[first])+alpha*(p0-mean(pricepaths))
    } else {
      qnorm<-log(base_q)+gamma*log(base_p/p0)
      eta<-mean(qnorm[first])-gamma*mean(log(pricepaths))
    }
    s0<-plogis(eta)
  }
  kmin<-if(structure=="stackelberg")1L else 0L
  stopifnot(length(k)==1,is.finite(k),k==as.integer(k),k>=kmin,k<=length(firms))
  quality<-if(demand=="logit") log((1-s0)*base_q/s0)+alpha*(base_p-p0) else ((1-s0)*base_q/s0)*(base_p/p0)^gamma
  # Exact ordered-k estimation fixes one firm ranking before HMC and shares it
  # across the complete market panel. Reproduce that observable rule here:
  # equal-period mean firm conditional quantity share for Logit and conditional
  # revenue share for CES, followed by a firm-label tie break.
  size_paths<-matrix(NA_real_,8,length(firms),dimnames=list(NULL,firms))
  for(tt in 1:8) {
    fs<-c(.34,.22+.07*sin(2*pi*tt/8),.18+.025*cos(2*pi*tt/8),.14,NA_real_)
    fs[5]<-1-sum(fs[1:4]); qc<-fs[fi]*split; p<-pricepaths[tt,]
    rc<-qc*p/sum(qc*p)
    size_paths[tt,]<-if(demand=="logit") fs else
      as.numeric(tapply(rc,factor(owner,levels=firms),sum))
  }
  fixed_order<-order(-colMeans(size_paths),firms)
  fixed_core<-seq_along(firms)%in%head(fixed_order,k)
  records<-list(); truth<-list()
  for(tt in 1:8) {
    fs<-c(.34,.22+.07*sin(2*pi*tt/8),.18+.025*cos(2*pi*tt/8),.14,NA_real_)
    fs[5]<-1-sum(fs[1:4]); qc<-fs[fi]*split
    p<-pricepaths[tt,]
    expenditure_index<-s0*p0+(1-s0)*sum(qc*p)
    r0<-s0*p0/expenditure_index; rc<-qc*p/sum(qc*p)
    size<-size_paths[tt,]
    # Public hard-boundary lottery, conditioned on a nonempty leader sector.
    # The positive-k support revision precedes every new posterior fit.
    boundary<-cutoff+tau*qlogis(runif(1)*plogis((max(size)-cutoff)/tau))
    lf<-if(partition=="cutoff")size>boundary else fixed_core
    followers<-!lf
    R<-if(demand=="logit") fs*(1-s0) else size*(1-r0)
    if(demand=="logit") {
      if(structure=="fringe") {
        if(conduct=="bertrand")h<-ifelse(lf,1/(1-R),1)
        else h<-ifelse(lf,1+R/s0,1)
      } else if(conduct=="bertrand") {
        phi<-1/(1-sum(R[followers]^2/(1-R[followers]+R[followers]^2)))
        h<-ifelse(lf,1/(1-phi*R),1/(1-R))
      } else h<-ifelse(lf,1+R/(s0+sum(R[followers])),1+R/s0)
      m<-h[fi]/(alpha*p)
    } else {
      h<-gamma-1
      H<-1+h*r0; b<-h/(gamma*H)
      if(structure=="fringe") {
        if(conduct=="bertrand")mf<-ifelse(lf,1/(gamma-h*R),1/gamma)
        else mf<-ifelse(lf,1/gamma+b*R,1/gamma)
      } else if(conduct=="bertrand") {
        t<-1+h*(1-R); af<-h*R/(1-R+h*(1-R+R^2))
        phi<-1/(1-sum(R[followers]*af[followers]))
        mf<-ifelse(lf,1/(gamma-h*phi*R),1/t)
      } else {
        af<-h*(H-2*R+gamma*h*r0*R/H)/(H*(H+(gamma-2)*R))
        C<-sum(R[followers]*af[followers])
        mf<-ifelse(lf,1/gamma+b*R/(1+C),1/gamma+b*R)
      }
      m<-mf[fi]
    }
    noise<-rt(n,df=5)*sqrt(3/5)*noise_sd
    measured_margin<-m*exp(-noise)
    stopifnot(all(p>0),all(qc>0),all(m>0&m<1),all(measured_margin>0&measured_margin<1),abs(sum(qc)-1)<1e-14)
    costs<-p*(1-m)
    records[[tt]]<-data.frame(prices=p,quantities=1000*qc,margins=measured_margin,market="synthetic",year=as.character(tt),group=owner,product=paste0("p",seq_len(n)),stringsAsFactors=FALSE)
    truth[[tt]]<-list(prices=p,quantity_share=qc*(1-s0),revenue_share=rc*(1-r0),outside_quantity=s0,outside_revenue=r0,size=size,firm_ids=firms,core=firms[lf],leaders=if(structure=="stackelberg")firms[lf] else character(),boundary=if(partition=="cutoff")boundary else NA_real_,cutoff=if(partition=="cutoff")cutoff else NA_real_,k_hard=if(partition=="cutoff")sum(size>cutoff) else as.integer(k),k_realized=sum(lf),margins=m,costs=costs,quality=quality,fixed_rank=firms[fixed_order])
  }
  dat<-do.call(rbind,records)
  bq<-if(demand=="logit")quality+alpha*(p0-mean(dat$prices)) else log(quality)-gamma*mean(log(dat$prices))
  centered_error<-abs(mean(bq[!duplicated(owner)]))
  bq_centered<-bq-mean(bq[!duplicated(owner)])
  share_pred<-qlogis(s0)+if(demand=="logit")rep(bq_centered,8)-alpha*(dat$prices-mean(dat$prices)) else
    rep(bq_centered,8)-gamma*(log(dat$prices)-mean(log(dat$prices)))+gamma*log(p0)
  share_error<-max(abs(share_pred-log(as.vector(t(qpaths)))))
  if(normalize_bayes)stopifnot(share_error<1e-13,centered_error<1e-13)
  list(data=dat,truth=truth,demand=demand,conduct=conduct,structure=structure,seed=seed,layout=layout,
    parameters=list(alpha=alpha,gamma=gamma,s0=s0,outside_price=p0,
      cutoff=if(partition=="cutoff")cutoff else NA_real_,cutoff_scale=tau,
      k=if(partition=="ordered_k")as.integer(k) else NA_integer_,partition=partition,
      normalize_bayes=normalize_bayes,margin_noise_sd=noise_sd,margin_noise_df=5),
    normalization=list(share_equation_error=share_error,firm_center_error=centered_error))
}

if(sys.nframe()==0L) {
  args<-commandArgs(TRUE); out<-if(length(args))args[1] else "doc/stackelberg"
  bayes_mode<-length(args)>1&&args[2]=="bayes"
  configs<-expand.grid(demand=c("logit","ces"),conduct=c("bertrand","cournot"),seed=c(1103L,2207L,3301L),layout=c("single","multi_leader","multi_follower"),stringsAsFactors=FALSE)
  configs$partition<-"cutoff"
  if(bayes_mode) {
    other<-configs;other$partition<-"ordered_k";configs<-rbind(configs,other)
  }
  rows<-lapply(seq_len(nrow(configs)),function(i) {
    cc<-configs[i,]; z<-synthetic_stackelberg_panel(cc$demand,cc$conduct,cc$seed,cc$layout,
      normalize_bayes=bayes_mode,partition=cc$partition)
    # Independently confirm equal within-firm quality and exact outside-normalized demand.
    err<-vapply(z$truth,function(t) {
      pp<-t$prices; qq<-if(z$demand=="logit") exp(t$quality-z$parameters$alpha*(pp-z$parameters$outside_price)) else t$quality*(pp/z$parameters$outside_price)^(-z$parameters$gamma)
      uq<-qq/(1+sum(qq)); max(abs(uq-t$quantity_share))
    },numeric(1))
    qcomp<-vapply(split(z$truth[[1]]$quality,z$data$group[seq_along(z$truth[[1]]$quality)]),function(x) diff(range(x)),numeric(1))
    stopifnot(max(err)<1e-13,max(abs(qcomp))<1e-13)
    data.frame(cc,outside_quantity=z$parameters$s0,
      bayes_share_error=z$normalization$share_equation_error,
      firm_center_error=z$normalization$firm_center_error,
      min_cost=min(unlist(lapply(z$truth,`[[`,"costs"))),max_demand_error=max(err),min_k=min(vapply(z$truth,`[[`,integer(1),"k_realized")),max_k=max(vapply(z$truth,`[[`,integer(1),"k_realized")))
  })
  write.csv(do.call(rbind,rows),file.path(out,if(bayes_mode)"synthetic_bayes_normalization_validation.csv" else "synthetic_dgp_validation.csv"),row.names=FALSE)
  cat("PASS:",length(rows),"independent structural panel DGPs with fixed firm qualities.\n")
}
