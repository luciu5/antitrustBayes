# Independent research oracle. This file is not package production code.
# Run: Rscript doc/stackelberg/math_check.R [output-directory]
suppressPackageStartupMessages(library(numDeriv))
suppressPackageStartupMessages(library(nleqslv))
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "doc/stackelberg"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  single = list(owner = c("A", "B", "C", "D"), leaders = "A", r = c(.25,.20,.15,.10)),
  two = list(owner = c("A", "B", "C", "D"), leaders = c("A","B"), r = c(.25,.20,.15,.10)),
  multi_leader = list(owner = c("A","A","B","C","D"), leaders = "A", r = c(.15,.10,.20,.15,.10)),
  multi_follower = list(owner = c("A","B","B","C","D"), leaders = "A", r = c(.25,.12,.08,.15,.10)),
  all = list(owner = c("A","B","C","D"), leaders = c("A","B","C","D"), r = c(.25,.20,.15,.10)),
  tiny = list(owner = c("A","B","C","D"), leaders = "A", r = c(.25,1e-4,2e-4,3e-4)),
  none = list(owner = c("A","B","C","D"), leaders = character(), r = c(.25,.20,.15,.10))
)

make_market <- function(case, demand, conduct, input = FALSE) {
  n <- length(case$r); p <- 2 + seq_len(n)*.13; r <- case$r
  o <- case$owner; firms <- unique(o); R <- as.numeric(tapply(r, factor(o, levels=firms), sum))
  lf <- firms %in% case$leaders; f <- !lf; fi <- which(o %in% firms[f]); li <- which(o %in% firms[lf])
  a <- 2.3; g <- if(is.null(case$gamma)) 4 else case$gamma
  r0 <- 1-sum(r); p0 <- 1.7; sign <- if(input) -1 else 1
  # 'sign' multiplies (p-c)q; for input c denotes marginal value, not cost.
  if (demand == "logit") {
    beta <- -sign*a; delta <- log(r/r0)-beta*(p-p0)
    qfun <- function(pp) { x <- exp(delta+beta*(pp-p0)); x/(1+sum(x)) }
    pfun <- function(qq) { if(any(qq<=0)||sum(qq)>=1) stop("domain"); p0+(log(qq/(1-sum(qq)))-delta)/beta }
    jfun <- function(pp) { q <- qfun(pp); beta*(diag(q)-tcrossprod(q)) }
    kfun <- function(qq) { (diag(1/qq)+matrix(1/(1-sum(qq)),n,n))/beta }
    if(conduct == "bertrand") {
      bf <- sum(R[f]^2/(1-R[f]+R[f]^2)); phi <- 1/(1-bf)
      h <- ifelse(lf, 1/(1-phi*R), 1/(1-R))
    } else {
      h <- ifelse(lf, 1+R/(r0+sum(R[f])), 1+R/r0)
    }
  } else {
    stopifnot(!input)
    av <- r/r0*(p/p0)^(g-1); av0 <- p0^(1-g)
    qfun <- function(pp) av*pp^(-g)/(av0+sum(av*pp^(1-g)))
    pfun <- function(qq) {
      if(any(qq<=0)) stop("domain")
      bb <- av^(1/g)*qq^(-1/g)
      fun <- function(tt) sum(tt*bb*qq)+av0*tt^g-1
      hi <- 1; while(fun(hi)<0) hi <- hi*2
      tt <- uniroot(fun,c(0,hi),tol=1e-13)$root; tt*bb
    }
    jfun <- function(pp) { q <- qfun(pp); rr <- pp*q; -g*diag(q/pp)+(g-1)*outer(q,rr/pp) }
    kfun <- function(qq) { pp <- pfun(qq); rr0 <- 1-sum(pp*qq); H <- 1+(g-1)*rr0; -diag(pp/qq)/g-(g-1)/(g*H)*tcrossprod(pp) }
    if(conduct == "bertrand") {
      t <- g-(g-1)*R; af <- (g-1)^2*R/(t*(t-1)+(g-1)^2*R)
      bf <- sum(R[f]*af[f]); phi <- 1/(1-bf)
      m <- ifelse(lf,1/(g-(g-1)*phi*R),1/t)
    } else {
      H <- 1+(g-1)*r0; b <- (g-1)/(g*H); mf <- 1/g+b*R; A <- 1-mf
      af <- g*b*(A-b*R+R*(g-1)^2*r0/H^2)/(A+(g-1)*b*R)
      cf <- sum(R[f]*af[f]); m <- ifelse(lf,1/g+b*R/(1+cf),mf)
    }
  }
  m <- if(demand == "logit") h[match(o,firms)]/a/p else m[match(o,firms)]
  cst <- p*(1-sign*m); q <- qfun(p); z <- if(conduct=="bertrand") p else q
  raw <- function(zz) {
    if(conduct=="bertrand") { pp<-zz; qq<-qfun(pp); J<-jfun(pp) }
    else { qq<-zz; pp<-pfun(qq); J<-kfun(qq) }
    ans<-numeric(n)
    for(k in seq_len(n)) {
      own<-which(o==o[k])
      ans[k]<-if(conduct=="bertrand") sign*(qq[k]+sum((pp[own]-cst[own])*J[own,k]))
              else sign*(pp[k]-cst[k]+sum(qq[own]*J[own,k]))
    }
    ans
  }
  profit <- function(zz, firm) { pp<-if(conduct=="bertrand") zz else pfun(zz); qq<-if(conduct=="bertrand") qfun(pp) else zz; ix<-o==firm; sum(sign*(pp[ix]-cst[ix])*qq[ix]) }
  solve_f <- function(zleader, start=z[fi]) {
    zz<-z; zz[li]<-zleader
    if(length(fi)) {
      transform<-function(x) if(conduct=="cournot") exp(x) else x
      ff<-function(x) { zz[fi]<-transform(x); scale<-if(conduct=="bertrand") qfun(zz)[fi] else pfun(zz)[fi]; raw(zz)[fi]/scale }
      sol<-nleqslv(if(conduct=="cournot") log(start) else start,ff,control=list(ftol=1e-15,xtol=1e-15,maxit=400))
      if(max(abs(sol$fvec))>2e-9) stop("Follower failed: ",max(abs(sol$fvec)))
      zz[fi]<-transform(sol$x)
    }
    zz
  }
  reaction <- matrix(0,length(fi),length(li))
  if(length(fi)&&length(li)) {
    if(demand=="logit"&&conduct=="bertrand") {
      reaction<-outer(R[match(o[fi],firms)]/(1-R[match(o[fi],firms)]+R[match(o[fi],firms)]^2)*phi, r[li])
    } else if(demand=="logit") reaction<- -outer(q[fi]/(r0+sum(R[f])),rep(1,length(li)))
    else if(conduct=="bertrand") reaction<-outer(p[fi]*af[match(o[fi],firms)]*phi, r[li]/p[li])
    else reaction<- -outer(q[fi]*af[match(o[fi],firms)]/(1+cf),p[li])
  }
  list(p=p,q=q,z=z,fi=fi,li=li,raw=raw,profit=profit,solve_f=solve_f,reaction=reaction,jfun=jfun,kfun=kfun,qfun=qfun,pfun=pfun,owner=o,leaders=case$leaders,cost=cst,m=m)
}

rows<-list(); optrows<-list(); ir<-0; io<-0
for(demand in c("logit","ces")) for(conduct in c("bertrand","cournot")) for(nm in names(cases)) {
  for(input in if(demand=="logit") c(FALSE,TRUE) else FALSE) {
    mk<-make_market(cases[[nm]],demand,conduct,input)
    zz<-mk$solve_f(mk$z[mk$li]); n<-length(zz)
    derivative_error<-0; implicit_error<-0; leader_error<-0; inverse_error<-max(abs(mk$kfun(mk$q)%*%mk$jfun(mk$p)-diag(n)))
    if(length(mk$fi)&&length(mk$li)) {
      fd<-jacobian(function(x) mk$solve_f(x)[mk$fi],mk$z[mk$li],method="Richardson",method.args=list(eps=1e-4,d=1e-4))
      J<-jacobian(mk$raw,mk$z)
      imp<- -solve(J[mk$fi,mk$fi,drop=FALSE],J[mk$fi,mk$li,drop=FALSE])
      derivative_error<-max(abs(fd-mk$reaction)); implicit_error<-max(abs(imp-mk$reaction))
    }
    for(l in unique(mk$owner)) {
      own<-which(mk$owner==l); positions<-match(own,mk$li)
      isleader<-l %in% mk$leaders
      prof<-if(isleader) function(x) { zl<-mk$z[mk$li]; zl[positions]<-x; mk$profit(mk$solve_f(zl),l) }
            else function(x) { zf<-mk$z; zf[own]<-x; mk$profit(zf,l) }
      gg<-grad(prof,mk$z[own],method.args=list(eps=1e-4,d=1e-4))
      if(isleader) leader_error<-max(leader_error,max(abs(gg)))
      eig<-max(eigen(hessian(prof,mk$z[own],method.args=list(eps=1e-4,d=1e-3)),symmetric=TRUE,only.values=TRUE)$values)
      # Independent bounded direct optimization of the reduced profit.
      if(conduct=="bertrand") {
        opt<-optim(mk$z[own]*1.025,function(x)-prof(x),method="L-BFGS-B",lower=mk$p[own]*.65,upper=mk$p[own]*1.5,control=list(factr=1e3,pgtol=1e-10,ndeps=rep(1e-6,length(own))))
      } else {
        # Keep inside the inverse-demand domain. Tiny changes are sufficient
        # for local-max verification; this is not a global existence proof.
        opt<-optim(mk$z[own]*1.025,function(x)-prof(x),method="L-BFGS-B",lower=mk$q[own]*.6,upper=mk$q[own]*1.4,control=list(factr=1e3,pgtol=1e-10,ndeps=rep(1e-7,length(own))))
      }
      io<-io+1; optrows[[io]]<-data.frame(demand,conduct,case=nm,input,firm=l,role=if(isleader) "leader" else "follower",action_error=max(abs(opt$par-mk$z[own])),profit_gain=-opt$value-prof(mk$z[own]),own_hessian_max=eig,convergence=opt$convergence)
    }
    # Verify ordinary follower Nash best responses independently as well.
    ferr<-if(length(mk$fi)) max(abs(mk$raw(zz)[mk$fi])) else 0
    ir<-ir+1; rows[[ir]]<-data.frame(demand,conduct,case=nm,input,follower_foc=ferr,reaction_fd=derivative_error,reaction_implicit=implicit_error,leader_gradient=leader_error,inverse_jacobian=inverse_error)
  }
}
tab<-do.call(rbind,rows); ot<-do.call(rbind,optrows)
write.csv(tab,file.path(out,"derivative_validation.csv"),row.names=FALSE)
write.csv(ot,file.path(out,"optimization_validation.csv"),row.names=FALSE)
print(aggregate(cbind(follower_foc,reaction_fd,reaction_implicit,leader_gradient,inverse_jacobian)~demand+conduct+input,tab,max),digits=5)
print(aggregate(cbind(action_error,profit_gain)~demand+conduct+input,ot,max),digits=5)
stopifnot(max(tab$follower_foc)<2e-9,max(tab$reaction_fd)<2e-7,max(tab$reaction_implicit)<2e-7,max(tab$leader_gradient)<2e-7,max(tab$inverse_jacobian)<2e-9,max(ot$action_error)<2e-5,max(ot$profit_gain)<1e-9,max(ot$own_hessian_max)<0)
cat("PASS: hard-game mathematical audit; local optimization and derivative gates.\n")
