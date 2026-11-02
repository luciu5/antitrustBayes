# Independent post-intervention primitive profit and follower-reaction checks.
# This research oracle uses no package demand, FOC, or follower solver code.
# Usage: R_LIBS=deps Rscript this.R coordination-path output-dir [logit|ces]
args <- commandArgs(TRUE)
path <- if(length(args)) args[1] else "/tmp/stackelberg-coordination"
out <- if(length(args)>1) args[2] else "doc/stackelberg"
demands <- if(length(args)>2) args[3] else c("logit","ces")
pkgload::load_all(path,quiet=TRUE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
script <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
oracle <- new.env(parent=globalenv())
for(e in parse(file.path(dirname(script),"math_check.R")))
  if(is.call(e) && identical(e[[1]],as.name("<-")) &&
     as.character(e[[2]]) %in% c("cases","make_market")) eval(e,oracle)
max0 <- function(x) if(length(x)) max(abs(x)) else 0
rows <- list(); ii <- 0L
for(demand in demands) for(conduct in c("bertrand","cournot"))
 for(input in if(demand=="logit") c(FALSE,TRUE) else FALSE)
  for(layout in c("single","two","multi_leader","multi_follower")) {
   cc <- oracle$cases[[layout]]; mk <- oracle$make_market(cc,demand,conduct,input)
   o <- coordination::stackelberg(mk$p,cc$r,mk$m,ownerPre=mk$owner,
     leadersPre=mk$leaders,demand=demand,conduct=conduct,output=!input,
     insideSize=1000*sum(cc$r),priceOutside=1.7,
     alpha=if(demand=="logit") 2.3 else NULL,gamma=if(demand=="ces") 4 else NULL,
     control.equ=list(implicitCheck=FALSE))
   for(shock in c("ownership","cost","combined")) {
    owner <- cc$owner; leaders <- cc$leaders
    delta <- numeric(length(owner))
    if(shock %in% c("cost","combined")) delta[owner=="B"] <- .1
    if(shock %in% c("ownership","combined")) {
      owner[owner %in% c("B","C")] <- "BC"
      leaders[leaders %in% c("B","C")] <- "BC"
      leaders <- unique(leaders)
    }
    s <- coordination::stackelberg_simulate(o,ownerPost=owner,
      leadersPost=leaders,mcDelta=delta)
    p <- as.numeric(s@pricePost); cst <- mk$cost*(1+delta)
    n <- length(p); sigma <- if(input) -1 else 1; M <- 1000; E <- 1000
    p0 <- 1.7; beta <- -sigma*2.3; g <- 4; h <- g-1
    if(demand=="logit") {
      quality <- log(cc$r/(1-sum(cc$r)))-beta*(mk$p-p0)
      qfun <- function(pp) {v<-exp(quality+beta*(pp-p0));M*v/(1+sum(v))}
      pfun <- function(qq) {
        if(any(qq<=0)||sum(qq)>=M) stop("quantity domain")
        p0+(log(qq/(M-sum(qq)))-quality)/beta
      }
      jfun <- function(pp) {qq<-qfun(pp);beta*(diag(qq)-tcrossprod(qq)/M)}
      kfun <- function(qq) (diag(1/qq)+matrix(1/(M-sum(qq)),n,n))/beta
    } else {
      quality <- cc$r/(1-sum(cc$r))*(mk$p/p0)^h
      qfun <- function(pp) E*quality*pp^(-g)/(p0^(-h)+sum(quality*pp^(-h)))
      pfun <- function(qq) {
        if(any(qq<=0)) stop("quantity domain")
        b <- (quality/qq)^(1/g)
        f <- function(w) p0^(-h)*w^g+w*sum(b*qq)-E
        hi <- 1; while(f(hi)<0) hi<-2*hi
        uniroot(f,c(0,hi),tol=1e-14)$root*b
      }
      jfun <- function(pp) {qq<-qfun(pp);r<-pp*qq/E;-g*diag(qq/pp)+h*outer(qq,r/pp)}
      kfun <- function(qq) {
        pp<-pfun(qq);r0<-1-sum(pp*qq)/E;b<-h/(g*(1+h*r0))
        -diag(pp/(g*qq))-b*tcrossprod(pp)/E
      }
    }
    q <- qfun(p); z <- if(conduct=="bertrand") p else q
    li <- which(owner %in% leaders); fi <- setdiff(seq_len(n),li)
    raw <- function(zz) {
      pp<-if(conduct=="bertrand") zz else pfun(zz)
      qq<-if(conduct=="bertrand") qfun(pp) else zz
      J<-if(conduct=="bertrand") jfun(pp) else kfun(qq)
      vapply(seq_len(n),function(j) {
        own<-which(owner==owner[j])
        sigma*(if(conduct=="bertrand") qq[j]+sum((pp[own]-cst[own])*J[own,j]) else
          pp[j]-cst[j]+sum(qq[own]*J[own,j]))
      },numeric(1))
    }
    profit <- function(zz,f) {
      pp<-if(conduct=="bertrand") zz else pfun(zz)
      qq<-if(conduct=="bertrand") qfun(pp) else zz
      sum(sigma*(pp[owner==f]-cst[owner==f])*qq[owner==f])
    }
    solve_f <- function(zl) {
      zz<-z;zz[li]<-zl
      if(length(fi)) {
        f<-function(v) {
          trial<-zz;trial[fi]<-exp(v)
          sc<-if(conduct=="bertrand") qfun(trial) else pfun(trial)
          raw(trial)[fi]/sc[fi]
        }
        sol<-nleqslv::nleqslv(log(z[fi]),f,method="Broyden",
          control=list(ftol=1e-14,xtol=1e-14,maxit=500))
        stopifnot(all(is.finite(sol$fvec)),max0(sol$fvec)<1e-11)
        zz[fi]<-exp(sol$x)
      }
      zz
    }
    Rfd <- if(length(fi)&&length(li)) numDeriv::jacobian(
      function(x) solve_f(x)[fi],z[li],method.args=list(eps=1e-4,d=1e-4)) else
      matrix(0,length(fi),length(li))
    Rp <- coordination::stackelberg_response(s,FALSE)
    leader_gradient <- 0
    for(firm in leaders) {
      own<-which(owner==firm);pos<-match(own,li)
      grad <- numDeriv::grad(function(x) {zl<-z[li];zl[pos]<-x;profit(solve_f(zl),firm)},
        z[own],method.args=list(eps=1e-4,d=1e-4))
      sc <- if(conduct=="bertrand") q[own] else p[own]
      leader_gradient <- max(leader_gradient,max0(grad/sc))
    }
    sc <- if(conduct=="bertrand") q else p
    ii<-ii+1L
    rows[[ii]] <- data.frame(demand,conduct,input,layout,shock,
      n_leaders=length(leaders),cost_error=max0(s@mcPost-cst),
      demand_quantity_error=max0(antitrust::calcQuantities(s,FALSE)-q),
      follower_foc=max0(raw(z)[fi]/sc[fi]),
      leader_gradient,independent_reaction_fd=max0(Rp-Rfd),
      independent_follower_actions=max0(solve_f(z[li])-z))
   }
  }
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,"policy_primitive_validation.csv"),row.names=FALSE)
print(aggregate(cbind(cost_error,demand_quantity_error,follower_foc,leader_gradient,
  independent_reaction_fd,independent_follower_actions)~demand+conduct+input,tab,max),digits=5)
stopifnot(max(tab$cost_error)<2e-9,max(tab$demand_quantity_error)<2e-8,
  max(tab$follower_foc)<2e-8,max(tab$leader_gradient)<2e-7,
  max(tab$independent_reaction_fd)<2e-7,max(tab$independent_follower_actions)<2e-7)
cat("PASS: independent post-merger, cost and combined-shock profit and reaction gates.\n")
