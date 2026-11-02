# Additional independent mathematical checks where Cournot reactions change sign.
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "doc/stackelberg"
script <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
oracle <- new.env(parent=globalenv())
suppressPackageStartupMessages(library(numDeriv))
suppressPackageStartupMessages(library(nleqslv))
for(e in parse(file.path(dirname(script),"math_check.R")))
 if(is.call(e)&&identical(e[[1]],as.name("<-"))&&
    as.character(e[[2]]) %in% c("cases","make_market"))eval(e,oracle)
rows <- list(); i <- 0L
for(g in c(1.1,1.5,2,6)) {
  cc <- list(owner=c("A","B","C","D"),leaders="A",
             r=c(.03,.82,.08,.04),gamma=g)
  mk <- oracle$make_market(cc,"ces","cournot")
  nr <- numDeriv::jacobian(function(q)mk$solve_f(q)[mk$fi],mk$q[mk$li],
                          method.args=list(eps=1e-4,d=1e-4))
  jj <- numDeriv::jacobian(mk$raw,mk$q)
  ri <- -solve(jj[mk$fi,mk$fi,drop=FALSE],jj[mk$fi,mk$li,drop=FALSE])
  objective <- function(q)mk$profit(mk$solve_f(q),"A")
  derivative <- numDeriv::grad(objective,mk$q[mk$li])
  curvature <- numDeriv::hessian(objective,mk$q[mk$li])[1,1]
  opt <- optimize(function(lq)-objective(exp(lq)),
                  log(mk$q[mk$li])+log(c(.5,2)),tol=1e-11)
  i <- i+1L
  rows[[i]] <- data.frame(gamma=g,dominant_follower_response=mk$reaction[1,1],
    finite_difference_error=max(abs(nr-mk$reaction)),
    implicit_error=max(abs(ri-mk$reaction)),leader_gradient=abs(derivative),
    leader_curvature=curvature,optimization_relative_action_error=
      abs(exp(opt$minimum)/mk$q[mk$li]-1))
}
tab<-do.call(rbind,rows)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(tab,file.path(out,"ces_negative_response_validation.csv"),row.names=FALSE)
print(tab,digits=5)
stopifnot(any(tab$dominant_follower_response>0),max(tab$finite_difference_error)<2e-7,
  max(tab$implicit_error)<2e-7,max(tab$leader_gradient)<2e-7,
  all(tab$leader_curvature<0),max(tab$optimization_relative_action_error)<2e-5)
cat("PASS: negative CES Cournot reaction coefficients, IFT and direct leader optimization.\n")
