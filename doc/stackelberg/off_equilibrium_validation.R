# Actual profit derivatives away from equilibrium, including large market size.
# Usage: R_LIBS=... Rscript this.R coordination-path out-dir [logit|ces]
argv<-commandArgs(TRUE)
package_path<-if(length(argv))argv[1] else "/tmp/stackelberg-coordination"
out<-if(length(argv)>1)argv[2] else "doc/stackelberg"
demands<-if(length(argv)>2)argv[3] else c("logit","ces")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all(package_path,quiet=TRUE)
suppressPackageStartupMessages(library(nleqslv))
suppressPackageStartupMessages(library(numDeriv))
script<-sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
oracle<-new.env(parent=globalenv())
for(e in parse(file.path(dirname(script),"math_check.R")))
  if(is.call(e)&&identical(e[[1]],as.name("<-"))&&
     as.character(e[[2]])%in%c("cases","make_market"))eval(e,oracle)
rows<-list(); ii<-0L; size_multiplier<-1000
for(d in demands) for(cn in c("bertrand","cournot"))
 for(input in if(d=="logit")c(FALSE,TRUE) else FALSE)
  for(layout in c("single","multi_leader","multi_follower")) {
    cc<-oracle$cases[[layout]]
    for(stage in c("nash","leaders")) {
      cc$leaders<-if(stage=="nash")character() else unique(cc$owner)[1:2]
      mk<-oracle$make_market(cc,d,cn,input)
      obj<-coordination::stackelberg(mk$p,cc$r,mk$m,ownerPre=mk$owner,
        leadersPre=mk$leaders,demand=d,conduct=cn,output=!input,
        insideSize=size_multiplier*sum(cc$r),priceOutside=1.7,
        alpha=if(d=="logit")2.3 else NULL,gamma=if(d=="ces")4 else NULL)
      qbase<-mk$q*size_multiplier
      zbase<-if(cn=="bertrand")mk$p else qbase
      p_from_z<-function(z)if(cn=="bertrand")z else mk$pfun(z/size_multiplier)
      q_from_z<-function(z)if(cn=="bertrand")size_multiplier*mk$qfun(z) else z
      payoff<-function(z,firm) {
        own<-mk$owner==firm; pp<-p_from_z(z); qq<-q_from_z(z)
        sum((if(input)-1 else 1)*(pp[own]-mk$cost[own])*qq[own])
      }
      ngrad<-numeric(length(zbase))
      if(stage=="nash") {
        zz<-zbase*(1+seq_along(zbase)*.002)
        reported<-if(cn=="bertrand")coordination::stackelberg_residuals(obj,prices=zz)$all else
          coordination::stackelberg_residuals(obj,quantities=zz)$all
        for(firm in unique(mk$owner)) {
          own<-which(mk$owner==firm)
          fn<-function(x){z<-zz;z[own]<-x;payoff(z,firm)}
          ngrad[own]<-numDeriv::grad(fn,zz[own])
        }
      } else {
        zl<-zbase[mk$li]*1.04
        ff<-coordination::stackelberg_followers(obj,zl)
        zz<-ff$choices
        reported<-if(cn=="bertrand")coordination::stackelberg_residuals(obj,prices=zz)$leader else
          coordination::stackelberg_residuals(obj,quantities=zz)$leader
        ngrad<-numeric(length(mk$li))
        for(firm in mk$leaders) {
          own<-which(mk$owner==firm); pos<-match(own,mk$li)
          fn<-function(x){a<-zl;a[pos]<-x;payoff(coordination::stackelberg_followers(obj,a)$choices,firm)}
          ngrad[pos]<-numDeriv::grad(fn,zl[pos],method.args=list(eps=1e-4,d=1e-4))
        }
      }
      ii<-ii+1L; rows[[ii]]<-data.frame(demand=d,conduct=cn,input,layout,stage,
        maximum_gradient=max(abs(ngrad)),derivative_error=max(abs(ngrad-reported)))
    }
  }
tab<-do.call(rbind,rows)
write.csv(tab,file.path(out,"off_equilibrium_validation.csv"),row.names=FALSE)
print(aggregate(cbind(maximum_gradient,derivative_error)~demand+conduct+input+stage,tab,max),digits=5)
stopifnot(max(tab$derivative_error)<2e-7,all(tab$maximum_gradient>1e-6))
cat("PASS: actual nonzero profit-gradient sign, units and candidate-state checks.\n")
