# Independent package validation against the previously signed research oracle.
# Usage: R_LIBS=/tmp/stackelberg-lib Rscript this.R coordination-path out-dir
argv <- commandArgs(TRUE)
package_path <- if(length(argv)) argv[1] else "/tmp/stackelberg-coordination"
out_dir <- if(length(argv)>1) argv[2] else "doc/stackelberg"
demand_filter <- if(length(argv)>2) argv[3] else c("logit","ces")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(nleqslv))
suppressPackageStartupMessages(library(numDeriv))
pkgload::load_all(package_path,quiet=TRUE)
script_path <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
oracle <- new.env(parent=globalenv())
# Load only the market definitions, without rerunning the signed audit suite.
for(e in parse(file.path(dirname(script_path),"math_check.R"))) {
  if(is.call(e)&&identical(e[[1]],as.name("<-"))&&
     as.character(e[[2]])%in%c("cases","make_market")) eval(e,oracle)
}
max0 <- function(x) if(length(x)) max(abs(x)) else 0
rows <- list(); ii <- 0L
for(demand in demand_filter) for(conduct in c("bertrand","cournot"))
  for(case_name in names(oracle$cases)) for(input in if(demand=="logit") c(FALSE,TRUE) else FALSE) {
    cc <- oracle$cases[[case_name]]
    mk <- oracle$make_market(cc,demand,conduct,input)
    obj <- coordination::stackelberg(prices=mk$p,shares=cc$r,margins=mk$m,
       ownerPre=mk$owner,leadersPre=mk$leaders,demand=demand,conduct=conduct,
       output=!input,insideSize=sum(cc$r),priceOutside=1.7,
       alpha=if(demand=="logit")2.3 else NULL,gamma=if(demand=="ces")4 else NULL)
    ra <- coordination::stackelberg_response(obj,method="analytic")
    ri <- coordination::stackelberg_response(obj,method="implicit")
    sf <- coordination::stackelberg_followers(obj,mk$z[mk$li])
    rf <- if(length(mk$fi)&&length(mk$li)) numDeriv::jacobian(
      function(x) coordination::stackelberg_followers(obj,x)$choices[mk$fi],
      mk$z[mk$li],method.args=list(eps=1e-4,d=1e-4)) else ra
    tmp <- obj; tmp@priceStart <- mk$p*1.015
    pg <- antitrust::calcPrices(tmp,method="implicit")
    resid <- coordination::stackelberg_residuals(obj)
    # Raw oracle profit gradients of actual leader choices, followers fully solved.
    lg <- 0
    for(firm in mk$leaders) {
      own <- which(mk$owner==firm); pos <- match(own,mk$li)
      fun <- function(x) {
        zl<-mk$z[mk$li]; zl[pos]<-x
        zz<-coordination::stackelberg_followers(obj,zl)$choices
        mk$profit(zz,firm)
      }
      lg <- max(lg,max0(numDeriv::grad(fun,mk$z[own],method.args=list(eps=1e-4,d=1e-4))))
    }
    ii<-ii+1L
    rows[[ii]]<-data.frame(demand,conduct,case=case_name,input,
       baseline_price=max0(obj@pricePre-mk$p),baseline_quantity=max0(antitrust::calcQuantities(obj)-mk$q),
       cost_inversion=max0(obj@mcPre-mk$cost),margin=max0(antitrust::calcMargins(obj)-mk$m),
       analytic_response=max0(ra-mk$reaction),implicit_response=max0(ri-mk$reaction),
       fully_solved_response_fd=max0(rf-ra),follower_actions=max0(sf$choices-mk$z),
       generic_equilibrium=max0(pg-mk$p),fully_solved_leader_gradient=lg)
  }
tab<-do.call(rbind,rows)
write.csv(tab,file.path(out_dir,"deterministic_validation.csv"),row.names=FALSE)
print(aggregate(.~demand+conduct+input,tab[setdiff(names(tab),"case")],max),digits=5)
stopifnot(max(tab$baseline_price)<2e-7,max(tab$baseline_quantity)<2e-7,
  max(tab$cost_inversion)<2e-9,max(tab$margin)<2e-8,
  max(tab$analytic_response)<2e-8,max(tab$implicit_response)<2e-7,
  max(tab$fully_solved_response_fd)<2e-7,max(tab$follower_actions)<2e-7,
  max(tab$generic_equilibrium)<2e-6,max(tab$fully_solved_leader_gradient)<2e-7)
cat("PASS: independent package equilibrium, cost, derivative and generic-route gates.\n")
