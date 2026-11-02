# Independent check of the exact new CES Bayes -> coordination normalization.
# Legacy CES applies its bounded logistic a second time to revenue outside
# share; this script quantifies that difference without changing legacy code.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else dirname(sub("^--file=", "",
  commandArgs()[grepl("^--file=", commandArgs())][1]))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(coordination))
suppressPackageStartupMessages(library(antitrust))
eps <- 1e-6
bounded <- function(x) eps+(1-2*eps)*plogis(x)
qc <- c(.23,.17,.32,.28)
owner <- c("A","A","B","C")
rows <- list(); i <- 0L
for (s0 in c(.005,.05,.3,.8)) for (p0 in c(.2,1.7,20))
  for (p in list(c(.7,1.8,2.3,4),c(7,18,23,40))) {
    avg <- sum(p*qc); insideQ <- 1234.5
    q <- insideQ*qc
    q0 <- insideQ*s0/(1-s0)
    eIn <- sum(p*q); E <- eIn+p0*q0
    physical <- p0*q0/E
    lp <- qlogis(s0)+log(p0)-log(avg)
    exact <- plogis(lp)
    legacy <- bounded(lp)
    r <- (1-exact)*p*qc/avg
    i <- i+1L
    rows[[i]] <- data.frame(s0=s0,p0=p0,avg_price=avg,
      exact_R0=exact,physical_error=abs(exact-physical),
      legacy_R0=legacy,legacy_difference=legacy-exact,
      expenditure_error=abs(eIn/(1-exact)/E-1),
      quantity_error=max(abs((E*r/p)/q-1)))
  }
tab <- do.call(rbind,rows)
stopifnot(max(tab$physical_error)<1e-14,
          max(tab$expenditure_error)<1e-13,max(tab$quantity_error)<1e-13,
          max(abs(tab$legacy_difference))>5e-7)
# Native reconstruction at the same unconditional revenue shares. All margin
# values are recovered from the audited FOC by a fixed-gamma calibration.
native <- list(); i <- 0L
for (cn in c("bertrand","cournot")) for (g in c(1.5,3.5)) {
  p <- c(1.7,2.1,2.9,2.5); p0 <- 1.7; s0 <- .13
  avg <- sum(p*qc); insideQ <- 1234.5
  r0 <- plogis(qlogis(s0)+log(p0)-log(avg))
  r <- (1-r0)*p*qc/avg
  obj <- coordination::stackelberg(prices=p,shares=r,margins=rep(NA_real_,4),
    ownerPre=owner,leadersPre="A",demand="ces",conduct=cn,gamma=g,
    insideSize=insideQ*avg,priceOutside=p0,
    labels=paste0("j",seq_along(p)),control.equ=list(implicitCheck=FALSE))
  i <- i+1L
  native[[i]] <- data.frame(conduct=cn,gamma=g,
    outside_error=abs((1-sum(obj@shares))-r0),
    expenditure_error=abs(obj@mktSize/(insideQ*avg/(1-r0))-1),
    quantity_error=max(abs(antitrust::calcQuantities(obj)/(insideQ*qc)-1)),
    baseline_price_relative=obj@diagnostics$baselineReproductionRelative)
}
native <- do.call(rbind,native)
stopifnot(max(native$outside_error)<1e-13,
          max(native$expenditure_error)<1e-13,max(native$quantity_error)<1e-12,
          max(native$baseline_price_relative)<2e-7)
write.csv(tab,file.path(out,"ces_outside_conversion_validation.csv"),row.names=FALSE)
write.csv(native,file.path(out,"ces_outside_native_validation.csv"),row.names=FALSE)
print(tab,digits=5);print(native,digits=5)
cat("PASS: 24 physical conversion cases and four native CES reconstructions.\n")
cat("Legacy revenue bounding is retained for old structures only; new games use exact physical R0.\n")
