# Compare timing endpoints directly with unchanged antitrust classes.
# Demand parameters and costs are fixed; this tests conduct rather than the
# different minimum-distance calibration criteria of the two packages.
args <- commandArgs(TRUE)
path <- if(length(args)) args[1] else "/tmp/stackelberg-coordination"
out <- if(length(args)>1) args[2] else "doc/stackelberg"
demands <- if(length(args)>2) args[3] else c("logit","ces")
pkgload::load_all(path, quiet=TRUE)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
oracle <- new.env(parent=globalenv())
for(e in parse(file.path(dirname(script),"math_check.R")))
  if(is.call(e) && identical(e[[1]],as.name("<-")) &&
     as.character(e[[2]]) %in% c("cases","make_market")) eval(e,oracle)
rows <- list(); i <- 0L
for(demand in demands) for(conduct in c("bertrand","cournot"))
 for(input in if(demand=="logit") c(FALSE,TRUE) else FALSE)
  for(layout in c("single","multi_leader","multi_follower"))
   for(endpoint in c("all","none")) {
    cc <- oracle$cases[[layout]]
    cc$leaders <- if(endpoint=="all") unique(cc$owner) else character()
    mk <- oracle$make_market(cc,demand,conduct,input)
    obj <- coordination::stackelberg(mk$p,cc$r,mk$m,ownerPre=mk$owner,
      leadersPre=mk$leaders,demand=demand,conduct=conduct,output=!input,
      insideSize=1000*sum(cc$r),priceOutside=1.7,
      alpha=if(demand=="logit") 2.3 else NULL,
      gamma=if(demand=="ces") 4 else NULL)
    ordinary_class <- if(demand=="logit") {
      if(conduct=="bertrand") "Logit" else "LogitCournot"
    } else if(conduct=="bertrand") "CES" else "CESCournot"
    ordinary <- methods::new(ordinary_class)
    for(nm in methods::slotNames(ordinary))
      methods::slot(ordinary,nm) <- methods::slot(obj,nm)
    ordinary@priceStart <- obj@prices*1.02
    suppressWarnings(pp <- antitrust::calcPrices(ordinary))
    mm <- antitrust::calcMargins(ordinary)
    # Recovered baseline primitives are identical at both timing endpoints.
    mc <- antitrust::calcMC(ordinary)
    i <- i+1L
    rows[[i]] <- data.frame(demand,conduct,input,layout,endpoint,
      margin=max(abs(mm-antitrust::calcMargins(obj))),
      cost=max(abs(mc-obj@mcPre)),
      price=max(abs(pp-obj@pricePre)))
   }
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,"ordinary_nesting_validation.csv"),row.names=FALSE)
print(aggregate(cbind(margin,cost,price)~demand+conduct+input,tab,max),digits=5)
stopifnot(max(tab$margin)<2e-9,max(tab$cost)<2e-9,max(tab$price)<2e-5)
cat("PASS: both timing endpoints agree with ordinary antitrust conduct.\n")
