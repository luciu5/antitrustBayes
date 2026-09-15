# Independent all-margin calibration and FOC-inverted-cost checks.
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
rows <- list(); ii <- 0L
for(demand in demands) for(conduct in c("bertrand","cournot"))
 for(input in if(demand=="logit") c(FALSE,TRUE) else FALSE)
  for(layout in c("single","multi_leader","multi_follower"))
   for(moment_case in c("noisy_partial","single_exact")) {
    cc <- oracle$cases[[layout]]; mk <- oracle$make_market(cc,demand,conduct,input)
    n <- length(mk$p); weights <- seq_len(n)/n
    observed <- mk$m * rep(c(1.025,.985,1.01,.99,1.015),length.out=n)
    if(moment_case=="noisy_partial") observed[3] <- NA_real_ else {
      observed <- rep(NA_real_,n); observed[2] <- mk$m[2]
    }
    usable <- is.finite(observed)
    if(demand=="logit") {
      x <- mk$m*2.3
      lambda <- sum(weights[usable]*x[usable]*observed[usable])/
        sum(weights[usable]*x[usable]^2)
      expected_param <- 1/lambda
      expected_margin <- x*lambda
      objective_gap <- 0
    } else {
      obj <- function(gamma) {
        trial <- cc; trial$gamma <- gamma
        pred <- oracle$make_market(trial,"ces",conduct)$m
        sum(weights[usable]*(log(1/observed[usable])-log(1/pred[usable]))^2)
      }
      opt <- optimize(function(z) obj(1+exp(z)),c(log(1e-7),log(99)),tol=1e-11)
      expected_param <- 1+exp(opt$minimum)
      trial <- cc; trial$gamma <- expected_param
      expected_margin <- oracle$make_market(trial,"ces",conduct)$m
      objective_gap <- opt$objective
    }
    o <- coordination::stackelberg(mk$p,cc$r,observed,ownerPre=mk$owner,
      leadersPre=mk$leaders,demand=demand,conduct=conduct,output=!input,
      weights=weights,insideSize=1000*sum(cc$r),priceOutside=1.7)
    got_param <- if(demand=="logit") abs(o@slopes$alpha) else o@slopes$gamma
    if(demand=="ces") objective_gap <- abs(obj(got_param)-objective_gap)
    implied <- antitrust::calcMargins(o)
    expected_cost <- mk$p*(1-if(input)-expected_margin else expected_margin)
    perturbed <- o; perturbed@priceStart <- o@prices*1.02
    reproduced <- antitrust::calcPrices(perturbed)
    residual <- coordination::stackelberg_residuals(o)$all
    scale <- if(conduct=="bertrand") antitrust::calcQuantities(o) else mk$p
    ii <- ii+1L
    rows[[ii]] <- data.frame(demand,conduct,input,layout,moment_case,
      parameter_error=abs(got_param-expected_param),objective_gap,
      implied_margin_error=max(abs(implied-expected_margin)),
      cost_error=max(abs(o@mcPre-expected_cost)),
      baseline_error=max(abs(reproduced-mk$p)),
      normalized_foc=max(abs(residual/scale)),
      noisy_cost_difference=max(abs(o@mcPre-mk$p*(1-if(input)-observed else observed)),na.rm=TRUE))
   }
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,"calibration_validation.csv"),row.names=FALSE)
print(aggregate(cbind(parameter_error,objective_gap,implied_margin_error,cost_error,
  baseline_error,normalized_foc)~demand+conduct+input,tab,max),digits=5)
stopifnot(max(tab$parameter_error)<3e-6,max(tab$objective_gap)<1e-10,
  max(tab$implied_margin_error)<2e-7,max(tab$cost_error)<5e-7,
  max(tab$baseline_error)<2e-7,max(tab$normalized_foc)<2e-8)
stopifnot(all(tab$noisy_cost_difference[tab$moment_case=="noisy_partial"]>1e-5))
cat("PASS: independent weighted calibration, missing/single margins, FOC costs and baseline reproduction.\n")
