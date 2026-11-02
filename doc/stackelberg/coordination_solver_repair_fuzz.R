library(antitrust)
library(coordination)
set.seed(91731)
ans <- list(); fixtures <- list()
add <- function(...) ans[[length(ans) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
for (demand in c("logit", "ces")) for (conduct in c("bertrand", "cournot")) for (z in seq_len(30)) {
  n <- sample(4:7, 1)
  ## Keep a positive outside share while exploring concentrated and tiny products.
  raw <- rexp(n, rate = runif(1, .5, 4))
  shares <- raw / sum(raw) * runif(1, .2, .94)
  if (z %% 10 == 0) shares[1:2] <- c(.48, .39) * runif(1, .8, .98)
  if (z %% 7 == 0) shares[n] <- 1e-8
  if (sum(shares) >= .995) shares <- shares / sum(shares) * .97
  prices <- exp(runif(n, log(.25), log(18)))
  firms <- paste0("f", sample(seq_len(max(3, ceiling(n / 2))), n, replace = TRUE))
  firmlev <- unique(firms)[sample(seq_along(unique(firms)), sample(0:length(unique(firms)), 1))]
  if (demand == "logit") {
    fit <- try(stackelberg(prices, shares, rep(NA_real_, n), ownerPre = firms,
                           leadersPre = firmlev, alpha = runif(1, .15, 4),
                           output = TRUE, conduct = conduct, insideSize = runif(1, 2, 500),
                           priceOutside = runif(1, 0, 4),
                           control.equ = list(implicitCheck = FALSE)), silent = TRUE)
  } else {
    fit <- try(stackelberg(prices, shares, rep(NA_real_, n), ownerPre = firms,
                           leadersPre = firmlev, gamma = exp(runif(1, log(.1), log(15))) + 1,
                           demand = "ces", conduct = conduct, insideSize = runif(1, 2, 500),
                           priceOutside = exp(runif(1, log(.1), log(8))),
                           control.equ = list(implicitCheck = FALSE)), silent = TRUE)
  }
  fixtures[[paste(demand,conduct,z)]] <- fit
  if (inherits(fit, "try-error")) {
    add(demand = demand, conduct = conduct, iter = z, pass = FALSE, type = "constructor", detail = as.character(fit))
    next
  }
  rr <- try(stackelberg_residuals(fit)$maxNormalized, silent = TRUE)
  finite <- is.finite(rr) && rr < 2e-10 && all(is.finite(fit@pricePre))
  add(demand = demand, conduct = conduct, iter = z, pass = finite, type = "equilibrium", detail = as.character(rr))
  ## Off-baseline follower solves at moderate deviations test exposed choice maps.
  li <- which(fit@firmOwnerPre %in% fit@leadersPre); fi <- which(!(fit@firmOwnerPre %in% fit@leadersPre))
  if (length(li) && length(fi)) {
    a <- if (conduct == "bertrand") fit@pricePre[li] else if (demand == "ces")
      coordination:::.ces_demand_from_prices(fit, fit@pricePre, rep(TRUE, n))$quantities[li]
      else fit@mktSize * fit@shares[li]
    for (sc in c(.75, 1.25)) {
      ff <- try(stackelberg_followers(fit, a * sc), silent = TRUE)
      finiteff <- !inherits(ff, "try-error") && all(is.finite(ff$prices[fi])) && all(is.finite(ff$quantities[fi]))
      add(demand = demand, conduct = conduct, iter = z, pass = finiteff, type = paste0("follower-", sc), detail = if (finiteff) "ok" else as.character(ff))
    }
  }
}
tab <- do.call(rbind, ans)
write.csv(tab, "/tmp/stackelberg-solver-repair-fuzz.csv", row.names = FALSE)
print(table(tab$demand, tab$conduct, tab$type, tab$pass))
bad <- tab[!tab$pass, , drop = FALSE]
if(nrow(bad)) {
 print(bad)
 for(i in seq_len(nrow(bad))) {
  row <- bad[i,]; stopifnot(row$demand=="logit",row$conduct=="bertrand",grepl("^follower-",row$type))
  fit<-fixtures[[paste(row$demand,row$conduct,row$iter)]]
  st<-coordination:::.sk_state(fit,TRUE);li<-which(st$owner%in%st$leaders);fi<-which(!st$owner%in%st$leaders)
  p<-fit@pricePre;p[li]<-p[li]*as.numeric(sub("follower-","",row$type))
  foc<-function(a){tmp<-fit;pp<-p;pp[fi]<-a;tmp@pricePre<-pp;s<-calcShares(tmp,TRUE,revenue=FALSE);mu<-pp-st$costs;H<-tapply(mu*s,st$owner,sum);1+fit@slopes$alpha*(mu[fi]-H[st$owner[fi]])}
  sol<-nleqslv::nleqslv(fit@pricePre[fi],foc,method="Newton",control=list(ftol=1e-12,xtol=1e-14,maxit=1000))
  stopifnot(max(abs(sol$fvec))<2e-10,any(sol$x<0))
  cat("CLASSIFIED unsupported positive-price-interior deviation: unconstrained followerprices ",sol$x," normalized FOC ",max(abs(sol$fvec)),"\n")
 }
}
cat("All 120 baseline constructions passed the unchanged 2e-10 normalized equilibrium gate.\n")
cat("Randomized numerical repair check complete; original unsupported deviation remains recorded in CSV.\n")
