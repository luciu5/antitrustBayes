#!/usr/bin/env Rscript
# Independent feasibility oracle only. This does not modify the Stan model.
# Integrate a Normal market intercept as well as centered Normal group quality.
options(stringsAsFactors = FALSE)
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[startsWith(args, "--file=")])
outdir <- if (length(script)) dirname(normalizePath(script[1])) else "doc"

log_gaussian <- function(z, V) {
  U <- chol(V)
  u <- forwardsolve(t(U), z)
  -.5 * (length(z) * log(2*pi) + 2 * sum(log(diag(U))) + sum(u*u))
}

# z excludes the prior mean of the intercept. Z is one-hot group incidence.
# Cov(z) = diag(v) + tau^2 Z P Z' + market_sd^2 11'.
# The positive-sum form below avoids subtracting large quadratic terms.
collapse_market <- function(z, v, group, tau, market_sd) {
  G <- max(group)
  if (G == 1L) tau <- 0 # centered quality is deterministically zero
  a <- cstat <- numeric(G)
  for (i in seq_along(z)) {
    a[group[i]] <- a[group[i]] + 1/v[i]
    cstat[group[i]] <- cstat[group[i]] + z[i]/v[i]
  }
  means <- cstat/a
  within <- sum((z-means[group])^2/v)
  w <- 1/(1+tau^2*a)
  aw <- a*w
  A <- sum(aw)
  cbar <- sum(aw*means)/A
  dc <- mean(w)
  den <- dc + market_sd^2*A
  logdet <- sum(log(v)) + sum(log1p(tau^2*a)) + log(den)
  quadratic <- within + sum(aw*(means-cbar)^2) + A*cbar^2/den
  list(log_density=-.5*(length(z)*log(2*pi)+logdet+quadratic),
       intercept_mean=market_sd^2*A*cbar/den,
       intercept_variance=market_sd^2*dc/den)
}

joint_structured <- function(rs, rm, group, paired, tau, ss, sm, rho, sb) {
  z <- rs
  v <- rep(ss^2, length(rs))
  z[paired] <- z[paired] - rho*ss/sm*rm
  v[paired] <- ss^2*(1-rho^2)
  ans <- collapse_market(z,v,group,tau,sb)
  ans$log_density <- ans$log_density + sum(dnorm(rm,0,sm,log=TRUE))
  ans
}

joint_covariance <- function(group, paired, tau, ss, sm, rho, sb) {
  N <- length(group); G <- max(group); Q <- length(paired)
  Z <- diag(G)[group,,drop=FALSE]
  P <- diag(G)-matrix(1/G,G,G)
  E <- diag(N)[,paired,drop=FALSE]
  rbind(cbind(tau^2*Z%*%P%*%t(Z) + diag(ss^2,N) + matrix(sb^2,N,N),
              rho*ss*sm*E),
        cbind(rho*ss*sm*t(E),diag(sm^2,Q)))
}

fixtures <- list(
  singleton=list(g=rep(1L,3),paired=c(1L,3L)),
  cross_section=list(g=1:6,paired=1:6),
  multiproduct=list(g=c(1L,1L,2L,3L,3L,3L),paired=c(1L,3L,5L)),
  unbalanced_panel=list(g=c(1L,1L,2L,3L,1L,2L,2L,3L,4L),
                        paired=c(1L,3L,4L,5L,8L,9L))
)
grid <- expand.grid(tau=c(0,.01,1.5),ss=c(.001,.25,.8),
                    sm=c(.002,.4,1.1),rho=c(-.94,0,.94),sb=c(0,.03,.8))
set.seed(20260920L)
density_rows <- list(); integral_rows <- list(); gradient_rows <- list()
for (name in names(fixtures)) {
  fx <- fixtures[[name]]; N <- length(fx$g); G <- max(fx$g)
  latent <- rnorm(G); latent <- latent-mean(latent)
  e <- rnorm(N); em <- rnorm(length(fx$paired)); market_z <- rnorm(1)
  for (j in seq_len(nrow(grid))) {
    p <- grid[j,]
    rs <- p$tau*latent[fx$g] + p$ss*e + p$sb*market_z
    rm <- p$sm*(p$rho*e[fx$paired]+sqrt(1-p$rho^2)*em)
    args <- c(list(rs=rs,rm=rm,group=fx$g,paired=fx$paired),as.list(p))
    st <- do.call(joint_structured,args)
    V <- do.call(joint_covariance,c(list(group=fx$g,paired=fx$paired),as.list(p)))
    dense <- log_gaussian(c(rs,rm),V)
    U <- chol(V)
    cv <- c(rep(p$sb^2,N),rep(0,length(rm)))
    invcv <- backsolve(U,forwardsolve(t(U),cv))
    mean_dense <- sum(invcv*c(rs,rm))
    var_dense <- p$sb^2-sum(cv*invcv)
    err <- abs(st$log_density-dense)
    merr <- abs(st$intercept_mean-mean_dense)
    verr <- abs(st$intercept_variance-var_dense)
    density_rows[[length(density_rows)+1L]] <- data.frame(
      fixture=name,p,density_error=err,posterior_mean_error=merr,
      posterior_variance_error=verr,tolerance=1e-6,
      passed=max(err,merr,verr)<1e-6)
  }
  # Integrate the original normalized likelihood times intercept prior directly.
  # Change to the conditional intercept's standard-normal coordinate so narrow
  # peaks are not missed by quadrature. Dense C is reference-only.
  for (ss in c(.02,.4)) for (rho in c(0,.8)) {
    tau <- .7; sm <- .5; sb <- .8
    rs <- tau*latent[fx$g]+ss*e+.35
    rm <- sm*(rho*e[fx$paired]+sqrt(1-rho^2)*em)
    st <- joint_structured(rs,rm,fx$g,fx$paired,tau,ss,sm,rho,sb)
    C <- joint_covariance(fx$g,fx$paired,tau,ss,sm,rho,0)
    post_sd <- sqrt(st$intercept_variance)
    integrand <- Vectorize(function(q) {
      b <- st$intercept_mean+post_sd*q
      exp(log_gaussian(c(rs-b,rm),C)+dnorm(b,0,sb,log=TRUE)+
            log(post_sd)-st$log_density)
    })
    val <- integrate(integrand,-10,10,rel.tol=1e-9)$value
    integral_rows[[length(integral_rows)+1L]] <- data.frame(
      fixture=name,ss=ss,rho=rho,normalized_integral=val,
      error=abs(val-1),passed=abs(val-1)<1e-7)
  }
  for (ss in c(.005,.25)) for (rho in c(-.8,.8)) {
    tau <- .7; sm <- .5; sb <- .8
    rs <- tau*latent[fx$g]+ss*e+.35
    rm <- sm*(rho*e[fx$paired]+sqrt(1-rho^2)*em)
    x <- c(log(tau),log(ss),log(sm),atanh(rho/.95),log(sb))
    eval_target <- function(x,dense=FALSE) {
      p <- list(tau=exp(x[1]),ss=exp(x[2]),sm=exp(x[3]),
                rho=.95*tanh(x[4]),sb=exp(x[5]))
      if (dense) return(log_gaussian(c(rs,rm),do.call(joint_covariance,
        c(list(group=fx$g,paired=fx$paired),p))))
      do.call(joint_structured,c(list(rs=rs,rm=rm,group=fx$g,
                                     paired=fx$paired),p))$log_density
    }
    a <- numDeriv::grad(eval_target,x)
    b <- numDeriv::grad(function(x) eval_target(x,TRUE),x)
    error <- max(abs(a-b)/(1+abs(a)))
    gradient_rows[[length(gradient_rows)+1L]] <- data.frame(
      fixture=name,ss=ss,rho=rho,max_relative_gradient_error=error,
      tolerance=2e-4,passed=error<2e-4)
  }
}
tables <- list(density=do.call(rbind,density_rows),
               integral=do.call(rbind,integral_rows),
               gradient=do.call(rbind,gradient_rows))
for (name in names(tables)) {
  x <- tables[[name]]
  write.csv(x,file.path(outdir,paste0("market_intercept_collapse_",name,"_checks.csv")),
            row.names=FALSE)
  cat(name, sum(x$passed),"/",nrow(x),"passed\n")
}
cat("maximum density error:",max(tables$density$density_error),"\n")
cat("maximum gradient error:",max(tables$gradient$max_relative_gradient_error),"\n")
stopifnot(all(vapply(tables,function(x) all(x$passed),logical(1))))
