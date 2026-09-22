#!/usr/bin/env Rscript
# Native Stan AD validation of the collapsed market-intercept likelihood.
# This script extracts the exact production helper definitions, compiles a
# small fixed-parameter harness, and compares normalized state and market
# mixture densities with a dense Gaussian oracle. It does not run HMC.

options(stringsAsFactors = FALSE, mc.cores = 1L)
suppressPackageStartupMessages(library(rstan))
rstan_options(auto_write = FALSE, threads_per_chain = 1L)
args <- commandArgs(FALSE)
script_arg <- sub("^--file=", "", args[startsWith(args, "--file=")])
outdir <- if (length(script_arg)) dirname(normalizePath(script_arg[1])) else "doc"
stan_root <- file.path(dirname(outdir), "inst", "stan")

extract_stan_function <- function(path, name) {
  src <- paste(readLines(path, warn = FALSE), collapse = "\n")
  hit <- regexpr(paste0("\\breal\\s+", name, "\\s*\\("), src, perl = TRUE)
  if (hit[1] < 0L) stop("Cannot find Stan function ", name, " in ", path)
  start <- as.integer(hit[1])
  rel <- regexpr("\\{", substring(src, start), perl = TRUE)[1]
  if (rel < 0L) stop("Cannot find body for ", name)
  open <- start + rel - 1L
  ch <- strsplit(substring(src, open), "", fixed = TRUE)[[1]]
  depth <- 0L; line_comment <- FALSE; block_comment <- FALSE
  string <- FALSE; escaped <- FALSE; close <- NA_integer_; i <- 1L
  while (i <= length(ch)) {
    x <- ch[i]; y <- if (i < length(ch)) ch[i + 1L] else ""
    if (line_comment) {
      if (x == "\n") line_comment <- FALSE
    } else if (block_comment) {
      if (x == "*" && y == "/") { block_comment <- FALSE; i <- i + 1L }
    } else if (string) {
      if (!escaped && x == '"') string <- FALSE
      if (!escaped && x == "\\") escaped <- TRUE else escaped <- FALSE
    } else if (x == "/" && y == "/") { line_comment <- TRUE; i <- i + 1L
    } else if (x == "/" && y == "*") { block_comment <- TRUE; i <- i + 1L
    } else if (x == '"') string <- TRUE
    else if (x == "{") depth <- depth + 1L
    else if (x == "}") {
      depth <- depth - 1L
      if (depth == 0L) { close <- open + i - 1L; break }
    }
    i <- i + 1L
  }
  if (is.na(close)) stop("Unbalanced function body for ", name)
  substring(src, start, close)
}

helper_names <- c("centered_intercept_gaussian_from_stats",
                  "centered_intercept_gaussian_density",
                  "collapsed_joint_from_local_intercept_stats")
kernel_files <- file.path(stan_root, c("logit_supply.stan", "ces_supply.stan",
  "logit_stackelberg.stan", "ces_stackelberg.stan"))
kernel_text <- lapply(kernel_files, function(path)
  paste(vapply(helper_names, function(n) extract_stan_function(path, n), ""),
        collapse = "\n\n"))
names(kernel_text) <- basename(kernel_files)
hash_text <- function(x) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(x, f, useBytes = TRUE)
  unname(tools::md5sum(f))
}
helper_hash <- vapply(kernel_text, hash_text, "")
representatives <- !duplicated(helper_hash)

model_code_for <- function(helper_code, model_name) paste0(
  "functions {\n", helper_code, "\n}\n",
  "data {\n",
  " int<lower=1> N; int<lower=1> M; int<lower=1> K; int<lower=1> Q;\n",
  " array[M] int<lower=1> G_market;\n",
  " array[M] int<lower=1> obs_start; array[M] int<lower=1> obs_n;\n",
  " array[N] int<lower=1,upper=N> market_obs_order;\n",
  " array[N] int<lower=1,upper=N> market_obs_position;\n",
  " array[N] int<lower=1> quality_local;\n",
  " array[M] int<lower=1> margin_start; array[M] int<lower=0> margin_n;\n",
  " array[Q] int<lower=1,upper=Q> margin_order;\n",
  " array[Q] int<lower=1,upper=N> margin_idx; matrix[M,K] log_k_prior;\n",
  "}\nparameters {\n",
  " vector[N] share_resid; matrix[K,Q] margin_resid;\n",
  " real mu_intercept;\n",
  " real log_sigma_share; real log_sigma_margin;\n",
  " real<lower=0> sigma_quality; real<lower=0> sigma_intercept;\n",
  " real rho_raw;\n}\n",
  "transformed parameters {\n",
  " real sigma_share=exp(log_sigma_share);\n",
  " real sigma_margin=exp(log_sigma_margin);\n",
  " real rho=0.95*tanh(rho_raw);\n",
  " vector[N] share_resid_centered=share_resid-rep_vector(mu_intercept,N);\n",
  " matrix[M,K] log_joint_market_k; vector[M] log_lik_market;\n",
  " for (m in 1:M) { vector[K] score;\n",
  "  for (k in 1:K) { vector[margin_n[m]] local_margin_resid;\n",
  "   for (q in 1:margin_n[m])\n",
  "    local_margin_resid[q]=margin_resid[k,margin_start[m]+q-1];\n",
  "   log_joint_market_k[m,k]=collapsed_joint_from_local_intercept_stats(\n",
  "    share_resid_centered,local_margin_resid,margin_order,margin_idx,margin_start[m],\n",
  "    market_obs_order,market_obs_position,obs_start[m],quality_local,\n",
  "    obs_n[m],margin_n[m],G_market[m],sigma_share,sigma_margin,rho,\n",
  "    sigma_quality,sigma_intercept);\n",
  "   score[k]=log_k_prior[m,k]+log_joint_market_k[m,k];\n",
  "  } log_lik_market[m]=log_sum_exp(score);\n",
  " }\n}\nmodel {\n",
  " target += sum(log_lik_market);\n",
  " target += normal_lpdf(mu_intercept | 0, 1);\n",
  " target += normal_lpdf(sigma_intercept | 0, 0.5)\n",
  "   - normal_lccdf(0 | 0, 0.5);\n",
  " target += normal_lpdf(sigma_quality | 1, 1)\n",
  "   - normal_lccdf(0 | 1, 1);\n",
  " target += normal_lpdf(log_sigma_share | log(0.5), 0.5);\n",
  " target += normal_lpdf(log_sigma_margin | log(0.75), 1);\n",
  " target += log1m(square(rho))\n",
  "   + log1m(square(tanh(rho_raw)))\n",
  "   - (lbeta(0.5, 2.0) + beta_lcdf(square(0.95) | 0.5, 2.0)\n",
  "      - log(0.95));\n",
  "}\n")

# Three markets have 4/3/5 observations and 2/3/3 quality groups. The middle
# market has no observed margins; the others have incomplete margin coverage.
N <- 12L; M <- 3L; K <- 3L; Q <- 5L
dat <- list(N=N,M=M,K=K,Q=Q,G_market=c(2L,3L,3L),
  obs_start=c(1L,5L,8L),obs_n=c(4L,3L,5L),
  market_obs_order=seq_len(N),market_obs_position=c(1:4,1:3,1:5),
  quality_local=c(1L,1L,2L,2L,1L,2L,3L,1L,2L,1L,3L,2L),
  margin_start=c(1L,3L,3L),margin_n=c(2L,0L,3L),
  margin_order=seq_len(Q),margin_idx=c(2L,4L,8L,10L,12L),
  log_k_prior=rbind(log(c(.52,.31,.17)),log(c(.22,.50,.28)),
                    log(c(.39,.26,.35))))

log_gaussian <- function(y,V) {
  if (!length(y)) return(0)
  U <- chol(V); z <- forwardsolve(t(U),y)
  -.5*(length(y)*log(2*pi)+2*sum(log(diag(U)))+sum(z*z))
}
dense_state <- function(sr,mr,m,ss,sm,rho,tau,sb) {
  obs <- seq.int(dat$obs_start[m],length.out=dat$obs_n[m])
  g <- dat$quality_local[obs]; G <- dat$G_market[m]
  Z <- matrix(0,length(obs),G); Z[cbind(seq_along(obs),g)] <- 1
  P <- diag(G)-matrix(1/G,G,G)
  qpos <- if(dat$margin_n[m]) seq.int(dat$margin_start[m],
    length.out=dat$margin_n[m]) else integer()
  pair <- dat$margin_idx[qpos]; local <- match(pair,obs); nq <- length(qpos)
  C <- matrix(0,length(obs),nq)
  if(nq) for(j in seq_len(nq)) C[local[j],j] <- rho*ss*sm
  Vss <- tau^2*Z%*%P%*%t(Z)+sb^2*matrix(1,length(obs),length(obs))+
    diag(ss^2,length(obs))
  V <- rbind(cbind(Vss,C),cbind(t(C),diag(sm^2,nq)))
  log_gaussian(c(sr[obs],mr[qpos]),V)
}
oracle <- function(pars) {
  sr <- pars$share_resid - pars$mu_intercept
  mr <- pars$margin_resid
  ss <- exp(pars$log_sigma_share)
  sm <- exp(pars$log_sigma_margin)
  sb <- pars$sigma_intercept
  tau <- pars$sigma_quality
  rho <- .95*tanh(pars$rho_raw)
  ll <- matrix(NA_real_,M,K)
  for(m in seq_len(M)) for(k in seq_len(K))
    ll[m,k] <- dense_state(sr,mr[k,],m,ss,sm,rho,tau,sb)
  mix <- sum(vapply(seq_len(M),function(m) {
    x <- dat$log_k_prior[m,]+ll[m,]; a <- max(x); a+log(sum(exp(x-a)))
  },numeric(1)))
  hyperprior <- dnorm(pars$mu_intercept,0,1,log=TRUE)+
    dnorm(sb,0,.5,log=TRUE)-pnorm(0,0,.5,lower.tail=FALSE,log.p=TRUE)+
    dnorm(tau,1,1,log=TRUE)-pnorm(0,1,1,lower.tail=FALSE,log.p=TRUE)+
    dnorm(pars$log_sigma_share,log(.5),.5,log=TRUE)+
    dnorm(pars$log_sigma_margin,log(.75),1,log=TRUE)+
    log1p(-rho^2)+log1p(-tanh(pars$rho_raw)^2)-
    (lbeta(.5,2)+pbeta(.95^2,.5,2,log.p=TRUE)-log(.95))
  list(target=mix+hyperprior,state=ll,scales=c(ss,sm,tau,sb),rho=rho,
       hyperprior=hyperprior)
}

set.seed(49271L)
base_share <- rnorm(N,0,.3); base_margin <- matrix(rnorm(K*Q,0,.2),K,Q)
cases <- list(
  ordinary=c(ss=.31,sm=.43,tau=.48,sb=.37,rho=.63),
  low_share_high_rho=c(ss=.012,sm=.34,tau=.52,sb=.41,rho=.939),
  low_both_scales=c(ss=.018,sm=.29,tau=.009,sb=.004,rho=-.91),
  low_intercept_scale=c(ss=.27,sm=.39,tau=.44,sb=.0015,rho=.72),
  rho_zero=c(ss=.25,sm=.36,tau=.31,sb=.22,rho=0))

target_rows <- list(); gradient_rows <- list(); hash_rows <- list()
for (src in names(kernel_text)[representatives]) {
  code <- model_code_for(kernel_text[[src]],
                         paste0("MarketInterceptCheck",gsub("\\W","",src)))
  smod <- stan_model(model_code=code,
    model_name=paste0("MarketInterceptCheck",gsub("\\W","",src)),
    auto_write=FALSE,verbose=FALSE)
  hash_rows[[src]] <- data.frame(source=src,helper_md5=helper_hash[[src]])
  for (case in names(cases)) {
    z <- cases[[case]]
    init <- list(share_resid=base_share,margin_resid=base_margin,
      mu_intercept=.13,
      log_sigma_share=log(z[["ss"]]),log_sigma_margin=log(z[["sm"]]),
      sigma_quality=z[["tau"]],sigma_intercept=z[["sb"]],
      rho_raw=atanh(z[["rho"]]/.95))
    fit <- sampling(smod,data=dat,chains=1L,iter=1L,warmup=0L,
      algorithm="Fixed_param",init=list(init),refresh=0L,
      pars=c("log_joint_market_k","log_lik_market"))
    u <- unconstrain_pars(fit,init)
    native_lp <- as.numeric(log_prob(fit,u,adjust_transform=FALSE))
    ad <- as.numeric(grad_log_prob(fit,u,adjust_transform=FALSE))
    draws <- extract(fit,pars=c("log_joint_market_k","log_lik_market"),
                     permuted=TRUE)
    native_state <- matrix(draws$log_joint_market_k[1,,,drop=FALSE],M,K)
    native_mix <- as.numeric(draws$log_lik_market[1,])
    # Decode through Stan rather than assume its unconstrained ordering for
    # matrix parameters; this also checks the fixed-state parameter roundtrip.
    constrained <- constrain_pars(fit,u)
    ref <- oracle(constrained)
    ref_mix <- vapply(seq_len(M),function(m) {
      x <- dat$log_k_prior[m,]+ref$state[m,]; a <- max(x); a+log(sum(exp(x-a)))
    },numeric(1))
    target_rows[[paste(src,case)]] <- data.frame(source=src,case=case,
      N=N,markets=M,k_states=K,margin_obs=Q,no_margin_markets=sum(dat$margin_n==0),
      sigma_share=ref$scales[1],sigma_margin=ref$scales[2],
      sigma_quality=ref$scales[3],sigma_intercept=ref$scales[4],rho=ref$rho,
      stan_log_target=native_lp,dense_oracle_log_target=ref$target,
      abs_target_error=abs(native_lp-ref$target),
      max_state_error=max(abs(native_state-ref$state)),
      max_market_mixture_error=max(abs(native_mix-ref_mix)),
      pass=max(abs(native_lp-ref$target),abs(native_state-ref$state),
               abs(native_mix-ref_mix))<2e-8)
    # Independent finite differences use the dense R oracle after decoding
    # each unconstrained state through Stan, so matrix ordering cannot be a
    # source of agreement or error.
    fd <- numeric(length(u))
    for(j in seq_along(u)) {
      h <- 1e-5*max(1,abs(u[j])); up <- um <- u
      up[j] <- up[j]+h; um[j] <- um[j]-h
      fd[j] <- (oracle(constrain_pars(fit,up))$target-
        oracle(constrain_pars(fit,um))$target)/(2*h)
    }
    pname <- names(u)
    if (is.null(pname)) pname <- paste0("u[",seq_along(u),"]")
    gradient_rows[[paste(src,case)]] <- data.frame(source=src,case=case,
      parameter=pname,native_ad=ad,oracle_central_difference=fd,
      abs_error=abs(ad-fd),relative_error=abs(ad-fd)/pmax(1,abs(ad),abs(fd)),
      pass=abs(ad-fd)<3e-5+3e-5*abs(fd))
  }
}

targets <- do.call(rbind,target_rows)
grads <- do.call(rbind,gradient_rows)
write.csv(targets,file.path(outdir,
  "market_intercept_collapse_native_validation.csv"),row.names=FALSE)
write.csv(grads,file.path(outdir,
  "market_intercept_collapse_native_validation_gradients.csv"),row.names=FALSE)
write.csv(do.call(rbind,hash_rows),file.path(outdir,
  "market_intercept_collapse_native_validation_helper_hashes.csv"),row.names=FALSE)
cat("Compiled extracted production-helper harnesses; no HMC was run.\n")
cat("Native likelihood target checks:",sum(targets$pass),"/",nrow(targets),
    "; max target error",format(max(targets$abs_target_error),digits=6),
    "; max state error",format(max(targets$max_state_error),digits=6),
    "; max k-mixture error",format(max(targets$max_market_mixture_error),digits=6),"\n")
cat("Native autodiff checks:",sum(grads$pass),"/",nrow(grads),
    "; max absolute gradient error",format(max(grads$abs_error),digits=6),
    "; max scaled error",format(max(grads$relative_error),digits=6),"\n")
stopifnot(all(targets$pass),all(grads$pass))
