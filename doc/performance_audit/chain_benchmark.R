options(warn = 1)
suppressPackageStartupMessages(library(antitrustBayes))
args <- commandArgs(TRUE)
fixture <- args[1L]; cores <- as.integer(args[2L]); seed <- as.integer(args[3L]); out <- args[4L]
f <- readRDS(fixture)
d <- antitrustBayes:::.bayes_bridge_canonicalize_data(f@stan_data)
if (is.null(d$market_s0_parameterization)) d$market_s0_parameterization <- 1L
kernel <- f@stanfit@model_name
model <- antitrustBayes:::stanmodels[[kernel]]
starts <- f@stanfit@inits
starts <- lapply(starts, function(x) {x$c_s0 <- numeric(0); x$c_s0_raw <- numeric(0); x})
settings <- list(chains = 4L, iter = 1000L, warmup = 500L, cores = cores, seed = seed,
                 init = starts, refresh = 0L,
                 control = list(adapt_delta = .95, max_treedepth = 13L, metric = 'diag_e'))
manifest <- list(fixture_hash = digest::digest(file = fixture, algo='sha256'),
                 data_hash = digest::digest(d, algo='sha256'), starts_hash=digest::digest(starts,algo='sha256'),
                 kernel = kernel, model_hash=digest::digest(model@model_code,algo='sha256'), settings=settings)
saveRDS(manifest,paste0(out,'.manifest.rds'))
warnings <- character()
time <- system.time({fit <- tryCatch(withCallingHandlers(
  do.call(rstan::sampling, c(list(object=model, data=d),settings)),
  warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)})
if (inherits(fit,'condition') || !methods::is(fit,'stanfit') || fit@mode != 0L || !length(fit@sim$samples)) {
  msg <- if(inherits(fit,'condition')) conditionMessage(fit) else 'No populated posterior stanfit'
  saveRDS(fit,paste0(out,'.failed.rds')); stop(msg)
}
saveRDS(fit,paste0(out,'.stanfit.rds'))
post <- rstan::get_sampler_params(fit,inc_warmup=FALSE)
allsp <- rstan::get_sampler_params(fit,inc_warmup=TRUE)
bfmi <- vapply(post,function(x) mean(diff(x[,'energy__'])^2)/var(x[,'energy__']),numeric(1))
keep <- intersect(c('alpha_structural_market','s0_prob','sigma_b_strat','sigma_share','sigma_margin'),fit@sim$pars_oi)
a <- rstan::extract(fit,pars=keep,permuted=FALSE,inc_warmup=FALSE)
sm <- posterior::summarise_draws(posterior::as_draws_array(a))
clock <- rstan::get_elapsed_time(fit)
row <- data.frame(fixture=basename(fixture),kernel=kernel,cores=cores,seed=seed,
 wall_seconds=unname(time['elapsed']),native_warmup_seconds=sum(clock[,1]),native_sampling_seconds=sum(clock[,2]),
 divergences=sum(vapply(post,function(x)sum(x[,'divergent__']),numeric(1))),
 treedepth_hits=sum(vapply(post,function(x)sum(x[,'treedepth__']>=13),numeric(1))),
 min_bfmi=min(bfmi),max_rhat=max(sm$rhat),min_bulk_ess=min(sm$ess_bulk),min_tail_ess=min(sm$ess_tail),
 structural_ess_per_second=min(sm$ess_bulk)/unname(time['elapsed']),
 leapfrog_evaluations=sum(vapply(allsp,function(x)sum(x[,'n_leapfrog__']),numeric(1))),
 posterior_draw_hash=digest::digest(rstan::extract(fit,permuted=FALSE,inc_warmup=FALSE),algo='sha256'),
 warning_count=length(warnings),warnings=paste(unique(warnings),collapse=' | '))
write.csv(row,paste0(out,'.csv'),row.names=FALSE)
print(row)
