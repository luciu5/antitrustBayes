#!/usr/bin/env Rscript
# Sourcing defines functions only. Run from a development-source session with
# compatible libraries/toolchains and explicit external cache options.
relative_validation_fixtures <- function() {
  source('doc/validation/nutpie_backend_validation.R', local = TRUE)
  source('doc/validation/single_market_intercept_validation.R', local = TRUE)
  base <- nutpie_validation_fixture('C')$arguments
  clean <- function(x) {
    x$market_intercept_parameterization <- NULL
    x$share_basis <- 'relative'
    x$quality_parameterization <- 'collapsed'
    x$use_rho <- TRUE
    x$use_hmt <- FALSE
    x
  }
  base <- clean(base)
  models <- c('MC','B','C','BF','CF','BS','CS')
  out <- stats::setNames(lapply(models, function(model) {
    x <- base
    x$supply_model <- if (model == 'MC') 'monopolistic' else
      if (model %in% c('C','CF','CS')) 'cournot' else 'bertrand'
    x$structure <- if (model %in% c('BF','CF')) 'fringe' else
      if (model %in% c('BS','CS')) 'stackelberg' else 'simultaneous'
    # Explicit rich structure determines the conduct-role flag. Supplying
    # use_cutoff=FALSE as well is an invalid contradictory public call.
    x$use_cutoff <- if (x$structure == 'simultaneous') FALSE else NULL
    x$partition_prior <- if (x$structure == 'simultaneous') 'cutoff' else 'ordered_k'
    x$k_support <- if (x$structure == 'simultaneous') NULL else 2L
    x$k_prior <- if (x$structure == 'simultaneous') NULL else 1
    x
  }), paste0('logit_', models))
  # Original conditional-k family support is fixed before observing results.
  for (k in c(1L,3L)) {
    x <- out$logit_CF; x$k_support <- k
    out[[paste0('logit_CF_k',k)]] <- x
  }
  for (model in c('C','CS')) {
    x <- out[[paste0('logit_',model)]]
    x$demand_model <- 'ces'; x$outside_price <- 1.2
    params <- bayes_params(gamma = 6, s0 = .2, supply_model = 3L,
      structure = x$structure, partition_prior = x$partition_prior,
      k_support = x$k_support, k_prior = x$k_prior)
    market <- bayes_market(x$prices, x$quantities, x$margins, group = x$group,
      outside_price = 1.2, insideSize = 1000, market_side = 'output')
    specified <- bayes_specify(params, market)
    model_object <- suppressWarnings(posterior_policy_models(specified, n_draws=1L))$models[[1L]][[1L]]
    if (methods::is(model_object,'AntitrustFit')) model_object <- model_object@model
    x$margins <- as.numeric(antitrust::calcMargins(model_object,preMerger=TRUE))
    stopifnot(all(is.finite(x$margins)),all(x$margins>0 & x$margins<1))
    out[[paste0('ces_',model)]] <- clean(x)
  }
  # Explicit M4T4: common q within each market, separate share contrasts per MT.
  design <- smi_design()
  fixture <- smi_fixture(data.frame(id='C_M4T4',conduct='C',markets=4L,periods=4L), design)
  out$logit_C_M4T4 <- clean(fixture$arguments)
  x <- out$ces_C; p <- x$prices; r <- x$quantities/sum(x$quantities)
  gamma <- 6; q <- log(r)-mean(log(r))+gamma*(log(p)-mean(log(p)))
  rows <- list()
  for (m in 1:4) for (tt in 1:4) {
    price <- p * exp(.02*sin(seq_along(p)+tt+m))
    v <- q-gamma*log(price); quantity <- exp(v-max(v)); quantity <- 1000*quantity/sum(quantity)
    market <- bayes_market(price,quantity,rep(.3,length(p)),group=x$group,
      outside_price=1.2,insideSize=1000,market_side='output')
    specified <- bayes_specify(bayes_params(gamma=gamma,s0=.2,supply_model=3L),market)
    object <- suppressWarnings(posterior_policy_models(specified,n_draws=1L))$models[[1L]][[1L]]
    if(methods::is(object,'AntitrustFit'))object<-object@model
    rows[[length(rows)+1L]] <- data.frame(prices=price,quantities=quantity,
      margins=as.numeric(antitrust::calcMargins(object,preMerger=TRUE)),
      group=x$group,market=paste0('m',m),year=paste0('t',tt))
  }
  data <- do.call(rbind,rows)
  x[names(data)] <- as.list(data); x$margin_observed <- rep(TRUE,nrow(data))
  x$s0_pooling <- 'market'; x$use_alpha_market_pool <- TRUE
  out$ces_C_M4T4 <- x
  out
}

relative_validation_run <- function(root, max_seconds=13200, run_policy=TRUE, resume=FALSE,
    include_ces_noiseless_stress=TRUE) {
  source('doc/validation/nutpie_backend_validation.R', local=TRUE)
  stopifnot(max_seconds>0,max_seconds<=14400)
  dir.create(root,recursive=TRUE,showWarnings=FALSE); root<-normalizePath(root)
  manifest_path<-file.path(root,'sampling-manifest.rds')
  prior_manifest<-if(file.exists(manifest_path))readRDS(manifest_path)else NULL
  if(!is.null(prior_manifest) && !isTRUE(resume))
    stop('Preserve existing validation; use explicit resume or select a fresh output directory.')
  preparation_warnings<-character()
  fixtures<-withCallingHandlers(relative_validation_fixtures(),warning=function(w){
    preparation_warnings<<-c(preparation_warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  controls<-list(chains=4L,iter=1000L,warmup=500L,cores=4L,refresh=0L,
    seed=452891L,control=list(adapt_delta=.99,max_treedepth=13L))
  manifest<-list(started=Sys.time(),deadline=Sys.time()+max_seconds,fixtures=fixtures,
    controls=controls,starting_commit=trimws(system2('git',c('rev-parse','HEAD'),stdout=TRUE)),
    kernel_hashes=vapply(list.files('inst/stan','[.]stan$',full.names=TRUE),
      digest::digest,character(1),file=TRUE,algo='sha256'),
    runtime=nutpie_validation_runtime(normalizePath('.')),status='running',
    preparation_warnings=preparation_warnings,attempts=list())
  if(!is.null(prior_manifest)) {
    if(!identical(prior_manifest$controls,controls) ||
       !identical(prior_manifest$kernel_hashes,manifest$kernel_hashes))
      stop('Cannot resume changed fixtures, sampler controls, or Stan source.')
    snapshot<-file.path(root,paste0('manifest-before-resume-',format(Sys.time(),'%Y%m%dT%H%M%S'),'.rds'))
    if(file.exists(snapshot))stop('Resume snapshot already exists.')
    saveRDS(prior_manifest,snapshot)
    changed<-names(fixtures)[!vapply(names(fixtures),function(id)
      identical(prior_manifest$fixtures[[id]],fixtures[[id]]),logical(1))]
    for(id in changed) {
      old<-prior_manifest$fixtures[[id]];corrected<-old
      corrected$use_cutoff<-NULL
      fitted<-vapply(c('cmdstanr','nutpie'),function(backend)
        file.exists(file.path(root,paste(id,backend,sep='_'),'fit.rds')),logical(1))
      if(any(fitted) || !isFALSE(old$use_cutoff) ||
         !old$structure %in% c('fringe','stackelberg') ||
         !identical(corrected,fixtures[[id]]))
        stop('Cannot resume changed scientific fixture: ',id)
    }
    manifest<-prior_manifest
    if(length(changed)) {
      manifest$call_corrections<-c(manifest$call_corrections,list(list(
        at=Sys.time(),cases=changed,original_manifest=snapshot,
        reason='Remove contradictory use_cutoff=FALSE from unfitted explicitly declared rich structures.')))
      manifest$fixtures<-fixtures
    }
    manifest$deadline<-Sys.time()+max_seconds
    manifest$status<-'running'
    manifest$resumed_at<-c(manifest$resumed_at,Sys.time())
    manifest$resume_snapshots<-c(manifest$resume_snapshots,snapshot)
  }
  save<-function()saveRDS(manifest,manifest_path)
  save()
  summaries<-diagnostics<-evidence<-policy<-list()
  keep<-function(expr)tryCatch(withCallingHandlers(expr,warning=function(w){
    manifest$warnings<<-c(manifest$warnings,conditionMessage(w));save();invokeRestart('muffleWarning')
  }),error=identity)
  scheduled<-names(fixtures)
  if(!isTRUE(include_ces_noiseless_stress))scheduled<-setdiff(scheduled,'ces_C_M4T4')
  manifest$scheduled_cases<-scheduled;save()
  for(id in scheduled)for(backend in c('cmdstanr','nutpie')) {
    if(Sys.time()>manifest$deadline){manifest$status<-'budget_exhausted';save();return(invisible(manifest))}
    key<-paste(id,backend,sep='_'); directory<-file.path(root,key)
    dir.create(directory,showWarnings=FALSE)
    args<-fixtures[[id]];args$backend<-backend;args$stan_options<-controls
    args$stan_options$seed<-controls$seed+match(id,names(fixtures))
    if(backend=='cmdstanr'){
      args$stan_options$output_dir<-directory
      args$stan_options$sig_figs<-18L
    }
    fit_path<-file.path(directory,'fit.rds')
    if(file.exists(fit_path)) {
      fit<-readRDS(fit_path)
      recorded<-manifest$attempts[[key]]
      if(!identical(digest::digest(file=fit_path,algo='sha256'),recorded$fit_sha256) ||
         !identical(nutpie_validation_identity(fit),recorded$identity) ||
         !identical(fit@metadata$stan_source_hash,
                    digest::digest(file=.cmdstan_kernel_file(fit@metadata$stan_kernel),algo='sha256')))
        stop('Saved fit identity differs; cannot resume: ',key)
      validate_fit_spec(fit)
      manifest$attempts[[key]]$fit_reused<-TRUE
      cat('REUSE verified fit',key,format(Sys.time()),'\n');flush.console()
    } else {
      manifest$attempts[[key]]<-list(started=Sys.time(),status='fitting');save()
      cat('FIT',key,format(Sys.time()),'\n');flush.console()
      start<-proc.time()[['elapsed']];fit<-keep(do.call(bayes_est,args))
      manifest$attempts[[key]]$seconds<-proc.time()[['elapsed']]-start
      if(inherits(fit,'error')){
        manifest$attempts[[key]]$status<-'failed'
        manifest$attempts[[key]]$reason<-conditionMessage(fit);save();next
      }
      saveRDS(fit,fit_path)
      manifest$attempts[[key]]$fit_sha256<-digest::digest(file=fit_path,algo='sha256')
      manifest$attempts[[key]]$identity<-nutpie_validation_identity(fit)
    }
    manifest$attempts[[key]]$status<-'fit_complete';save()
    parameter<-keep(nutpie_validation_parameters(fit))
    if(!inherits(parameter,'error'))summaries[[key]]<-cbind(case=id,backend,parameter$summary)
    diag<-keep(nutpie_validation_diagnostics(fit))
    if(!inherits(diag,'error'))diagnostics[[key]]<-cbind(case=id,backend,diag$summary)
    cat('BRIDGE',key,format(Sys.time()),'\n');flush.console()
    set.seed(452991L)
    evidence_path<-file.path(directory,'evidence.rds')
    comparison<-if(file.exists(evidence_path))readRDS(evidence_path)else
      keep(bayes_compare(candidate=fit,bridge=TRUE,bridge_repetitions=3L,bridge_method='normal'))
    if(!file.exists(evidence_path))saveRDS(comparison,evidence_path)
    if(!inherits(comparison,'error'))evidence[[key]]<-cbind(case=id,backend,comparison) else
      manifest$attempts[[key]]$evidence_error<-conditionMessage(comparison)
    # Full chains give chain-aware policy MCSE. Retain all attrition.
    if(run_policy && id %in% c('logit_C','logit_CF','logit_CS','logit_CF_k1','logit_CF_k3','ces_C','ces_CS')){
      cat('POLICY',key,format(Sys.time()),'\n');flush.console()
      policy_path<-file.path(directory,'policy.rds')
      pol<-if(file.exists(policy_path))readRDS(policy_path)else
        keep(nutpie_validation_policy(fit,seed=491871L+match(id,names(fixtures))+100L*(backend=="nutpie")))
      if(!file.exists(policy_path))saveRDS(pol,policy_path)
      if(!inherits(pol,'error'))policy[[key]]<-cbind(case=id,backend,pol$summary) else
        manifest$attempts[[key]]$policy_error<-conditionMessage(pol)
    }
    for(name in c('summaries','diagnostics','evidence','policy')){
      table<-get(name);if(length(table))write.csv(nutpie_validation_bind(table),file.path(root,paste0(name,'.csv')),row.names=FALSE)
    }
    manifest$attempts[[key]]$status<-'attempt_complete';save()
  }
  complete<-length(manifest$attempts)==2L*length(fixtures) &&
    all(vapply(manifest$attempts,function(x)identical(x$status,'attempt_complete'),logical(1)))
  manifest$status<-if(complete)'completed_attempts'else'incomplete_attempts'
  manifest$ended<-Sys.time();save()
  invisible(manifest)
}
