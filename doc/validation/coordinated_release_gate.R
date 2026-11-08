#!/usr/bin/env Rscript
# Extended tier: exactly two fresh spec-14 fits, no bridge sampling or retuning.
# Source this file to inspect helpers without sampling. From repository root:
# R_LIBS=/compatible/library Rscript --vanilla doc/validation/coordinated_release_gate.R --sample
# After sampling, select the final core library and run --analyze --policy.
# Durable output is intentionally outside this staged source tree.
release_gate_hash <- function(path) digest::digest(file = path, algo = 'sha256')
release_gate_sources <- function() {
  paths <- c('DESCRIPTION', sort(list.files('R', '[.]R$', full.names = TRUE)),
    sort(list.files('inst/stan', '[.]stan$', full.names = TRUE)),
    paste0('doc/validation/', c('coordinated_release_gate.R',
      'relative_share_backend_validation.R', 'relative_share_noisy_ces_panel.R',
      'single_market_intercept_validation.R', 'nutpie_backend_validation.R')))
  stats::setNames(vapply(paths, release_gate_hash, character(1)), paths)
}
release_gate_packages <- function() {
  lapply(stats::setNames(c('antitrust','coordination','cmdstanr','posterior','pkgload'),
    c('antitrust','coordination','cmdstanr','posterior','pkgload')), function(p) {
    path <- find.package(p)
    files <- c(file.path(path,'DESCRIPTION'), list.files(file.path(path,'R'),
      '[.](rdb|rdx)$',full.names=TRUE))
    list(version=as.character(utils::packageVersion(p)), path=path,
      files=stats::setNames(vapply(files,release_gate_hash,character(1)),files))
  })
}
release_gate_prepare <- function(root) {
  if(file.exists(file.path(root,'manifest.rds'))) stop('Existing gate output is preserved; use --analyze, never refit it.')
  source('doc/validation/relative_share_backend_validation.R', local=TRUE)
  source('doc/validation/relative_share_noisy_ces_panel.R', local=TRUE)
  warnings <- character()
  fixtures <- withCallingHandlers(relative_validation_fixtures(), warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')
  })
  noisy <- relative_noisy_ces_panel_arguments(fixtures$ces_C_M4T4)
  calls <- list(logit_C_M4T4=fixtures$logit_C_M4T4,
    ces_C_M4T4_finite_noise=noisy$arguments)
  # Recreate the current Logit truth before its prescribed observation errors.
  source('doc/validation/single_market_intercept_validation.R', local=TRUE)
  truth_logit <- suppressWarnings(smi_fixture(data.frame(id='C_M4T4',conduct='C',
    markets=4L,periods=4L),smi_design()))$truth
  truth <- list(logit_C_M4T4=truth_logit,
    ces_C_M4T4_finite_noise=as.data.frame(fixtures$ces_C_M4T4[c(
      'market','year','group','prices','quantities','margins')]))
  seeds <- c(logit_C_M4T4=726041L,ces_C_M4T4_finite_noise=726042L)
  for(id in names(calls)) {
    calls[[id]]$share_basis <- 'relative'; calls[[id]]$backend <- 'cmdstanr'
    calls[[id]]$stan_options <- list(chains=4L,iter=1000L,warmup=500L,cores=4L,
      refresh=100L,seed=unname(seeds[id]),sig_figs=18L,
      output_dir=file.path(root,id),control=list(adapt_delta=.99,max_treedepth=13L))
  }
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  manifest <- list(created=Sys.time(),repository=normalizePath('.'),
    package_version='0.2.0.9009',model_specification=14L,
    baseline_revision='fb1f026',source_hashes=release_gate_sources(),
    packages=release_gate_packages(),libraries=.libPaths(),R=R.version.string,
    cmdstan=list(path=cmdstanr::cmdstan_path(),version=as.character(cmdstanr::cmdstan_version())),
    calls=calls,calls_hash=digest::digest(calls,algo='sha256'),truth=truth,
    ces_noise=noisy,preparation_warnings=warnings,attempts=list(),status='prepared',
    thresholds=list(structural_coverage=.95,maximum_relative_median_error=.08,
      s0_minimum_coverage=3L,s0_markets=4L,s0_reference_year='t1',maximum_s0_mae=.06,
      maximum_rhat=1.05,minimum_bulk_ess=100,divergences=0L,
      noise_rms_ratio=c(.5,2),maximum_foc_residual=1e-6,
      fixed_policy_draw=local({set.seed(726051L);sample.int(2000L,1L)}),policy_plan_seed=726051L),
    scope=paste('Exactly two current relative-share panel fits. Structural recovery over all four markets;',
      'outside-share gate at the prespecified t1 reference in all four markets (all 16 cells reported).',
      'Noise-scale medians must lie between 0.5 and 2 times actual prescribed log-error RMS.',
      'All active sampled parameters must pass diagnostics; no derived constant quantities enter extrema.',
      'No conduct selection, bridge sampling, repeated-sample coverage or 124-task benchmark.'))
  saveRDS(manifest,file.path(root,'manifest.rds'))
  write.csv(noisy$observation_table,file.path(root,'ces-observation-errors.csv'),row.names=FALSE)
  manifest
}
release_gate_assert <- function(manifest) {
  stopifnot(identical(manifest$source_hashes,release_gate_sources()),
    identical(manifest$calls_hash,digest::digest(manifest$calls,algo='sha256')))
}
release_gate_sample <- function(root) {
  m <- release_gate_prepare(root)
  save <- function()saveRDS(m,file.path(root,'manifest.rds'))
  on.exit({if(!identical(m$status,'sampling_complete')) {m$status <- 'interrupted';save()}},add=TRUE)
  for(id in names(m$calls)) {
    release_gate_assert(m);dir.create(file.path(root,id),showWarnings=FALSE)
    m$attempts[[id]] <- list(status='fitting',started=Sys.time(),warnings=character());save()
    cat('FIT',id,format(Sys.time()),'\n');flush.console();started <- proc.time()[['elapsed']]
    fit <- tryCatch(withCallingHandlers(do.call(bayes_est,m$calls[[id]]),warning=function(w) {
      m$attempts[[id]]$warnings <<- c(m$attempts[[id]]$warnings,conditionMessage(w));save()
      invokeRestart('muffleWarning')
    }),error=identity)
    m$attempts[[id]]$seconds <- proc.time()[['elapsed']]-started
    path <- file.path(root,id,if(inherits(fit,'error'))'fit-error.rds' else 'fit.rds')
    saveRDS(fit,path);m$attempts[[id]]$artifact_sha256 <- release_gate_hash(path)
    m$attempts[[id]]$status <- if(inherits(fit,'error'))'failed' else 'fit_complete'
    if(inherits(fit,'error'))m$attempts[[id]]$reason <- conditionMessage(fit)
    save();cat(m$attempts[[id]]$status,id,format(Sys.time()),'\n');flush.console()
  }
  m$status <- 'sampling_complete';m$ended <- Sys.time();save();invisible(m)
}
release_gate_fit <- function(root,m,id) {
  release_gate_assert(m);path <- file.path(root,id,'fit.rds')
  stopifnot(file.exists(path),identical(release_gate_hash(path),m$attempts[[id]]$artifact_sha256))
  fit <- readRDS(path);getFromNamespace('validate_fit_spec','antitrustBayes')(fit)
  md <- fit@metadata;kernel <- paste0(m$calls[[id]]$demand_model,'_supply')
  stopifnot(identical(as.integer(md$model_spec_version),14L),identical(md$share_basis,'relative'),
    identical(md$backend,'cmdstanr'),identical(md$stan_kernel,kernel),
    identical(md$stan_source_hash,unname(m$source_hashes[paste0('inst/stan/',kernel,'.stan')])),
    identical(md$stan_data_hash,digest::digest(fit@stan_data[md$stan_data_fields],algo='sha256')),
    identical(dim(fit@stanfit$draws())[1:2],c(500L,4L)))
  for(field in c('prices','quantities','margins','market','year','group'))
    stopifnot(identical(as.character(fit@data[[field]]),as.character(m$calls[[id]][[field]])))
  fit
}
release_gate_active <- function(fit) {
  lines <- readLines(file.path('inst/stan',paste0(fit@metadata$stan_kernel,'.stan')))
  begin <- which(lines=='parameters {');end <- which(seq_along(lines)>begin & lines=='}')[1L]
  declarations <- strsplit(paste(sub('//.*$','',lines[(begin+1L):(end-1L)]),collapse=' '),';',fixed=TRUE)[[1L]]
  bases <- trimws(sub('.*[[:space:]]([A-Za-z][A-Za-z0-9_]*)[[:space:]]*$','\\1',trimws(declarations)))
  cols <- posterior::variables(fit@stanfit$draws())
  cols[sub('\\[.*$','',cols)%in%bases]
}
release_gate_analyze <- function(root) {
  m <- readRDS(file.path(root,'manifest.rds'));rows <- diagnostics <- gates <- noise <- list()
  add <- function(id,test,value,threshold,passed) gates[[length(gates)+1L]] <<-
    data.frame(case=id,test,value,threshold,passed)
  for(id in names(m$calls)) {
    fit <- tryCatch(release_gate_fit(root,m,id),error=identity)
    if(inherits(fit,'error')) {add(id,'fit_identity',conditionMessage(fit),'verified current fit',FALSE);next}
    add(id,'fit_identity','verified','verified current fit',TRUE)
    pars <- if(fit@demand_model=='ces')'gamma_structural_market' else 'alpha_structural_market'
    draws <- as.matrix(fit@stanfit$draws(c(pars,'s0_prob','sigma_share','sigma_margin'),format='draws_matrix'))
    summary <- lapply(seq_len(ncol(draws)),function(j) {
      q <- quantile(draws[,j],c(.025,.5,.975));data.frame(case=id,variable=colnames(draws)[j],
        q025=q[1L],median=q[2L],q975=q[3L])
    });tab <- do.call(rbind,summary);tab$truth <- NA_real_
    structural <- startsWith(tab$variable,paste0(pars,'['))
    target <- if(fit@demand_model=='ces')rep(6,4L) else 2.4*exp(.06*(0:3))
    tab$truth[structural] <- target
    cells <- unique(data.frame(cell=fit@stan_data$market_year_idx,market=fit@data$market,year=fit@data$year))
    cells <- cells[order(cells$cell),];true <- m$truth[[id]]
    s0_truth <- if(fit@demand_model=='ces')rep(.2,nrow(cells)) else vapply(seq_len(nrow(cells)),function(j)
      unique(true$s0[true$market==cells$market[j]&true$year==cells$year[j]]),numeric(1))
    s0rows <- match(paste0('s0_prob[',cells$cell,']'),tab$variable);tab$truth[s0rows] <- s0_truth
    tab$covered <- with(tab,truth>=q025&truth<=q975);tab$relative_error <- abs(tab$median/tab$truth-1)
    rows[[id]] <- tab
    add(id,'all_structural_truth_in_95pct',sum(tab$covered[structural]),'4/4',all(tab$covered[structural]))
    add(id,'max_structural_relative_median_error',max(tab$relative_error[structural]),'.08',max(tab$relative_error[structural])<=.08)
    ref <- s0rows[as.character(cells$year)=='t1'];stopifnot(length(ref)==4L)
    coverage <- sum(tab$covered[ref]);mae <- mean(abs(tab$median[ref]-tab$truth[ref]))
    add(id,'s0_t1_truth_in_95pct',coverage,'>=3/4',coverage>=3L)
    add(id,'s0_t1_median_mae',mae,'<=.06',mae<=.06)
    active <- release_gate_active(fit);stopifnot(length(active)>0L)
    diag <- as.data.frame(posterior::summarise_draws(fit@stanfit$draws(active)))
    diagnostics[[id]] <- cbind(case=id,diag)
    valid <- all(is.finite(diag$rhat)&is.finite(diag$ess_bulk))
    rh <- if(valid)max(diag$rhat)else Inf;ess <- if(valid)min(diag$ess_bulk)else 0
    div <- sum(fit@stanfit$sampler_diagnostics()[,,'divergent__'])
    add(id,'max_active_parameter_rhat',rh,'<=1.05',rh<=1.05)
    add(id,'min_active_parameter_bulk_ess',ess,'>=100',ess>=100)
    add(id,'divergences',div,'0',div==0)
    args <- m$calls[[id]]
    if(fit@demand_model=='ces') {
      se <- m$ces_noise$observation_table$prescribed_log_quantity_error
      me <- m$ces_noise$observation_table$prescribed_log_margin_error
    } else {
      j <- rep(1:7,16);mi <- rep(rep(1:4,each=4),each=7);tt <- rep(rep(1:4,4),each=7)
      se <- .01*sin(j+2*mi+tt);me <- .01*cos(j+mi+2*tt)
    }
    # The relative likelihood removes one common log-share direction per cell.
    cell <- interaction(args$market,args$year,drop=TRUE)
    centered <- se-ave(se,cell,FUN=mean)
    rms <- c(sigma_share=sqrt(sum(centered^2)/(length(se)-nlevels(cell))),sigma_margin=sqrt(mean(me^2)))
    for(v in names(rms)) {
      med <- tab$median[tab$variable==v];ratio <- med/rms[v]
      noise[[paste(id,v)]] <- data.frame(case=id,variable=v,actual_error_rms=unname(rms[v]),posterior_median=med,ratio=unname(ratio))
      add(id,paste0(v,'_median_over_actual_rms'),unname(ratio),'[.5,2]',is.finite(ratio)&&ratio>=.5&&ratio<=2)
    }
  }
  tables <- list(recovery=rows,active_parameter_diagnostics=diagnostics,noise=noise,gates=gates)
  tables <- lapply(tables,function(x)if(length(x))do.call(rbind,x)else data.frame())
  for(n in names(tables))write.csv(tables[[n]],file.path(root,paste0(n,'.csv')),row.names=FALSE)
  saveRDS(tables,file.path(root,'analysis.rds'));tables
}
# Independent FOCs from primitive demand, costs and ownership: no calcMC,
# calcMargins, calcShares or package elasticity helper is used as an oracle.
# CES follows the legacy Cournot equation: calibrated outside-revenue alpha
# remains fixed over a policy, as required by the S4 behavioral contract.
release_gate_raw_foc <- function(object,pre=TRUE) {
  p <- if(pre)object@pricePre else object@pricePost
  cost <- if(pre)object@mcPre else object@mcPost
  own <- if(pre)object@ownerPre else object@ownerPost
  stopifnot(is.matrix(own),all(is.finite(p)),all(p>0),all(is.finite(cost)))
  if(methods::is(object,'CES')) {
    g <- object@slopes$gamma;w <- object@slopes$meanval
    weights <- w*p^(1-g);revenue <- weights/(object@priceOutside^(1-g)+sum(weights))
    a <- object@slopes$alpha
    derivative_loss <- p/g+(g-1)*(1+a)/(g*(1+g*a))*p*as.vector(own%*%revenue)
  } else {
    a <- object@slopes$alpha;utility <- object@slopes$meanval+a*(p-object@priceOutside)
    weights <- exp(utility);s0 <- 1/(1+sum(weights));q <- weights*s0
    # Inverse-demand derivative D_ji=(delta_ji/q_j + 1/s0)/a.
    inverse <- (diag(1/q)+matrix(1/s0,length(p),length(p)))/a
    derivative_loss <- -rowSums(own*t(inverse)*matrix(q,length(p),length(p),byrow=TRUE))
  }
  p-cost-derivative_loss
}
release_gate_policy <- function(root) {
  m <- readRDS(file.path(root,'manifest.rds'));release_gate_assert(m)
  rows <- list();artifacts <- list();warnings <- character()
  quiet <- function(expr)withCallingHandlers(expr,warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  })
  for(id in names(m$calls)) {
    fit <- release_gate_fit(root,m,id);truth <- m$truth[[id]];args <- m$calls[[id]]
    plan <- posterior_policy_draws(fit,n_draws=1L,seed=m$thresholds$policy_plan_seed)
    stopifnot(identical(plan$draw_indices,m$thresholds$fixed_policy_draw))
    prepared <- quiet(posterior_policy_models(fit,draw_plan=plan,year='t1',insideSize=1000,
      demand_anchoring='latent_quality'))
    artifacts[[id]] <- list(plan=plan,posterior=prepared)
    for(market in paste0('m',1:4))for(state in c('truth','posterior_fixed_draw')) {
      result <- tryCatch({
        if(state=='truth') {
          ix <- which(as.character(truth$market)==market&as.character(truth$year)=='t1')
          if(args$demand_model=='logit') {
            p <- truth$price[ix];q <- truth$share[ix];margin <- truth$margin[ix]
            owner <- as.character(truth$product[ix]);a <- unique(truth$alpha[ix]);s0 <- unique(truth$s0[ix])
            params <- bayes_params(alpha=a,s0=s0,supply_model=3L)
          } else {
            p <- truth$prices[ix];q <- truth$quantities[ix];margin <- truth$margins[ix]
            owner <- as.character(truth$group[ix]);params <- bayes_params(gamma=6,s0=.2,supply_model=3L)
          }
          market_data <- bayes_market(p,q,margin,group=owner,insideSize=1000,
            outside_price=if(is.null(args$outside_price))1 else args$outside_price,market_side='output')
          object <- quiet(posterior_policy_models(bayes_specify(params,market_data),n_draws=1L))$models[[1L]][[1L]]
        } else {
          uid <- as.character(prepared$units$unit_id[prepared$units$market==market])
          stopifnot(length(uid)==1L);object <- prepared$models[[uid]][[1L]]
          owner <- as.character(args$group[as.character(args$market)==market&as.character(args$year)=='t1'])
        }
        source_object <- object
        if(methods::is(object,'AntitrustFit'))object <- object@model
        stopifnot(!is.null(object))
        baseline <- max(abs(release_gate_raw_foc(object,TRUE)))
        owner_post <- owner;firms <- unique(owner);owner_post[owner_post==firms[2L]] <- firms[1L]
        solved <- quiet(antitrust::simulate(source_object,antitrust::counterfactual(ownership=owner_post)))
        if(methods::is(solved,'AntitrustFit'))solved <- solved@model
        stopifnot(any(rowSums(solved@ownerPost)>1L))
        merger <- max(abs(release_gate_raw_foc(solved,FALSE)))
        artifacts[[id]][[paste(market,state,sep='_')]] <- list(baseline=object,merger=solved)
        data.frame(case=id,market,state,regime=c('baseline','merger'),residual=c(baseline,merger),reason='')
      },error=function(e)data.frame(case=id,market,state,regime=c('baseline','merger'),residual=Inf,reason=conditionMessage(e)))
      rows[[length(rows)+1L]] <- result
    }
  }
  table <- do.call(rbind,rows);table$threshold <- 1e-6;table$passed <- is.finite(table$residual)&table$residual<=table$threshold
  write.csv(table,file.path(root,'policy-focs.csv'),row.names=FALSE)
  saveRDS(list(created=Sys.time(),rows=table,artifacts=artifacts,warnings=warnings,
    packages=release_gate_packages(),sources=release_gate_sources()),file.path(root,'policy.rds'))
  table
}
if(sys.nframe()==0L) {
  pkgload::load_all('.',quiet=TRUE)
  root <- Sys.getenv('AB_RELEASE_GATE_OUTPUT','/mnt/d/Projects/antitrustBayes/.scratch/coordinated-release-gate')
  options(antitrustBayes.cmdstan_cache=file.path(root,'stan-cache'))
  flags <- commandArgs(TRUE)
  stopifnot(length(flags)>0L,all(flags%in%c('--sample','--analyze','--policy')))
  if('--sample'%in%flags)release_gate_sample(root)
  if('--analyze'%in%flags) {
    tables <- release_gate_analyze(root)
    print(tables$gates,row.names=FALSE)
  }
  if('--policy'%in%flags) {policy <- release_gate_policy(root);print(policy,row.names=FALSE)}
  if('--analyze'%in%flags && (!all(tables$gates$passed)||('--policy'%in%flags&&!all(policy$passed))))quit(status=1L)
}
