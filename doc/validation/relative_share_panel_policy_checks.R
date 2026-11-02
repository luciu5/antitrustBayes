# Bounded no-HMC structural checks using saved final M4T4 posterior fits.
# Sourcing defines functions only. Select compatible libraries and load current
# development source before calling relative_panel_policy_checks().
relative_panel_policy_checks <- function(root, output_dir, n_draws=4L,
    cases=c("logit_C_M4T4","ces_C_M4T4")) {
  rows<-list();warnings<-character();plans<-list();realization_failures<-list()
  add<-function(case,backend,unit,draw,test,error,tolerance=1e-8,reason='') {
    rows[[length(rows)+1L]]<<-data.frame(case,backend,unit,draw,test,error,tolerance,
      status=if(is.finite(error)&&error<=tolerance)'passed'else'failed',reason)
  }
  quiet<-function(expr)withCallingHandlers(expr,warning=function(w){
    warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  stopifnot(length(cases)>0L,all(cases %in% c('logit_C_M4T4','ces_C_M4T4')))
  for(case in cases)for(backend in c('cmdstanr','nutpie')) {
    path<-file.path(root,paste(case,backend,sep='_'),'fit.rds')
    if(!file.exists(path)){add(case,backend,'all',NA_integer_,'saved_fit_available',Inf,0,'Missing fit');next}
    fit<-readRDS(path);validate_fit_spec(fit)
    current<-digest::digest(file=.cmdstan_kernel_file(fit@metadata$stan_kernel),algo='sha256')
    if(!identical(current,fit@metadata$stan_source_hash))stop('Panel fit source mismatch.')
    plan<-posterior_policy_draws(fit,n_draws=n_draws,seed=453217L+(backend=='nutpie'))
    plans[[paste(case,backend)]]<-plan
    posterior<-.posterior_extract(fit,'policy_plan')
    v<-posterior$relative_share_predictor[plan$draw_indices,,drop=FALSE]
    if(fit@stan_data$quality_parameterization==2L)
      v<-v+plan$quality_draws[,as.integer(fit@stan_data$obs_to_quality),drop=FALSE]
    for(position in seq_along(plan$draw_indices))for(cell in seq_len(fit@stan_data$N_market_year)) {
      index<-which(as.integer(fit@stan_data$market_year_idx)==cell)
      utility<-v[position,index]+plan$derived_market_year_intercept_draws[position,cell]
      s0<-posterior$s0_prob[plan$draw_indices[position],cell]
      outside<-1/(1+sum(exp(utility)))
      add(case,backend,as.character(cell),plan$draw_indices[position],
        'same_quality_outside_calibration',abs(outside-s0),1e-12)
    }
    selected_years<-unique(as.character(fit@data$year))
    observed<-quiet(posterior_policy_models(fit,draw_plan=plan,year=selected_years,
      insideSize=1000,demand_anchoring='observed'))
    latent<-quiet(posterior_policy_models(fit,draw_plan=plan,year=selected_years,
      insideSize=1000,demand_anchoring='latent_quality'))
    realization_failures[[paste(case,backend)]]<-list(observed=observed$failures,latent=latent$failures)
    add(case,backend,'all',NA_integer_,'sixteen_market_period_units',abs(nrow(latent$units)-16L),0)
    for(unit in seq_len(nrow(latent$units)))for(position in seq_along(plan$draw_indices)) {
      uid<-as.character(latent$units$unit_id[unit]);draw<-plan$draw_indices[position]
      source<-latent$models[[uid]][[position]];anchor<-observed$models[[uid]][[position]]
      if(is.null(source)||is.null(anchor)) {
        failed<-latent$failures
        why<-if(nrow(failed))paste(unique(failed$message),collapse=';')else'Baseline realization unavailable'
        add(case,backend,uid,draw,'baseline_realization',Inf,0,why);next
      }
      model<-if(methods::is(source,'AntitrustFit'))source@model else source
      original<-if(methods::is(anchor,'AntitrustFit'))anchor@model else anchor
      add(case,backend,uid,draw,'physical_cost_preserved',max(abs(model@mcPre-original@mcPre)),1e-10)
      add(case,backend,uid,draw,'baseline_supply_foc',
        max(abs(as.numeric(quiet(antitrust::calcMC(model,preMerger=TRUE)))-model@mcPre)),1e-6)
      index<-which(as.character(fit@data$market)==as.character(latent$units$market[unit]) &
        as.character(fit@data$year)==as.character(latent$units$year[unit]))
      cell<-unique(as.integer(fit@stan_data$market_year_idx[index]))
      reference<-model;reference@pricePre<-as.numeric(fit@data$prices[index])
      target<-exp(v[position,index]-max(v[position,index]));target<-target/sum(target)
      actual<-as.numeric(quiet(antitrust::calcShares(reference,preMerger=TRUE)))
      actual<-actual/sum(actual)
      add(case,backend,uid,draw,'reference_relative_shares',max(abs(actual-target)),1e-8)
      if(uid %in% c('m1_t1','m4_t4')) {
        owner<-as.character(fit@data$group[index]);firms<-unique(owner)
        for(regime in c('merger','tariff25')) {
          check<-tryCatch({
            promoted<-trade::as_trade_fit(source,policy='tariff',tariffPre=0,
              cost_basis='effective',margin_basis='net_revenue')
            owner_post<-owner;if(regime=='merger')owner_post[owner_post==firms[2L]]<-firms[1L]
            solved<-quiet(trade::simulate(promoted,ownerPost=owner_post,
              tariffPost=if(regime=='merger')0 else ifelse(owner %in% firms[1:2],0,.25/1.25)))
            pre<-trade::tariff_accounts(solved,preMerger=TRUE)
            post<-trade::tariff_accounts(solved,preMerger=FALSE)
            if(any(!is.finite(post$price)))stop('Nonfinite counterfactual price')
            max(abs(pre$physical_cost-post$physical_cost))
          },error=identity)
          add(case,backend,uid,draw,paste0(regime,'_physical_cost'),
            if(inherits(check,'error'))Inf else check,1e-10,
            if(inherits(check,'error'))conditionMessage(check)else'')
        }
      }
    }
  }
  table<-do.call(rbind,rows)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(table,file.path(output_dir,'relative_share_panel_policy_checks.csv'),row.names=FALSE)
  result<-list(checks=table,warnings=unique(warnings),plans=plans,
    realization_failures=realization_failures,
    created=Sys.time(),source_hashes=vapply(list.files('inst/stan','[.]stan$',full.names=TRUE),
      digest::digest,character(1),file=TRUE,algo='sha256'))
  saveRDS(result,file.path(output_dir,'panel-policy-checks.rds'))
  invisible(result)
}
