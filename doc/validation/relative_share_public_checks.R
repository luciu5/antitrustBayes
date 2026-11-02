# Saved-fit public comparison and cache-restoration checks; no HMC.
relative_public_mixed_comparison <- function(root, output_dir) {
  models<-c('MC','B','C','BF','CF','BS','CS')
  backends<-rep(c('nutpie','cmdstanr'),length.out=length(models))
  fits<-stats::setNames(lapply(seq_along(models),function(i)
    readRDS(file.path(root,paste('logit',models[i],backends[i],sep='_'),'fit.rds'))),models)
  hypotheses<-list(MC='MC',B='B',C='C',BF=list(fit='BF',k='interior'),
    CF=list(fit='CF',k='interior'),BS=list(fit='BS',k='interior'),CS=list(fit='CS',k='interior'))
  universe<-bayes_model_universe(hypotheses,prior_weights=rep(1/7,7))
  set.seed(453823L)
  comparison<-bayes_compare(fits,universe=universe,bridge=TRUE,bridge_repetitions=1L)
  stopifnot(all(is.finite(comparison$logml)),all(is.finite(comparison$post_prob)),
    abs(sum(comparison$post_prob)-1)<1e-12)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(comparison,file.path(output_dir,'public-mixed-comparison.rds'))
  write.csv(comparison,file.path(output_dir,'public_mixed_comparison.csv'),row.names=FALSE)
  invisible(comparison)
}

relative_saved_nutpie_checks <- function(root, output_dir) {
  path<-file.path(root,'logit_CF_nutpie','fit.rds')
  fit<-readRDS(path);validate_fit_spec(fit)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  serialized<-file.path(output_dir,'serialized-fit.rds')
  saveRDS(fit,serialized);reloaded<-readRDS(serialized)
  stopifnot(identical(.bayes_backend_draws(fit@stanfit),.bayes_backend_draws(reloaded@stanfit)))
  before<-posterior_policy_draws(fit,n_draws=4L,seed=453917L)
  after<-posterior_policy_draws(reloaded,n_draws=4L,seed=453917L)
  stopifnot(identical(before,after))
  original_fit_hash<-digest::digest(file=path,algo='sha256')
  initial<-.bayes_bridge_restore_context(fit)
  metadata<-fit@stanfit$metadata
  library<-metadata$lib_path
  expected<-metadata$lib_hash
  stopifnot(file.exists(library),identical(digest::digest(file=library,algo='sha256'),expected))
  backup<-paste0(library,'.relative-validation-backup')
  if(file.exists(backup))stop('An earlier cache backup exists; preserve it.')
  stopifnot(file.rename(library,backup))
  on.exit({
    if(file.exists(library))unlink(library)
    if(!file.rename(backup,library))stop('Could not restore original validation library.')
  },add=TRUE)
  restored<-.bayes_bridge_restore_context(fit)
  stopifnot(isTRUE(restored$migration),
    identical(initial$constrained_parameter_names,restored$constrained_parameter_names),
    identical(initial$samples,restored$samples))
  rows<-unique(round(seq(1,nrow(initial$samples),length.out=8L)))
  original<-vapply(rows,function(i)initial$log_posterior(initial$samples[i,],NULL),numeric(1))
  reconstructed<-vapply(rows,function(i)restored$log_posterior(restored$samples[i,],NULL),numeric(1))
  error<-max(abs(original-reconstructed))
  stopifnot(is.finite(error),error<1e-8,
    identical(digest::digest(file=path,algo='sha256'),original_fit_hash))
  result<-list(serialized_posterior_identical=TRUE,serialized_plan_identical=TRUE,
    maximum_reconstructed_density_error=error,original_context=initial$restoration,
    reconstructed_context=restored$restoration,original_fit_sha256=original_fit_hash,
    created=Sys.time())
  saveRDS(result,file.path(output_dir,'saved-fit-cache-checks.rds'))
  invisible(result)
}
