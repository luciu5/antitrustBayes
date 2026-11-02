#' Carry rich-game firm roles through an ownership change
#'
#' Map pre-policy leaders or core firms to their post-policy owners without
#' reranking products. A merged firm retains a role when any constituent
#' pre-policy firm held it. The result can be passed directly to
#' `trade::simulate()` together with `ownerPost`.
#'
#' @param model A native `coordination` Stackelberg or exact core-fringe model,
#'   or a result from [bayes_calibrated_market()].
#' @param ownerPost Post-policy firm owner for each product in model order.
#' @return A one-element list named `leadersPost` or `corePost`.
#' @export
bayes_policy_roles <- function(model, ownerPost) {
  if (inherits(model, "bayes_calibrated_market")) model <- model$truth$model
  stack <- methods::is(model, "StackelbergLogit") ||
    methods::is(model, "StackelbergCES")
  fringe <- methods::is(model, "CoreFringeLogit") ||
    methods::is(model, "CoreFringeCES")
  if (stack || fringe) {
    pre_roles <- if (stack) model@leadersPre else model@corePre
    if (any(!as.character(pre_roles) %in% as.character(model@firmOwnerPre)))
      stop("Pre-policy firm roles do not match the model's ownership map.",
           call. = FALSE)
  }
  if (stack) {
    return(list(leadersPost = .stack_native_scenario(
      model, scenario = ownerPost)$leadersPost))
  }
  if (fringe) {
    return(list(corePost = .core_native_scenario(
      model, scenario = ownerPost)$corePost))
  }
  stop("'model' must be a native Stackelberg or exact core-fringe model.",
       call. = FALSE)
}
