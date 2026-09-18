#' Realize reusable structural models from a posterior policy draw plan
#'
#' Constructs one baseline structural model per selected draw and market-period
#' unit. The returned models can be passed to a deterministic policy package
#' without fitting again or simulating an unchanged counterfactual.
#'
#' @param fit A \code{BayesFit} for a supported output-market model.
#' @param draw_plan A plan returned by \code{posterior_policy_draws()}. If
#'   omitted, a new plan is constructed once.
#' @param n_draws Number of draws when constructing a plan; must agree with
#'   the supplied plan when both are specified.
#' @param seed Optional seed for a new plan. Cannot accompany a supplied plan.
#' @param market,year Simulation-unit selectors. Panel fits require an explicit
#'   year selection, as in \code{simulate()}.
#' @param demand_anchoring \code{"observed"} or \code{"latent_quality"}.
#'   Latent anchoring uses the plan's conditional quality realization and
#'   preserves the original structural cost primitive. Solving the resulting
#'   latent baseline may be necessary.
#' @param insideSize,priceOutside,labels Optional simulation-unit overrides,
#'   following the conventions of \code{simulate()}.
#' @return A \code{bayes_policy_models} list with \code{models} nested by
#'   unit and draw position, \code{draw_indices}, \code{units}, the reusable
#'   \code{draw_plan}, and a \code{failures} data frame. Failed models occupy
#'   their original list positions as \code{NULL}.
#' @details Initial support includes simultaneous Bertrand, Cournot and
#'   monopolistic competition, and native exact fringe/Stackelberg models under
#'   output Logit/CES demand. Ordinary models are \code{AntitrustFit} objects;
#'   exact role models are native \pkg{coordination} objects. No policy wedge
#'   or tariff equation is evaluated here. Attributes retain the draw index,
#'   local quality units, realized k and plan fingerprint. Quality and roles
#'   are never resampled during realization.
#'   Second-score output Logit auctions use the cost-based value index of
#'   \pkg{antitrust}. Their winning probabilities respond to effective costs
#'   rather than bid prices; realized quality is held fixed across policies.
#' @seealso \code{posterior_policy_draws()}, \code{posterior_quality()}
#' @export
posterior_policy_models <- function(fit, draw_plan = NULL, n_draws = NULL,
                                    seed = NULL, market = NULL, year = NULL,
                                    demand_anchoring = c("observed", "latent_quality"),
                                    insideSize = NULL, priceOutside = NULL,
                                    labels = NULL) {
  if (!methods::is(fit, "BayesFit"))
    stop("'fit' must be a BayesFit object.", call. = FALSE)
  demand_anchoring <- match.arg(demand_anchoring)
  specified <- methods::is(fit, "BayesSpecifiedFit")
  auction_logit <- identical(as.character(fit@demand_model), "logit") &&
    identical(as.integer(fit@supply_model), 4L)
  if (!identical(fit@market_side, "output") ||
      !fit@demand_model %in% c("logit", "ces") ||
      (!as.integer(fit@supply_model) %in% c(1L, 2L, 3L) && !auction_logit))
    stop("posterior_policy_models() supports output Logit/CES Bertrand, Cournot and MonCom games, and second-score output Logit auctions.",
         call. = FALSE)
  native_structure <- .native_ordered_structure(fit)
  structure_name <- fit@metadata[["structure", exact = TRUE]] %||% if (specified)
    fit@params[["structure", exact = TRUE]] %||% "simultaneous" else
    attr(fit@stan_data, "structure", exact = TRUE) %||%
      if (isTRUE(fit@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous"
  if (is.null(native_structure) && !identical(structure_name, "simultaneous"))
    stop("Smooth fringe is not supported by posterior_policy_models().", call. = FALSE)
  if (!is.null(draw_plan) && !is.null(seed))
    stop("'seed' applies only when creating a new draw_plan.", call. = FALSE)
  legacy_plan <- is.list(draw_plan) && is.list(draw_plan$fit_signature) &&
    (is.null(draw_plan$fit_signature$fit_signature_version) ||
     identical(as.integer(draw_plan$fit_signature$fit_signature_version), 1L))
  draws <- if (specified) NULL else if (legacy_plan) rstan::extract(fit@stanfit)
    else .posterior_extract(fit, "policy_plan")
  plan <- if (is.null(draw_plan))
    .lq_make_policy_plan(fit, draws, n_draws = n_draws,
                         seed = .lq_validate_seed(seed)) else
    .lq_validate_plan(fit, draw_plan, n_draws = n_draws, draws = draws)
  setup <- if (specified) .specified_simulation_units(fit, market, year) else
    .estimated_simulation_units(fit, market, year)
  indices <- as.integer(plan$draw_indices)
  models <- stats::setNames(vector("list", length(setup$info)), setup$units$unit_id)
  failures <- list()
  incidence <- .lq_quality_incidence(fit)
  for (u in seq_along(setup$info)) {
    unit <- setup$info[[u]]
    idx <- unit$idx
    df <- setup$data[idx, , drop = FALSE]
    uid <- unit$unit_id
    choose <- function(value, default, scalar = FALSE, name) {
      if (is.null(value)) default else .select_unit_arg(
        value, idx, uid, unit$market_id, nrow(setup$data),
        scalar = scalar, arg = name)
    }
    default_labels <- if (specified) fit@market$labels else {
      columns <- intersect(c("labels", "product", "product_id"), names(df))
      if (length(columns)) as.character(df[[columns[[1L]]]]) else
        paste0("product_", seq_len(nrow(df)))
    }
    default_labels <- make.unique(as.character(default_labels))
    market_object <- bayes_market(
      prices = as.numeric(df$prices), quantities = as.numeric(df$quantities),
      margins = as.numeric(df$margins), group = as.character(df$group),
      outside_price = choose(priceOutside,
        if (specified) fit@market$outside_price else unique(df$outside_price),
        scalar = TRUE, name = "priceOutside"),
      market_side = fit@market_side,
      insideSize = choose(insideSize, if (specified) fit@market$insideSize else 1,
                           scalar = TRUE, name = "insideSize"),
      labels = choose(labels, default_labels, name = "labels"))
    bridge_unit <- c(list(data = df, market = unit$market_id,
                         year = unit$year_id), unit)
    q_index <- if (is.null(incidence)) match(as.character(df$group), plan$group_labels)
      else incidence$obs_to_quality[idx]
    local_q <- unique(q_index)
    if (anyNA(local_q)) stop("Unit quality map is not in the draw plan.", call. = FALSE)
    values <- stats::setNames(vector("list", length(indices)), as.character(indices))
    for (pos in seq_along(indices)) {
      draw <- indices[[pos]]
      stage <- "baseline_construction"
      result <- tryCatch({
        params <- if (specified) fit@params else
          .estimated_draw_params(fit, draws, draw, unit)
        k <- if (is.null(native_structure)) NULL else plan$k_cache[[pos]][[uid]]
        if (!is.null(native_structure)) {
          if (length(k) != 1L || !is.finite(k))
            stop("The draw plan has no role realization for this unit.")
          constructor <- if (native_structure == "stackelberg")
            .stack_bridge_model_draw else .core_fringe_bridge_model_draw
          model <- constructor(
            fit, unit = bridge_unit, draw = draw, role_mapping = "map",
            fixed_k = k, role_offset = 0L, insideSize = market_object$insideSize,
            priceOutside = market_object$outside_price, labels = market_object$labels,
            draws = draws, dots = list())
        } else {
          model <- as_antitrust_fit(build_conduct_object(
            params, if (specified) draw else 1L, market_object,
            ownerPre = market_object$group, ownerPost = market_object$group,
            mcDelta = rep(0, market_object$n_prod),
            subset = rep(TRUE, market_object$n_prod)))
        }
        latent_metadata <- NULL
        if (demand_anchoring == "latent_quality") {
          stage <- "latent_pre_equilibrium"
          score <- .lq_reference_score(fit, draws, draw, bridge_unit, plan,
                                       market_object, params)
          latent <- if (is.null(native_structure)) .lq_apply_to_fit(model, score, plan)
            else .lq_apply_to_native(model, score, plan)
          model <- if (is.null(native_structure)) latent$fit else latent$native
          latent_metadata <- latent$metadata
          latent_metadata$quality <- latent_metadata$quality[local_q]
          attr(model, "bayes_latent_quality") <- latent_metadata
        }
        attr(model, "bayes_policy_draw") <- list(
          draw_index = draw, unit_id = uid, market = unit$market_id,
          year = unit$year_id, k = k, realized_k = k,
          demand_anchoring = demand_anchoring,
          quality_scope = plan$quality_scope,
          quality_units = plan$quality_units[local_q],
          quality = as.numeric(plan$quality_draws[pos, local_q]),
          plan_fingerprint = plan$fit_signature$latent_state_fingerprint)
        list(model = model, error = NULL)
      }, error = function(e) list(model = NULL, error = e))
      values[pos] <- list(result$model)
      if (!is.null(result$error)) failures[[length(failures) + 1L]] <- data.frame(
        unit_id = uid, market = unit$market_id, year = unit$year_id,
        draw_index = draw, stage = stage,
        category = class(result$error)[[1L]], message = conditionMessage(result$error),
        stringsAsFactors = FALSE)
    }
    models[u] <- list(values)
  }
  failure_df <- if (length(failures)) do.call(rbind, failures) else data.frame(
    unit_id = character(), market = character(), year = character(),
    draw_index = integer(), stage = character(), category = character(),
    message = character(), stringsAsFactors = FALSE)
  if (nrow(failure_df)) warning(nrow(failure_df), " structural draw realizations failed; inspect $failures.",
                               call. = FALSE)
  structure(list(models = models, draw_indices = indices, units = setup$units,
                 draw_plan = plan, demand_anchoring = demand_anchoring,
                 failures = failure_df), class = c("bayes_policy_models", "list"))
}
