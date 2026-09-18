## Deterministic simulation domain; this never modifies a fitted posterior.
.validate_simulation_price_domain <- function(object, price_domain) {
  price_domain <- match.arg(price_domain, c("positive", "real"))
  if (identical(price_domain, "real") &&
      (!identical(as.character(object@demand_model), "logit") ||
       !identical(as.character(object@market_side), "input") ||
       !identical(as.integer(object@supply_model), 3L) ||
       is.null(.native_ordered_structure(object)))) {
    stop("price_domain = 'real' currently supports only input Logit Cournot ",
         "native exact fringe or Stackelberg counterfactuals.", call. = FALSE)
  }
  price_domain
}

.require_real_rate_coordination <- function(price_domain) {
  if (identical(price_domain, "real") &&
      (!requireNamespace("coordination", quietly = TRUE) ||
       utils::packageVersion("coordination") < "0.0.0.9004")) {
    stop("price_domain = 'real' requires coordination >= 0.0.0.9004.",
         call. = FALSE)
  }
  invisible(NULL)
}

.bayes_sim_failure_detail <- function(error, stage, draw_index = NA_integer_) {
  field <- function(name, fallback = NULL) {
    value <- error[[name]]
    if (is.null(value)) fallback else value
  }
  list(stage = stage, draw_index = as.integer(draw_index),
       category = field("category", "unclassified_simulation_failure"),
       message = conditionMessage(error), condition_class = class(error),
       price_domain = field("price_domain", field("domain")),
       minimum_rate = field("minimum_rate"),
       minimum_product = field("minimum_product"),
       preMerger = field("preMerger"))
}

.model_uses_real_rate_domain <- function(model) {
  methods::is(model, "Logit") && !isTRUE(model@output) &&
    "control.equ" %in% methods::slotNames(model) &&
    identical(model@control.equ$price_domain, "real")
}
