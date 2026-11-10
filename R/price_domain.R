## Deterministic simulation domain; this never modifies a fitted posterior.
.validate_simulation_price_domain <- function(object, price_domain) {
  price_domain <- match.arg(price_domain, c("positive", "real"))
  if (identical(price_domain, "real")) {
    supply <- as.integer(object@supply_model)
    native <- .native_ordered_structure(object)
    structure_name <- object@metadata$structure
    if (is.null(structure_name)) structure_name <- "simultaneous"
    eligible <- identical(as.character(object@demand_model), "logit") &&
      identical(as.character(object@market_side), "input") &&
      ((is.null(native) && identical(structure_name, "simultaneous") &&
        supply %in% 1:5) ||
       (!is.null(native) && supply %in% c(1L, 3L)))
    if (!eligible) {
      stop("price_domain = 'real' supports simultaneous input Logit games or ",
           "native exact Bertrand/Cournot fringe and Stackelberg games.",
           call. = FALSE)
    }
  }
  price_domain
}

.require_real_rate_coordination <- function(price_domain, supply_model = NULL) {
  minimum <- if (identical(as.integer(supply_model), 1L))
    "0.0.0.9006" else "0.0.0.9004"
  if (identical(price_domain, "real") &&
      (!requireNamespace("coordination", quietly = TRUE) ||
       utils::packageVersion("coordination") < minimum)) {
    stop("price_domain = 'real' requires coordination >= ", minimum, ".",
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

.check_atomic_price_domain <- function(result, price_domain) {
  states <- if (methods::is(result, "CounterfactualPath")) result@results else
    list(result)
  for (model in states) {
    for (pre in c(TRUE, FALSE)) {
      prices <- as.numeric(if (pre) model@pricePre else model@pricePost)
      active <- if (pre) attr(model, "antitrust_subset_pre", exact = TRUE) else
        as.logical(model@subset)
      if (is.null(active)) active <- rep(TRUE, length(prices))
      if (length(active) != length(prices) || anyNA(active) || !any(active)) {
        stop("Simulated input Logit model has an invalid active-product map.",
             call. = FALSE)
      }
      bad <- which(active & (!is.finite(prices) |
        (identical(price_domain, "positive") & prices <= 0)))
      if (!length(bad)) next
      i <- bad[[1L]]
      finite <- prices[active & is.finite(prices)]
      minimum <- if (length(finite)) min(finite) else NA_real_
      category <- if (is.finite(prices[[i]]))
        "positive_domain_violation" else "invalid_price_solution"
      stop(structure(list(
        message = if (identical(category, "positive_domain_violation"))
          "Input Logit equilibrium has a nonpositive rate; use price_domain = 'real' to retain signed rates."
        else "Input Logit equilibrium has a nonfinite active rate.",
        call = NULL, category = category, price_domain = price_domain,
        minimum_rate = minimum,
        minimum_product = as.character(model@labels[[i]]),
        preMerger = pre
      ), class = c("bayes_price_domain_error", "error", "condition")))
    }
  }
  invisible(result)
}
