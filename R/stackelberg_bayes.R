# Exact hard-partition support for the Bayesian noncooperative Stackelberg
# models.  This file contains no smooth conduct interpolation: every
# component is one of the n + 1 hard leader partitions.

.stack_structure_code <- c(simultaneous = 0L, fringe = 1L, stackelberg = 2L)

.stack_resolve_structure <- function(structure = NULL, use_cutoff = FALSE,
                                     use_cutoff_supplied = FALSE) {
  structure_supplied <- !is.null(structure)
  if (is.null(structure)) {
    structure <- if (isTRUE(use_cutoff)) "fringe" else "simultaneous"
  } else {
    if (length(structure) != 1L || is.na(structure)) {
      stop("'structure' must be one of 'simultaneous', 'fringe', or 'stackelberg'.")
    }
    if (is.numeric(structure) && length(structure) == 1L &&
        is.finite(structure) && structure %in% 0:2) {
      structure <- names(.stack_structure_code)[match(as.integer(structure),
                                                       .stack_structure_code)]
    }
    structure <- tolower(as.character(structure))
    aliases <- c(simultaneous = "simultaneous", fringe = "fringe",
                 stackelberg = "stackelberg", stack = "stackelberg")
    structure <- unname(aliases[structure])
    if (is.na(structure)) {
      stop("'structure' must be one of 'simultaneous', 'fringe', or 'stackelberg'.")
    }
    if (isTRUE(use_cutoff_supplied)) {
      if (isTRUE(use_cutoff) && structure == "simultaneous") {
        stop("use_cutoff = TRUE conflicts with structure = 'simultaneous'.")
      }
      if (!isTRUE(use_cutoff) && structure %in% c("fringe", "stackelberg")) {
        stop("use_cutoff = FALSE conflicts with structure = '", structure, "'.")
      }
    }
  }
  list(name = structure, code = unname(.stack_structure_code[[structure]]),
       cutoff = structure %in% c("fringe", "stackelberg"),
       supplied = structure_supplied)
}

.stack_resolve_scale <- function(structure, cutoff_scale = NULL) {
  default <- if (identical(structure, "stackelberg")) 0.01 else 0.2
  if (is.null(cutoff_scale)) return(default)
  if (length(cutoff_scale) != 1L || !is.finite(cutoff_scale) || cutoff_scale <= 0) {
    stop("'cutoff_scale' must be one positive finite prespecified scalar.")
  }
  as.numeric(cutoff_scale)
}

.stack_log_inv_logit <- function(x) {
  ## log(sigmoid(x)), stable for both tails.
  out <- numeric(length(x))
  hi <- x >= 0
  out[hi] <- -log1p(exp(-x[hi]))
  out[!hi] <- x[!hi] - log1p(exp(x[!hi]))
  out
}

.stack_log1m_inv_logit <- function(x) {
  ## log(1 - sigmoid(x)), stable for both tails.
  out <- numeric(length(x))
  hi <- x >= 0
  out[hi] <- -x[hi] - log1p(exp(-x[hi]))
  out[!hi] <- -log1p(exp(x[!hi]))
  out
}

.stack_log1mexp <- function(x) {
  ## x <= 0; preserve -Inf for exact ties.
  out <- rep(-Inf, length(x))
  ok <- is.finite(x) & x < 0
  far <- ok & x < -0.6931471805599453
  out[far] <- log1p(-exp(x[far]))
  ii <- ok & x >= -0.6931471805599453
  out[ii] <- log(-expm1(x[ii]))
  out
}

#' Exact cutoff induced hard-partition log probabilities.
#' @param sizes Decreasing conditional firm sizes.
#' @param cutoff Scalar common boundary.
#' @param scale Positive logistic transition scale.
#' @return A numeric vector for k = 0, ..., length(sizes).
#' @keywords internal
stackelberg_partition_logprob <- function(sizes, cutoff, scale = 0.01) {
  sizes <- as.numeric(sizes)
  if (!length(sizes) || any(!is.finite(sizes)) || any(diff(sizes) > 0)) {
    stop("'sizes' must be a non-empty finite non-increasing vector.")
  }
  if (length(cutoff) != 1L || !is.finite(cutoff) ||
      length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("'cutoff' must be finite and 'scale' must be positive and finite.")
  }
  n <- length(sizes)
  x <- (sizes - cutoff) / scale
  out <- numeric(n + 1L)
  out[1L] <- .stack_log1m_inv_logit(x[1L])
  out[n + 1L] <- .stack_log_inv_logit(x[n])
  if (n > 1L) {
    gap <- (sizes[-n] - sizes[-1L]) / scale
    ## For an exact tie, the k component has zero mass by definition.
    mid <- .stack_log_inv_logit(x[-n]) +
      .stack_log1m_inv_logit(x[-1L]) + .stack_log1mexp(-gap)
    mid[sizes[-1L] == sizes[-n]] <- -Inf
    out[seq.int(2L, n)] <- mid
  }
  out
}

.stack_validate_k_support <- function(k_support, structure,
                                      n_by_market = NULL) {
  if (is.null(k_support)) return(NULL)
  if (!is.numeric(k_support) || !length(k_support) || anyNA(k_support) ||
      any(!is.finite(k_support)) || any(k_support != floor(k_support)) ||
      anyDuplicated(k_support)) {
    stop("'k_support' must be a non-empty vector of unique finite integers.")
  }
  lower <- if (identical(structure, "fringe")) 0L else 1L
  if (any(k_support < lower)) {
    if (identical(structure, "fringe")) {
      stop("Fringe 'k_support' must contain non-negative integers.")
    }
    stop("Stackelberg 'k_support' must contain strictly positive integers.")
  }
  if (!is.null(n_by_market) && length(n_by_market) &&
      any(k_support > min(n_by_market))) {
    stop("Every value in 'k_support' must be feasible in every selected market.")
  }
  k_support
}

.stack_validate_public_k <- function(k_support, k_prior, structure) {
  if (identical(structure, "simultaneous")) return(invisible(NULL))
  ks <- .stack_validate_k_support(k_support, structure)
  if (!is.null(k_prior) && !is.numeric(k_prior) && !is.list(k_prior)) {
    stop("'k_prior' must be a numeric vector or a list of numeric vectors.")
  }
  if (is.numeric(k_prior)) {
    if (anyNA(k_prior) || any(!is.finite(k_prior)) || any(k_prior < 0)) {
      stop("'k_prior' must contain finite non-negative probabilities.")
    }
    if (!is.null(ks) && length(k_prior) != length(ks)) {
      stop("'k_prior' must have the same length as 'k_support'.")
    }
    if (!is.null(ks) && abs(sum(k_prior) - 1) > 1e-8) {
      stop("'k_prior' must sum to one within numerical tolerance.")
    }
  }
  invisible(NULL)
}

.stack_validate_k_prior <- function(k_prior, supports, units, k_scope,
                                    k_support = NULL) {
  ## The public form is positional: k_prior[j] belongs to k_support[j].
  ## Named count vectors/lists remain a compatibility path for old fits when
  ## no explicit support was supplied.
  nunit <- length(supports)
  if (is.null(k_prior)) {
    return(lapply(supports, function(k) {
      if (!length(k)) return(numeric())
      rep(1 / length(k), length(k))
    }))
  }
  positional <- is.numeric(k_prior) &&
    (is.null(names(k_prior)) || !length(names(k_prior)) || !is.null(k_support))
  if (is.numeric(k_prior) && positional) {
    if (anyNA(k_prior) || any(!is.finite(k_prior)) || any(k_prior < 0) ||
        abs(sum(k_prior) - 1) > 1e-8) {
      stop("'k_prior' must contain finite non-negative values and sum to one.")
    }
    if (is.null(k_support)) {
      common <- unique(unlist(supports, use.names = FALSE))
      if (length(k_prior) != length(common)) {
        stop("Positional 'k_prior' requires the same length as the common support.")
      }
      k_support <- common
    } else if (length(k_prior) != length(k_support)) {
      stop("'k_prior' must have the same length as 'k_support'.")
    }
    out <- vector("list", nunit)
    for (u in seq_len(nunit)) {
      support <- supports[[u]]
      bad <- k_support[k_prior > 0 & !(k_support %in% support)]
      if (length(bad)) {
        stop("'k_prior' assigns positive mass to unsupported leader count(s): ",
             paste(bad, collapse = ", "))
      }
      ans <- numeric(length(support)); names(ans) <- as.character(support)
      keep <- match(support, k_support)
      ans[] <- ifelse(is.na(keep), 0, k_prior[keep])
      out[[u]] <- ans
    }
    return(out)
  }
  if (is.numeric(k_prior)) {
    if (is.null(names(k_prior)) || any(!nzchar(names(k_prior)))) {
      stop("'k_prior' must be positional or a named numeric vector by count.")
    }
    choices <- rep(list(k_prior), nunit)
  } else if (is.list(k_prior)) {
    if (is.null(names(k_prior)) || any(!nzchar(names(k_prior)))) {
      stop("'k_prior' list must be named by market.")
    }
    choices <- vector("list", nunit)
    for (u in seq_len(nunit)) {
      nm <- as.character(units$market[u])
      pick <- k_prior[[nm]]
      if (is.null(pick) && !is.null(units$market_label)) {
        pick <- k_prior[[units$market_label[u]]]
      }
      if (is.null(pick)) stop("'k_prior' has no entry for market '", nm, "'.")
      choices[[u]] <- pick
    }
  } else stop("'k_prior' must be numeric or a named list by market.")
  out <- vector("list", nunit)
  for (u in seq_len(nunit)) {
    z <- choices[[u]]
    if (!is.numeric(z) || is.null(names(z)) || any(!nzchar(names(z))) ||
        anyNA(z) || any(!is.finite(z)) || any(z < 0)) {
      stop("Each legacy 'k_prior' must contain finite non-negative named probabilities.")
    }
    kk <- suppressWarnings(as.integer(names(z)))
    if (anyNA(kk) || any(names(z) != as.character(kk))) {
      stop("'k_prior' names must be integer leader counts.")
    }
    if (any(kk < 0L) || anyDuplicated(kk)) {
      stop("'k_prior' contains invalid or duplicate leader counts.")
    }
    if (abs(sum(z) - 1) > 1e-8) stop("Each 'k_prior' must sum to one.")
    support <- supports[[u]]
    bad <- kk[z > 0 & !(kk %in% support)]
    if (length(bad)) {
      stop("'k_prior' assigns positive mass to unsupported leader count(s): ",
           paste(bad, collapse = ", "))
    }
    ans <- numeric(length(support)); names(ans) <- as.character(support)
    keep <- kk %in% support
    if (any(keep)) ans[match(kk[keep], support)] <- z[keep]
    out[[u]] <- ans
  }
  out
}

.stack_unique_firm_mky <- function(mky_int, group_int, shareIn, revShare,
                                   market_int, group_levels, metric = c("logit", "ces")) {
  metric <- match.arg(metric)
  mky_int <- as.integer(mky_int); group_int <- as.integer(group_int)
  market_int <- as.integer(market_int)
  if (!length(mky_int) || length(group_int) != length(mky_int) ||
      length(shareIn) != length(mky_int) || length(revShare) != length(mky_int) ||
      length(market_int) != length(mky_int)) {
    stop("Firm-market ranking inputs must have one value per observation.")
  }
  n_mky <- max(mky_int); N <- length(mky_int); n_market <- max(market_int)
  rows <- vector("list", n_mky)
  obs_to <- integer(N)
  cursor <- 0L

  ## Rank the complete market panel once.  Each firm-period is first reduced
  ## to one aggregate share, absent firm-periods contribute zero, and periods
  ## receive equal weight.  The resulting rank is inherited by every
  ## firm-market-year row below.
  market_rows <- vector("list", n_market)
  obs_to_market <- integer(N)
  fm_cursor <- 0L
  for (mm in seq_len(n_market)) {
    mky_mm <- sort(unique(mky_int[market_int == mm]))
    gs <- sort(unique(group_int[market_int == mm]))
    period_in <- matrix(0, nrow = length(gs), ncol = length(mky_mm))
    period_rev <- matrix(0, nrow = length(gs), ncol = length(mky_mm))
    for (tt in seq_along(mky_mm)) {
      ii <- which(mky_int == mky_mm[[tt]])
      for (gg in seq_along(gs)) {
        jj <- ii[group_int[ii] == gs[[gg]]]
        if (length(jj)) {
          period_in[gg, tt] <- sum(shareIn[jj])
          period_rev[gg, tt] <- sum(revShare[jj])
        }
      }
    }
    vals_in <- rowMeans(period_in)
    vals_rev <- rowMeans(period_rev)
    metric_v <- if (metric == "logit") vals_in else vals_rev
    firm_labels <- as.character(group_levels[gs])
    ord <- order(-metric_v, firm_labels)
    gs <- gs[ord]; vals_in <- vals_in[ord]; vals_rev <- vals_rev[ord]
    n <- length(gs)
    jj <- seq.int(fm_cursor + 1L, fm_cursor + n); fm_cursor <- fm_cursor + n
    market_rows[[mm]] <- list(group = gs, shareIn = vals_in, revShare = vals_rev,
                              rank = seq_len(n), start = jj[1L], n = n,
                              market = mm)
    ii_mm <- which(market_int == mm)
    obs_to_market[ii_mm] <- jj[match(group_int[ii_mm], gs)]
  }

  for (m in seq_len(n_mky)) {
    ii <- which(mky_int == m)
    mm <- market_int[ii[[1L]]]
    mkt <- market_rows[[mm]]
    present <- mkt$group[mkt$group %in% unique(group_int[ii])]
    gs <- present[order(match(present, mkt$group))]
    vals_in <- vapply(gs, function(g) sum(shareIn[ii[group_int[ii] == g]]), numeric(1))
    vals_rev <- vapply(gs, function(g) sum(revShare[ii[group_int[ii] == g]]), numeric(1))
    n <- length(gs); jj <- seq.int(cursor + 1L, cursor + n); cursor <- cursor + n
    rows[[m]] <- list(group = gs, shareIn = vals_in, revShare = vals_rev,
                      rank = mkt$rank[match(gs, mkt$group)], start = jj[1L], n = n,
                      market = mm, mky = m)
    obs_to[ii] <- jj[match(group_int[ii], gs)]
  }
  all <- do.call(rbind, lapply(rows, function(z) data.frame(
    group = z$group, shareIn = z$shareIn, revShare = z$revShare,
    market = z$market, mky = z$mky, rank = z$rank,
    stringsAsFactors = FALSE)))
  list(
    N_firm_mky = nrow(all),
    firm_mky_shareIn = as.numeric(all$shareIn),
    firm_mky_revShare = as.numeric(all$revShare),
    firm_mky_market = as.integer(all$market),
    firm_mky_market_year = as.integer(all$mky),
    firm_mky_group = as.integer(all$group),
    firm_mky_rank = as.integer(all$rank),
    firm_mky_start = as.integer(vapply(rows, function(z) z$start, integer(1))),
    firm_mky_n = as.integer(vapply(rows, function(z) z$n, integer(1))),
    obs_to_firm_mky = as.integer(obs_to),
    firm_mky_group_label = as.character(group_levels[all$group]),
    N_firm_market = fm_cursor,
    firm_market_shareIn = as.numeric(unlist(lapply(market_rows, `[[`, "shareIn"))),
    firm_market_revShare = as.numeric(unlist(lapply(market_rows, `[[`, "revShare"))),
    firm_market_market = as.integer(unlist(lapply(market_rows, function(z) {
      rep.int(z$market, z$n)
    }))),
    firm_market_group = as.integer(unlist(lapply(market_rows, `[[`, "group"))),
    firm_market_rank = as.integer(unlist(lapply(market_rows, `[[`, "rank"))),
    firm_market_start = as.integer(vapply(market_rows, function(z) z$start, integer(1))),
    firm_market_n = as.integer(vapply(market_rows, function(z) z$n, integer(1))),
    firm_market_group_label = as.character(group_levels[
      as.integer(unlist(lapply(market_rows, `[[`, "group")))
    ]),
    obs_to_firm_market = as.integer(obs_to_market)
  )
}

.stack_partition_arrays <- function(firm, mky_to_market, market_levels,
                                    k_prior = NULL, k_scope = "market",
                                    partition_prior = "cutoff", k_support = NULL,
                                    structure = c("stackelberg", "fringe")) {
  k_scope <- match.arg(k_scope, c("market", "period"))
  partition_prior <- match.arg(partition_prior, c("cutoff", "ordered_k"))
  structure <- match.arg(structure)
  ## A continuous common boundary is a market-year object.  The ordered-k
  ## prior is the only choice for which market scope shares one count over
  ## all observed periods.
  effective_scope <- if (partition_prior == "cutoff") "period" else k_scope
  n_mky <- length(firm$firm_mky_n)
  n_market <- length(market_levels)
  if (effective_scope == "market") {
    units <- data.frame(market = seq_len(n_market), stringsAsFactors = FALSE)
    units$market_label <- market_levels
    supports <- lapply(seq_len(n_market), function(mm) {
      nmax <- firm$firm_market_n[mm]
      lower <- if (structure == "fringe") 0L else 1L
      if (is.null(k_support)) seq.int(lower, nmax) else as.integer(k_support)
    })
    mky_to_unit <- as.integer(mky_to_market)
  } else {
    units <- data.frame(market = mky_to_market, mky = seq_len(n_mky),
                        market_label = market_levels[mky_to_market],
                        stringsAsFactors = FALSE)
    supports <- lapply(seq_len(n_mky), function(mm) {
      nmax <- firm$firm_mky_n[mm]
      lower <- if (structure == "fringe") 0L else 1L
      if (is.null(k_support)) seq.int(lower, nmax) else as.integer(k_support)
    })
    mky_to_unit <- seq_len(n_mky)
  }
  if (any(vapply(supports, length, integer(1)) == 0L)) {
    stop("Each selected market must contain at least one feasible partition count.")
  }
  if (!is.null(k_support)) {
    .stack_validate_k_support(k_support, structure,
                              n_by_market = firm$firm_market_n)
  }
  probs <- if (partition_prior == "ordered_k") {
    .stack_validate_k_prior(k_prior, supports, units, k_scope, k_support)
  } else lapply(supports, function(k) rep(NA_real_, length(k)))
  comp_unit <- integer(); comp_k <- integer(); comp_prior <- numeric()
  unit_start <- integer(length(supports)); unit_count <- integer(length(supports))
  pos <- 0L
  for (u in seq_along(supports)) {
    unit_start[u] <- pos + 1L; unit_count[u] <- length(supports[[u]])
    comp_unit <- c(comp_unit, rep.int(u, length(supports[[u]])))
    comp_k <- c(comp_k, supports[[u]])
    comp_prior <- c(comp_prior, if (partition_prior == "ordered_k") log(probs[[u]]) else rep(0, length(supports[[u]])))
    pos <- pos + length(supports[[u]])
  }
  list(
    N_partition_unit = length(supports),
    partition_unit_market = as.integer(units$market),
    partition_unit_market_year = if (effective_scope == "period") as.integer(units$mky) else as.integer(vapply(seq_len(nrow(units)), function(u) {
      which(mky_to_market == units$market[u])[1L]
    }, integer(1))),
    partition_unit_start = as.integer(unit_start),
    partition_unit_count = as.integer(unit_count),
    partition_unit_k_min = as.integer(vapply(supports, min, integer(1))),
    partition_unit_k_max = as.integer(vapply(supports, max, integer(1))),
    mky_to_partition_unit = as.integer(mky_to_unit),
    N_partition_component = length(comp_k),
    partition_component_unit = as.integer(comp_unit),
    partition_component_k = as.integer(comp_k),
    partition_prior_logprob = as.numeric(comp_prior),
    partition_prior_probability = as.numeric(exp(comp_prior)),
    partition_prior = partition_prior,
    k_scope = effective_scope,
    requested_k_scope = k_scope,
    k_prior_by_unit = probs,
    partition_supports = supports,
    supplied_k_support = k_support,
    supplied_k_prior = k_prior,
    firm_market = firm
  )
}

.stack_or <- function(x, y) if (is.null(x) || !length(x)) y else x[[1L]]
.stack_whole_or <- function(x, y) if (is.null(x) || !length(x)) y else x

.stack_append_data <- function(stan_data, mky_int, market_int, group_int,
                               group_levels, shares, revShare, structure,
                               demand_model = "logit",
                               partition_prior = "cutoff", k_prior = NULL,
                               k_scope = "market", cutoff_scale = NULL,
                               k_support = NULL) {
  firm <- .stack_unique_firm_mky(mky_int, group_int, shares, revShare,
                                 market_int, group_levels,
                                 metric = if (identical(demand_model, "ces")) "ces" else "logit")
  p <- .stack_partition_arrays(firm, stan_data$mky_to_market,
                               attr(stan_data, "market_levels"), k_prior,
                               k_scope, partition_prior, k_support,
                               structure = .stack_or(structure$name, structure))
  ## Return only names declared in the Stackelberg Stan data block.  The
  ## descriptive labels and nested prior metadata stay in fit metadata; Stan
  ## rejects unknown list elements during data validation.
  array1 <- function(x) array(x, dim = length(x))
  market_obs <- lapply(seq_len(stan_data$N_market), function(mm) {
    which(as.integer(market_int) == mm)
  })
  market_mky <- lapply(seq_len(stan_data$N_market), function(mm) {
    which(as.integer(stan_data$mky_to_market) == mm)
  })
  market_obs_n <- lengths(market_obs)
  market_mky_n <- lengths(market_mky)
  out <- c(
    list(
      N_firm_mky = firm$N_firm_mky,
      firm_mky_shareIn = array1(firm$firm_mky_shareIn),
      firm_mky_revShare = array1(firm$firm_mky_revShare),
      firm_mky_market = array1(firm$firm_mky_market),
      firm_mky_market_year = array1(firm$firm_mky_market_year),
      firm_mky_group = array1(firm$firm_mky_group),
      firm_mky_rank = array1(firm$firm_mky_rank),
      firm_mky_start = array1(firm$firm_mky_start),
      firm_mky_n = array1(firm$firm_mky_n),
      obs_to_firm_mky = array1(firm$obs_to_firm_mky),
      N_firm_market = firm$N_firm_market,
      firm_market_shareIn = array1(firm$firm_market_shareIn),
      firm_market_revShare = array1(firm$firm_market_revShare),
      firm_market_market = array1(firm$firm_market_market),
      firm_market_group = array1(firm$firm_market_group),
      firm_market_rank = array1(firm$firm_market_rank),
      firm_market_start = array1(firm$firm_market_start),
      firm_market_n = array1(firm$firm_market_n),
      obs_to_firm_market = array1(firm$obs_to_firm_market),
      market_obs_order = array1(unlist(market_obs, use.names = FALSE)),
      market_obs_start = array1(cumsum(c(1L, utils::head(market_obs_n, -1L)))),
      market_obs_n = array1(market_obs_n),
      market_mky_order = array1(unlist(market_mky, use.names = FALSE)),
      market_mky_start = array1(cumsum(c(1L, utils::head(market_mky_n, -1L)))),
      market_mky_n = array1(market_mky_n),
      N_partition_unit = p$N_partition_unit,
      partition_unit_market = array1(p$partition_unit_market),
      partition_unit_market_year = array1(p$partition_unit_market_year),
      partition_unit_start = array1(p$partition_unit_start),
      partition_unit_count = array1(p$partition_unit_count),
      partition_unit_k_min = array1(p$partition_unit_k_min),
      partition_unit_k_max = array1(p$partition_unit_k_max),
      mky_to_partition_unit = array1(p$mky_to_partition_unit),
      N_partition_component = p$N_partition_component,
      partition_component_unit = array1(p$partition_component_unit),
      partition_component_k = array1(p$partition_component_k),
      partition_prior_logprob = array1(p$partition_prior_logprob),
      partition_prior_probability = array1(p$partition_prior_probability),
      structure = as.integer(
        .stack_structure_code[[.stack_or(structure$name, structure)]]
      ),
      partition_prior_code = as.integer(identical(partition_prior, "ordered_k")),
      k_scope_code = as.integer(identical(p$k_scope, "market")),
      cutoff_scale = .stack_resolve_scale(.stack_or(structure$name, structure), cutoff_scale)
    )
  )
  attr(out, "firm_market_ranking") <- list(
    rule = "equal_period_mean_after_firm_period_aggregation",
    metric = if (identical(demand_model, "ces")) {
      "conditional_inside_revenue_share"
    } else {
      "conditional_inside_quantity_share"
    },
    tie_break = "descending_share_then_ascending_firm_label",
    N_firm_market = firm$N_firm_market,
    firm_market_shareIn = firm$firm_market_shareIn,
    firm_market_revShare = firm$firm_market_revShare,
    firm_market_market = firm$firm_market_market,
    firm_market_group = firm$firm_market_group,
    firm_market_rank = firm$firm_market_rank,
    firm_market_group_label = firm$firm_market_group_label,
    obs_to_firm_market = firm$obs_to_firm_market
  )
  attr(out, "partition_metadata") <- p
  out
}

## Build the same fixed ranking and exact partition metadata for a one-market
## bayes_specify() object.  Specified fits have no Stan data list, so this
## information lives in their fit metadata instead.
.stack_specified_metadata <- function(params, market) {
  exact <- identical(params$partition_prior, "ordered_k") &&
    params$structure %in% c("fringe", "stackelberg")
  if (!exact) return(list(ranking = NULL, partition = NULL))
  labels <- sort(unique(as.character(market$group)))
  group_int <- match(as.character(market$group), labels)
  q <- as.numeric(market$quantities)
  share <- q / sum(q)
  revenue <- as.numeric(market$prices) * q
  rev_share <- revenue / sum(revenue)
  firm <- .stack_unique_firm_mky(
    mky_int = rep.int(1L, market$n_prod), group_int = group_int,
    shareIn = share, revShare = rev_share,
    market_int = rep.int(1L, market$n_prod), group_levels = labels,
    metric = if (identical(params$demand_model, "ces")) "ces" else "logit"
  )
  p <- .stack_partition_arrays(
    firm, mky_to_market = 1L, market_levels = "unit_1",
    k_prior = params$k_prior, k_scope = "market",
    partition_prior = params$partition_prior,
    k_support = params$k_support, structure = params$structure
  )
  list(
    ranking = list(
      rule = "equal_period_mean_after_firm_period_aggregation",
      metric = if (identical(params$demand_model, "ces")) {
        "conditional_inside_revenue_share"
      } else {
        "conditional_inside_quantity_share"
      },
      tie_break = "descending_share_then_ascending_firm_label",
      N_firm_market = firm$N_firm_market,
      firm_market_shareIn = firm$firm_market_shareIn,
      firm_market_revShare = firm$firm_market_revShare,
      firm_market_market = firm$firm_market_market,
      firm_market_group = firm$firm_market_group,
      firm_market_rank = firm$firm_market_rank,
      firm_market_group_label = firm$firm_market_group_label,
      obs_to_firm_market = firm$obs_to_firm_market
    ),
    partition = p
  )
}
