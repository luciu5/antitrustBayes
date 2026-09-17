## Independent numerical reference for latent-quality policy propagation.
## Base R algebra is checked first.  The final block, when the antitrust
## package is installed, also checks its actual calcShares/calcMC methods.

coord_lib <- "/tmp/cga-Rlib"
if (dir.exists(file.path(coord_lib, "coordination"))) {
  .libPaths(c(coord_lib, .libPaths()))
}

logsumexp <- function(x) {
  stopifnot(length(x) > 0L, all(is.finite(x)))
  z <- max(x)
  z + log(sum(exp(x - z)))
}

softmax <- function(x) exp(x - logsumexp(x))

## x is the fitted non-s0 Logit linear predictor at the observed price.
## Structural outside utility is normalized to zero.  To preserve s0, the
## inside utility odds must sum to (1-s0)/s0, hence the minus sign below.
latent_logit <- function(x, s0, alpha, price, price_outside,
                         reanchor = TRUE) {
  stopifnot(length(x) == length(price), length(price) > 0L,
            is.finite(s0), s0 > 0, s0 < 1,
            length(alpha) == 1L, is.finite(alpha),
            length(price_outside) == 1L, is.finite(price_outside))
  c <- if (reanchor) -log(s0 / (1 - s0)) - logsumexp(x) else 0
  delta <- x + c
  t <- softmax(x)
  s <- (1 - s0) * t
  meanval <- delta - alpha * (price - price_outside)
  raw <- exp(meanval + alpha * (price - price_outside))
  denom <- 1 + sum(raw)
  list(c = c, delta = delta, t = t, s = s, s0 = 1 / denom,
       meanval = meanval, predicted = raw / denom,
       naive_s0 = 1 / (1 + sum(exp(x))))
}

## Current output CES convention: quantity odds q_i/q_0 =
## A_i * (p_i/p_0)^(-gamma).  x is the fitted non-s0 log quantity-odds
## predictor at the observed prices, including gamma*log(price_outside).
latent_ces_output <- function(x, s0, gamma, price, price_outside,
                              reanchor = TRUE) {
  stopifnot(length(x) == length(price), length(price) > 0L,
            is.finite(s0), s0 > 0, s0 < 1,
            length(gamma) == 1L, is.finite(gamma), gamma > 0,
            length(price_outside) == 1L, is.finite(price_outside),
            price > 0, price_outside > 0)
  c <- if (reanchor) -log(s0 / (1 - s0)) - logsumexp(x) else 0
  qodds <- exp(x + c)
  t <- qodds / sum(qodds)
  s <- (1 - s0) * t
  meanval <- qodds * (price / price_outside)^gamma
  revenue_den <- s0 * price_outside + sum(price * s)
  revenue <- price * s / revenue_den
  R0 <- s0 * price_outside / revenue_den
  raw_revenue <- meanval * price^(1 - gamma)
  raw_outside <- price_outside^(1 - gamma)
  revenue_check <- raw_revenue / (raw_outside + sum(raw_revenue))
  quantity_check <- (raw_revenue / price) /
    (raw_outside / price_outside + sum(raw_revenue / price))
  list(c = c, qodds = qodds, t = t, s = s, s0 = s0,
       meanval = meanval, revenue = revenue, R0 = R0,
       revenue_check = revenue_check, quantity_check = quantity_check,
       naive_s0 = 1 / (1 + sum(exp(x))))
}

## Dollar markup and marginal cost for output Logit Bertrand.  O[i,j] = 1
## when j's ownership enters product i's FOC.
logit_bertrand_output <- function(price, share, alpha, ownership) {
  stopifnot(length(price) == length(share), nrow(ownership) == length(price),
            ncol(ownership) == length(price), all(price > 0),
            all(share > 0), alpha != 0)
  conduct <- 1 - as.vector(ownership %*% share)
  markup <- 1 / (abs(alpha) * conduct)
  margin <- markup / price
  mc <- price - markup
  list(conduct = conduct, markup = markup, margin = margin, mc = mc,
       foc_residual = price - markup - mc)
}

map_panel_quality <- function(group, group_levels, b) {
  stopifnot(length(group) > 0L, length(group_levels) == length(b),
            !anyDuplicated(group_levels), all(group %in% group_levels))
  unname(b[match(group, group_levels)])
}

## --- Logit normalization and the two cost estimands ---------------------
p <- c(1.0, 1.7, 2.3)
x <- c(-0.15, 0.42, -0.75)
s0 <- 0.21
alpha <- -2.0                         # antitrust slope for output
p0 <- 0.8
lg <- latent_logit(x, s0, alpha, p, p0)
stopifnot(max(abs(sum(lg$t) - 1)) < 1e-12,
          max(abs(sum(lg$s) - (1 - s0))) < 1e-12,
          abs(lg$s0 - s0) < 1e-12,
          max(abs(lg$predicted - lg$s)) < 1e-12,
          abs(lg$naive_s0 - s0) > 1e-3)

## Design 1 keeps the observed-anchor MC while changing demand.  Design 2
## re-inverts at p, so it changes MC whenever the latent shares change.
observed_t <- c(0.50, 0.30, 0.20)
observed_s <- (1 - s0) * observed_t
O <- diag(length(p))
mc_observed <- logit_bertrand_output(p, observed_s, alpha, O)$mc
mc_latent <- logit_bertrand_output(p, lg$s, alpha, O)$mc
stopifnot(max(abs(logit_bertrand_output(p, lg$s, alpha, O)$foc_residual)) < 1e-12,
          max(abs(mc_latent - mc_observed)) > 1e-6)

## --- CES output conversion -----------------------------------------------
gamma <- 2.1
cx <- c(0.25, -0.35, 0.10)
ces <- latent_ces_output(cx, s0, gamma, p, p0)
stopifnot(max(abs(sum(ces$t) - 1)) < 1e-12,
          max(abs(sum(ces$s) - (1 - s0))) < 1e-12,
          abs(ces$R0 + sum(ces$revenue) - 1) < 1e-12,
          max(abs(ces$revenue_check - ces$revenue)) < 1e-12,
          max(abs(ces$quantity_check - ces$s)) < 1e-12)

## --- Shared panel quality mapping ----------------------------------------
groups <- c("bank_a", "bank_b", "bank_a", "bank_c", "bank_b")
levels <- c("bank_a", "bank_b", "bank_c")
b <- c(0.35, -0.20, -0.15)
panel_b <- map_panel_quality(groups, levels, b)
stopifnot(abs(sum(b)) < 1e-12,
          panel_b[1] == panel_b[3], panel_b[2] == panel_b[5])

## --- Check the actual installed antitrust share and cost methods ----------
if (requireNamespace("antitrust", quietly = TRUE)) {
  owner_post <- O[c(2, 1, 3), , drop = FALSE]
  make_common <- function(class_name, shares, slopes, share_inside) {
    methods::new(
      class_name, prices = p, margins = rep(0.2, length(p)),
      priceStart = p, normIndex = NA_real_, shareInside = share_inside,
      weights = rep(1, length(p)), priceOutside = p0, insideSize = 1,
      mktSize = 1, shares = shares, mcDelta = rep(0, length(p)),
      slopes = slopes, subset = rep(TRUE, length(p)),
      ownerPre = O, ownerPost = owner_post, pricePre = p,
      output = TRUE, labels = paste0("product_", seq_along(p))
    )
  }
  logit_obj <- make_common(
    "Logit", lg$s,
    list(alpha = alpha, meanval = lg$meanval), 1 - s0
  )
  got_logit <- antitrust::calcShares(logit_obj)
  stopifnot(max(abs(got_logit - lg$s)) < 1e-10,
            max(abs(antitrust::calcMC(logit_obj) -
                    logit_bertrand_output(p, lg$s, alpha, O)$mc)) < 1e-10)

  ces_obj <- make_common(
    "CES", ces$revenue,
    list(gamma = gamma, meanval = ces$meanval,
         alpha = ces$R0 / (1 - ces$R0)), 1 - ces$R0
  )
  got_ces_q <- antitrust::calcShares(ces_obj, revenue = FALSE)
  got_ces_r <- antitrust::calcShares(ces_obj, revenue = TRUE)
  stopifnot(max(abs(got_ces_q - ces$s)) < 1e-10,
            max(abs(got_ces_r - ces$revenue)) < 1e-10)
  cat("antitrust calcShares/calcMC checks passed\n")
} else {
  cat("antitrust not installed; package-method checks skipped\n")
}

## Coordination's native classes are optional.  These checks exercise the
## slot-level fixed-cost route used by the proposed latent bridge: construct
## once from the observed anchor, replace meanval only, solve a new pre price,
## and verify the stored MC and native FOC residuals.
if (requireNamespace("coordination", quietly = TRUE)) {
  p_native <- c(1.0, 1.7, 2.3)
  q_inside <- c(0.40, 0.24, 0.15)
  own_native <- c("a", "b", "c")
  shifts <- c(0.30, -0.20, 0.10)

  sl <- coordination::stackelberg(
    p_native, q_inside, margins = rep(NA_real_, 3),
    ownerPre = own_native, leadersPre = "a", demand = "logit",
    conduct = "bertrand", output = TRUE, insideSize = 1,
    priceOutside = 0.8, alpha = 2,
    control.equ = list(implicitCheck = FALSE)
  )
  sl_lat <- sl
  sl_lat@slopes$meanval <- sl@slopes$meanval + shifts
  sl_lat@priceStart <- p_native
  sl_lat@pricePre <- coordination::calcPrices(sl_lat, TRUE)
  stopifnot(max(abs(sl_lat@mcPre - sl@mcPre)) < 1e-12,
            max(abs(coordination::stackelberg_residuals(sl_lat, TRUE)$all)) < 1e-8)

  cf <- coordination::core_fringe(
    p_native, q_inside, margins = rep(NA_real_, 3),
    ownerPre = own_native, corePre = "a", demand = "logit",
    conduct = "bertrand", output = TRUE, insideSize = 1,
    priceOutside = 0.8, alpha = 2
  )
  cf_lat <- cf
  cf_lat@slopes$meanval <- cf@slopes$meanval + shifts
  cf_lat@priceStart <- p_native
  cf_lat@pricePre <- coordination::calcPrices(cf_lat, TRUE)
  stopifnot(max(abs(cf_lat@mcPre - cf@mcPre)) < 1e-12,
            max(abs(coordination::core_fringe_residuals(cf_lat, TRUE)$all)) < 1e-8)

  p0_native <- 0.8
  rev_inside <- p_native * q_inside /
    (p0_native * s0 + sum(p_native * q_inside))
  sl_ces <- coordination::stackelberg(
    p_native, rev_inside, margins = rep(NA_real_, 3),
    ownerPre = own_native, leadersPre = "a", demand = "ces",
    conduct = "bertrand", output = TRUE, insideSize = 1,
    priceOutside = p0_native, gamma = 2.1,
    control.equ = list(implicitCheck = FALSE)
  )
  sl_ces_lat <- sl_ces
  sl_ces_lat@slopes$meanval <- sl_ces@slopes$meanval * exp(shifts)
  sl_ces_lat@priceStart <- p_native
  sl_ces_lat@pricePre <- coordination::calcPrices(sl_ces_lat, TRUE)
  stopifnot(max(abs(sl_ces_lat@mcPre - sl_ces@mcPre)) < 1e-12,
            max(abs(coordination::stackelberg_residuals(sl_ces_lat, TRUE)$all)) < 1e-8)

  cf_ces <- coordination::core_fringe(
    p_native, rev_inside, margins = rep(NA_real_, 3),
    ownerPre = own_native, corePre = "a", demand = "ces",
    conduct = "bertrand", output = TRUE, insideSize = 1,
    priceOutside = p0_native, gamma = 2.1
  )
  cf_ces_lat <- cf_ces
  cf_ces_lat@slopes$meanval <- cf_ces@slopes$meanval * exp(shifts)
  cf_ces_lat@priceStart <- p_native
  cf_ces_lat@pricePre <- coordination::calcPrices(cf_ces_lat, TRUE)
  stopifnot(max(abs(cf_ces_lat@mcPre - cf_ces@mcPre)) < 1e-12,
            max(abs(coordination::core_fringe_residuals(cf_ces_lat, TRUE)$all)) < 1e-8)
  cat("coordination fixed-MC latent pre-solve checks passed\n")
} else {
  cat("coordination not installed; native fixed-MC checks skipped\n")
}

cat("latent-quality reference checks passed\n")
cat(sprintf("Logit reanchor c = %.8f; naive outside = %.8f; supplied outside = %.8f\n",
            lg$c, lg$naive_s0, s0))
cat(sprintf("CES outside revenue share = %.8f; latent-vs-observed MC max diff = %.8f\n",
            ces$R0, max(abs(mc_latent - mc_observed))))
