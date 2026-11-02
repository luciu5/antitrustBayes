# Deterministic checks for exact_k_reference.R. No Stan compilation or sampling.
here <- dirname(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1]))
if (!length(here) || is.na(here)) here <- "doc/stackelberg"
source(file.path(here, "exact_k_reference.R"))

firm <- c("A", "A", "B", "C", "D")
rank <- c(1L, 1L, 2L, 3L, 4L)
size <- c(.30, .30, .22, .13, .08)
price <- c(2.0, 2.8, 2.3, 1.9, 2.6)
s0 <- .27
alpha <- 1.8
nfirm <- length(unique(firm))

for (conduct in c("bertrand", "cournot")) {
  mc <- exact_k_logit_inverse_lerner(price, firm, size, rank, s0, alpha, 0,
                                      "fringe", conduct)
  stopifnot(max(abs(mc - alpha * price)) < 1e-14)
  fringe_all <- exact_k_logit_inverse_lerner(
    price, firm, size, rank, s0, alpha, nfirm, "fringe", conduct
  )
  stack_all <- exact_k_logit_inverse_lerner(
    price, firm, size, rank, s0, alpha, nfirm, "stackelberg", conduct
  )
  stopifnot(max(abs(fringe_all - stack_all)) < 1e-14)
}

q0 <- .21
p0 <- 1.3
pin <- 2.1
gamma <- 2.7
for (conduct in c("bertrand", "cournot")) {
  mc <- exact_k_ces_inverse_lerner(firm, size, rank, q0, p0, pin, gamma, 0,
                                    "fringe", conduct)
  stopifnot(max(abs(mc - gamma)) < 1e-14)
  fringe_all <- exact_k_ces_inverse_lerner(
    firm, size, rank, q0, p0, pin, gamma, nfirm, "fringe", conduct
  )
  stack_all <- exact_k_ces_inverse_lerner(
    firm, size, rank, q0, p0, pin, gamma, nfirm, "stackelberg", conduct
  )
  stopifnot(max(abs(fringe_all - stack_all)) < 1e-14)
}

# The two products of A share the same dollar markup under Logit and the same
# Lerner margin under CES.
z_logit <- exact_k_logit_inverse_lerner(price, firm, size, rank, s0, alpha, 2,
                                         "stackelberg", "bertrand")
stopifnot(abs(price[1] / z_logit[1] - price[2] / z_logit[2]) < 1e-14)
z_ces <- exact_k_ces_inverse_lerner(firm, size, rank, q0, p0, pin, gamma, 2,
                                     "stackelberg", "cournot")
stopifnot(abs(z_ces[1] - z_ces[2]) < 1e-14)

# One k is shared by all observations in a market. Compare the stable helper
# with direct finite-state summation.
ell <- matrix(c(
  -.4, -1.1, -1.8,
  -.8,  -.3, -1.2,
  -.2,  -.7,  -.5,
  -.9,  -.6,  -.4,
  -.1,  -.5, -1.0
), nrow = 5L, byrow = TRUE)
market <- c("m1", "m1", "m1", "m2", "m2")
support <- c(1L, 2L, 4L)
prior <- c(.2, .5, .3)
mix <- exact_k_market_mixture(ell, market, support, prior)
for (mm in unique(market)) {
  ii <- market == mm
  direct_score <- log(prior) + colSums(ell[ii, , drop = FALSE])
  direct <- log(sum(exp(direct_score)))
  stopifnot(abs(mix$log_mix[[mm]] - direct) < 1e-14)
}
stopifnot(max(abs(rowSums(mix$posterior) - 1)) < 1e-14)

cat("PASS: exact ordered-k R reference boundary and mixture checks.\n")
