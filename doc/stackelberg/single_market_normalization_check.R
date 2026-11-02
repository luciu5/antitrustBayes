#!/usr/bin/env Rscript
## Independent algebra/numerical check of the proposed Bayes single-market
## normalization.  No repository code is imported or modified.
set.seed(913)

p <- c(2.1, 2.7, 3.4, 4.2, 5.1)
qc <- c(.08, .17, .23, .29, .23) # conditional inside shares, sums to one
p0 <- 1.3
meanp <- mean(p)
meanlp <- mean(log(p))
firms <- c("A", "A", "B", "C", "C")
mean_firm <- function(x) mean(tapply(x, firms, mean))

## Existing single-market predictor is logit(s0) + b_quality plus a centered
## price term.  b_quality is centered over unique firms by b_group_raw-mean().
check_logit <- function(a, sign) {
  ## sign=-1 output; sign=+1 input. This is the signed market coefficient.
  dnorm <- log(qc) - sign * a * (p - p0)
  ## Centering gives eta0 = mean_firm(dnorm) - sign*a*(p0-meanp).
  ## For output sign=-1 this is the audit's stated +a*(p0-meanp).
  eta0 <- mean_firm(dnorm) - sign * a * (p0 - meanp)
  s0 <- plogis(eta0)
  bq <- dnorm - mean_firm(dnorm)
  pred <- eta0 + bq + sign * a * (p - meanp)
  ## The output-side formula in the audit is sign=-1, hence its plus a term.
  c(sign = sign, implied_s0 = s0, max_log_share_error = max(abs(pred - log(qc))),
    centered_bq_mean = mean_firm(bq))
}

check_ces <- function(gamma, outside_price = p0, price_scale = 1) {
  pp <- price_scale * p
  pp0 <- price_scale * outside_price
  mean_log_pp <- mean(log(pp))
  dnorm <- log(qc) + gamma * log(pp / pp0)
  eta0 <- mean_firm(dnorm) - gamma * mean_log_pp
  s0 <- plogis(eta0)
  bq <- dnorm - mean_firm(dnorm)
  pred <- eta0 + bq - gamma * (log(pp) - mean_log_pp) + gamma * log(pp0)
  c(gamma = gamma, price_scale = price_scale, implied_s0 = s0,
    max_log_share_error = max(abs(pred - log(qc))),
    centered_bq_mean = mean_firm(bq))
}

cat("Logit output/input exact identities:\n")
print(rbind(check_logit(0.8, -1), check_logit(0.8, 1)), digits = 12)
cat("CES exact identity and outside-price scaling:\n")
print(rbind(check_ces(2.4, 1.3), check_ces(2.4, 1.3, 10)), digits = 12)

## The centered single-market model cannot generally retain an arbitrary s0:
## centering fixes logit(s0) from q, p and the signed slope.  This example
## yields a value different from .30.
implied <- check_logit(.8, -1)[["implied_s0"]]
cat(sprintf("Output Logit implied s0=%.12f versus arbitrary s0=0.30\n", implied))
stopifnot(max(abs(check_logit(.8, -1)[c("max_log_share_error", "centered_bq_mean")])) < 1e-12,
          max(abs(check_logit(.8, 1)[c("max_log_share_error", "centered_bq_mean")])) < 1e-12,
          max(abs(check_ces(2.4, 1.3)[c("max_log_share_error", "centered_bq_mean")])) < 1e-12,
          abs(implied - .30) > 1e-3,
          max(abs(check_ces(2.4, 1.3, 10)[c("max_log_share_error", "centered_bq_mean")])) < 1e-12,
          abs(check_ces(2.4, 1.3)[["implied_s0"]] -
              check_ces(2.4, 1.3, 10)[["implied_s0"]]) > 1e-3)
cat("PASS: single-market normalization identities and unit-scaling caveat.\n")
