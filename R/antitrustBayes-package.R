#' Bayesian Merger Simulation for Antitrust Analysis
#'
#' Bayesian estimation of differentiated products demand and supply
#' models using Stan, with merger simulation via the \pkg{antitrust}
#' package. Provides posterior distributions of merger price effects,
#' consumer welfare changes, and other competition metrics.
#'
#' The package provides two main functions:
#' \describe{
#'   \item{\code{\link{bayes_est}}}{Estimate a Bayesian supply-and-demand
#'     model using Stan MCMC sampling.}
#'   \item{\code{\link{bayes_sim}}}{Simulate merger effects using posterior
#'     draws from \code{bayes_est} and the \pkg{antitrust} package
#'     infrastructure.}
#' }
#'
#' @name antitrustBayes-package
#' @aliases antitrustBayes
#' @keywords internal
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling extract
#' @importFrom rstan summary
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom bridgesampling bridge_sampler
#' @importFrom digest digest
#' @importFrom loo loo waic loo_compare loo_model_weights extract_log_lik relative_eff
#' @importFrom antitrust calcMC calcMargins calcPrices calcShares ownerToMatrix elast diversion CV HypoMonTest
#' @importClassesFrom antitrust AntitrustFit StructuralFit Logit CES
#' @importFrom stats plogis quantile sd weighted.mean
#' @importFrom utils capture.output
#' @useDynLib antitrustBayes, .registration = TRUE
"_PACKAGE"
