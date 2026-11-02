#' Bayesian Counterfactual Simulation for Antitrust Analysis
#'
#' Bayesian estimation of differentiated products demand and supply
#' models using Stan, with counterfactual simulation via the \pkg{antitrust}
#' package. Provides posterior distributions of price effects,
#' consumer welfare changes, and other competition metrics.
#'
#' The package provides three main entry points:
#' \describe{
#'   \item{\code{\link{bayes_est}}}{Estimate a Bayesian supply-and-demand
#'     model using Stan MCMC sampling.}
#'   \item{\code{\link{bayes_specify}}}{Construct a Bayesian structural fit
#'     from directly supplied parameter draws.}
#'   \item{\code{\link[stats]{simulate}}}{Apply an \pkg{antitrust}
#'     counterfactual across posterior draws. \code{\link{bayes_sim}} is a
#'     convenience wrapper around this S4 method.}
#' }
#'
#' @name antitrustBayes-package
#' @aliases antitrustBayes
#' @keywords internal
#'
#' @import methods
#' @importFrom bridgesampling bridge_sampler
#' @importFrom digest digest
#' @importFrom loo loo waic loo_compare loo_model_weights relative_eff
#' @importFrom antitrust calcMC calcMargins calcPrices calcShares ownerToMatrix elast diversion CV HypoMonTest
#' @importClassesFrom antitrust AntitrustFit StructuralFit Logit CES
#' @importMethodsFrom antitrust simulate
#' @importFrom stats plogis quantile sd weighted.mean
#' @importFrom utils capture.output
"_PACKAGE"
