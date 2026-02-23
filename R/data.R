#' Philadelphia National Bank Merger Data
#'
#' Data from the Philadelphia National Bank merger case (1960), one of the
#' landmark antitrust cases in U.S. banking. Contains deposit-market data
#' for 39 banks in the Philadelphia metropolitan area. The proposed merger
#' was between Philadelphia National Bank (group \code{"298816"}) and Girard
#' Trust Corn Exchange Bank (group \code{"889018"}).
#'
#' @format A data frame with 39 rows and 5 columns:
#' \describe{
#'   \item{bank}{Character. Bank name.}
#'   \item{group}{Character. Holding company identifier (top-holder RSSD ID).}
#'   \item{prices}{Numeric. Interest rate on deposits (percent).}
#'   \item{quantities}{Numeric. Total deposits (thousands of dollars).}
#'   \item{margins}{Numeric. Lerner index (p - c) / p.}
#' }
#'
#' @details
#' The merging parties are:
#' \itemize{
#'   \item \strong{Acquirer:} Philadelphia National Bank (group \code{"298816"},
#'     ~21\% share)
#'   \item \strong{Target:} Girard Trust Corn Exchange Bank (group
#'     \code{"889018"}, ~14\% share)
#' }
#'
#' @source Derived from FDIC Summary of Deposits and Call Report data for
#'   the Philadelphia banking market, circa 1960. See Taragin (2026) for
#'   details.
#'
#' @examples
#' data(pnb)
#' head(pnb)
#' sum(pnb$quantities)  # total deposits
#'
#' @name pnb
#' @docType data
#' @keywords datasets
NULL
