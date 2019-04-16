#' S4 class for constrained linear regression models
#' @description Inherits from \code{\link[rstan]{stanfit}}
#' @slot formula model formula (class formula)
#' @slot hasIntercept logical: Does the model formula include an intercept?
#' @slot scaleCenter numeric: location scale of betas 
#' @slot scaleScale numeric: scale scale of betas 
#' @export
ConstrainedLinReg <- setClass("ConstrainedLinReg",
                              slots = list(formula = "formula",
                                           hasIntercept = "logical",
                                           scaleCenter = "numeric",
                                           scaleScale = "numeric"),
                              contains = "stanfit")
