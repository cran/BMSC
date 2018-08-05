#' S4 class for constrained linear regression models
#' @description Inherits from \code{\link[rstan]{stanfit}}
#' @slot formula model formula (class formula)
#' @slot hasIntercept logical: Does the model formula include an intercept?
#' @export
ConstrainedLinReg <- setClass("ConstrainedLinReg",
                              slots = list(formula = "formula",
                                           hasIntercept = "logical"),
                              contains = "stanfit")
