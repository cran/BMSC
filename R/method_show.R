#' Print constraint estimation model
#' @description Prints the model formula and estimates as well as sigma with
#' the corresponding 95% credibility intervals.
#' @param object Model of class "MPIconstraintModel"
#' @export
setMethod("show",
          "ConstrainedLinReg",
          function(object) {
            
            betaMatrix <- getBetaMatrix(object, object@hasIntercept)
            sigmas <- extract(object)$sigma
            log_lik <- extract(object)$log_lik %>% rowSums %>% mean()
            rsq <- extract(object)$rsq %>% mean
            
            coefTable <- data.frame(
              Estimate = round(colMeans(betaMatrix), 3),
              Cred_Interval_0.95 =
                paste0("[",
                       extractQuantile(betaMatrix, 0.025),
                       ", ",
                       extractQuantile(betaMatrix, 0.975),
                       "]")
            )
            
            cat(displayFormula(object))
            print(coefTable)
            cat("\n")
            cat(displaySigmas(sigmas))
            cat("\n")
            cat(displayRsq(rsq))
            cat("\n")
            cat(displayLogLik(log_lik))
          }
)


#' Print constraint estimation model
#' @param x model object of class \code{\link{ConstrainedLinReg}}
#' @param ... arguments passed from or to other methods
#' @export
print.ConstrainedLinReg <- function(x, ...) {
  show(x)
}


#' Extract beta matrix from \code{\link{ConstrainedLinReg}} model
#' @description Extracts matrix of beta estimates
#' @param model model object: Model of class \code{\link{ConstrainedLinReg}}
#' @param hasIntercept logical: Does the model formula include an intercept?
#' @return matrix of estimates
getBetaMatrix <- function(model, hasIntercept) {
  
  betaMatrix <- extract(model)$beta
  colnames(betaMatrix) <- attr(terms(as.formula(attributes(model)$formula)),
                               'term.labels')
  
  if (hasIntercept) {
    betaMatrix <- do.call(cbind, list(extract(model)$beta0, betaMatrix))
    colnames(betaMatrix)[1] <- "Intercept"
  }
  
  return(betaMatrix)
}


extractQuantile <- function(betaMatrix, quant, digits = 3) {
  apply(betaMatrix,
        MARGIN = 2,
        function(x) round(quantile(x, quant), digits))
}


displayFormula <- function(modelObj) {
  paste("Model formula:",
        Reduce(paste, deparse(modelObj@formula)),
        "\n\n")
}


displaySigmas <- function(sigmas) {
  paste0("Sigma: ",
         round(mean(sigmas), 3),
         " (Cred_Interval_0.95: [",
         round(quantile(sigmas, 0.025), 3),
         ", ",
         round(quantile(sigmas, 0.975), 3),
         "])")
}

displayRsq <- function(rsq) {
  paste("R-squared:",
        round(rsq, 3))
}

displayLogLik <- function(log_lik) {
  paste("Log-likelihood:",
        round(log_lik, 3))
}
