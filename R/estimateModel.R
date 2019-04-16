#' Model selection algorithm for constrained estimation
#' @param formula formula object: formula object without exponents or
#' interactions. If \code{formula} is not of class \code{formula}, it is turned
#' into one.
#' @param mustInclude character vector: variables to include in any case; use ":" for interactions and "I(..)" for powers, e.g.: "I(x1^2):I(x2^3)".
#' @param maxExponent positive integer: highest exponent included in the
#' formula. Default is 1, e.g., only linear effects.
#' @param interactionDepth positive integer: maximum order of interaction.
#' Default is 1, e.g., only main effects (no interactions).
#' @param intercept logical: Should the intercept be included in the estimation or not?
#' @param constraint_1 logical: Should the all beta variables add up to 1?
#' @param data data.frame: dataset
#' @param yUncertainty numeric vector: optional, uncertainties in y variable
#' given in standard deviations
#' @param xUncertainty data.frame: optional, uncertainties in x variables. variable names must match with names in formula
#' @param maxNumTerms positive integer: maximum number of variables to include
#' @param scale logical: should the variables be scaled to mean 0 and sd 1?
#' @param chains positive integer: number of chains for MCMC sampling
#' @param iterations positive integer: number of iterations per chain for MCMC sampling
#' @return A list of potential models
#' @examples
#' \dontrun{
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(44)
#' n <- 80
#' x1 <- rnorm(n, sd = 1)
#' x2 <- rnorm(n, sd = 1)
#' x3 <- rnorm(n, sd = 1)
#' y <- 0.4 + 0.3 * x1 + 0.3 * x1 * x3 + 0.4 * x1 ^ 2 * x2 ^ 3 + rnorm(n, sd = 0.3)
#' yUncertainty <- rexp(n, 10) * 0.01
#' #optional (slow)
#' #xUncertainty <- data.frame(x3 = rep(0.1, n), x1 = rep(0.1, n), x2 = rep(1, n))
#' data <- data.frame(x1, x2, x3, y, yUncertainty)
#' models <- constrSelEst(y ~ x1 + x2 + x3, mustInclude = "x1", maxExponent = 3,
#'                        interactionDepth = 3, intercept = TRUE,
#'                        constraint_1 = TRUE, data = data,
#'                        yUncertainty = yUncertainty,
#'                        xUncertainty = NULL,
#'                        maxNumTerms = 10)
#' plotModelFit(models)
#' bestModel <- getBestModel(models, thresholdSE = 2)
#' print(bestModel)
#' }
#' @export
constrSelEst <- function(formula,
                         data,
                         mustInclude = "",
                         maxExponent = 1,
                         interactionDepth = 1,
                         intercept = TRUE,
                         constraint_1 = FALSE,
                         yUncertainty = rep(0, nrow(data)),
                         xUncertainty = NULL,
                         maxNumTerms = 10,
                         scale = FALSE,
                         chains = 4,
                         iterations = 2000) {
  
  withoutMissings <- handleMissingData(data, formula, yUncertainty)
  data <- withoutMissings$data
  yUncertainty <- withoutMissings$yUncertainty
  
  variableData <- selectModel(formula = formula,
                              data = data,
                              mustInclude = mustInclude,
                              maxExponent = maxExponent,
                              interactionDepth = interactionDepth,
                              intercept = intercept,
                              constraint_1 = constraint_1,
                              yUncertainty = yUncertainty,
                              scale = scale,
                              chains = chains,
                              iterations = iterations)

  models <- reAssessModels(formula = formula,
                           data = data,
                           variableData = variableData,
                           intercept = intercept,
                           constraint_1 = constraint_1,
                           yUncertainty = yUncertainty,
                           xUncertainty = xUncertainty,
                           mustInclude = mustInclude,
                           maxNumTerms = maxNumTerms,
                           scale = scale,
                           chains = chains,
                           iterations = iterations)
  return(models)
}


selectModel <- function(formula,
                        data,
                        maxExponent,
                        interactionDepth,
                        intercept,
                        constraint_1,
                        yUncertainty,
                        mustInclude,
                        scale,
                        chains,
                        iterations) {
  
  if (any(yUncertainty < 0))
    stop("uncertainties in y must be larger or equal 0")
  data$yUncertainty <- yUncertainty
  
  regFormula <- createFormula(formula, maxExponent, interactionDepth, intercept)
  mustIncludeFormula <- paste(paste(all.vars(regFormula)[1], "~",
                                    paste(mustInclude, collapse = "+")), "- 1")
  
  modelCompiled <-
    if (constraint_1) stanmodels$linRegHorseHoe else stanmodels$linRegHorseHoeUnConstr
  
  model <- estimateBayesianModel(data = data,
                                 regFormula = regFormula,
                                 mustIncludeFormula = mustIncludeFormula,
                                 intercept = intercept,
                                 modelCompiled = modelCompiled,
                                 selectVars = TRUE,
                                 scale = scale,
                                 chains = chains,
                                 iterations = iterations)
  
  X <- model.matrix(as.formula(regFormula), data = data)
  X2 <- model.matrix(as.formula(mustIncludeFormula), data = data)
  
  if (ncol(X2) > 0) {
    X <- X[, -which(duplicated(t(cbind(X, X2)), fromLast = TRUE))]
  }
  K <- unlist(lapply(1:ncol(extract(model)$lambda),
                     function(x) {
                       mean(1 - 1 / (1 + (extract(model)$lambda[, x] ^ 2 *
                                            extract(model)$tau ^ 2 /
                                            rowMeans(extract(model)$sigma ^ 2) *
                                            nrow(data))))
                     }))
  variableData <- data.frame(variable = c(colnames(X2),
                                          setdiff(colnames(X),
                                                  colnames(X)[grep("Intercept",
                                                                   colnames(X))])),
                             K = c(rep(1, ncol(X2)), round(K, 5)))
  variableData <- variableData[order(variableData$K, decreasing = TRUE), ]
  return(variableData)
}

reAssessModels <- function(formula,
                           data,
                           variableData,
                           intercept,
                           constraint_1,
                           yUncertainty,
                           xUncertainty,
                           mustInclude,
                           maxNumTerms,
                           Threshold = 0.01,
                           scale,
                           chains,
                           iterations) {
  data$yUncertainty <- yUncertainty
  
  if (constraint_1 & is.null(xUncertainty)) modelCompiled <- stanmodels$linReg 
  if (constraint_1 & !is.null(xUncertainty)) modelCompiled <- stanmodels$linRegXU
  if (!constraint_1 & is.null(xUncertainty)) modelCompiled <- stanmodels$linRegUnConstr
  if (!constraint_1 & !is.null(xUncertainty)) modelCompiled <- stanmodels$linRegUnConstrXU
  
  mustIncludeFormula <- paste(paste(all.vars(formula)[1], "~",
                                    paste(mustInclude, collapse = "+")), "- 1")
  Nterms <- max(max(1, length(mustInclude)),
                min(maxNumTerms,
                    length(variableData[variableData$K > Threshold,
                                        "variable"])))
  regFormulas <- lapply(1 : Nterms, function(x) {
    paste(all.vars(formula)[1], "~",
          paste(variableData[1 : x, "variable"], collapse = "+"))
  })

  rangeVars <- apply(model.matrix(as.formula(regFormulas[[length(regFormulas)]]),
                                    data = data)[, -1, drop = FALSE], 2, sd)
  
  if ((max(rangeVars) / min(rangeVars)) > 5000 & scale == FALSE) {
    warning("Predictor variables on very different scale. 
            Slow computation and unreliable results. 
            Please rescale your variables (or set scale = TRUE) or reduce max
            number of interactions and powers")
  }
  if (!intercept) {
    regFormulas <- lapply(1:length(regFormulas),
                          function(x) {paste(regFormulas[[x]], "- 1")})
  }
  modelList <- lapply(regFormulas, function(formula) {
    estimateBayesianModel(data = data,
                          regFormula = formula,
                          intercept = intercept,
                          xUncertainty = xUncertainty,
                          modelCompiled = modelCompiled,
                          mustIncludeFormula = mustIncludeFormula,
                          selectVars = FALSE,
                          scale = scale,
                          chains = chains,
                          iterations = iterations)
  })
  names(modelList) <- regFormulas
  modelList
}


estimateBayesianModel <- function(data,
                                  regFormula,
                                  mustIncludeFormula,
                                  intercept,
                                  xUncertainty = NULL,
                                  modelCompiled,
                                  selectVars,
                                  mc = TRUE,
                                  scale,
                                  chains,
                                  iterations) {
  N <- nrow(data)

  X <- model.matrix(as.formula(regFormula), data = data)
  
  X2 <- model.matrix(as.formula(mustIncludeFormula), data = data)
  
  if (ncol(X2) > 0) {
    X <- X[, -which(duplicated(t(cbind(X, X2)), fromLast = TRUE))]
  }
  
  if (!selectVars) {
    X <- cbind(X, X2)
  }
  
  xUnc <- 0
  xUncertaintyMatrix <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  if(!is.null(xUncertainty)){
    xUncertaintyMatrix <- getUncertaintyModelMatrix(data, xUncertainty, regFormula, X)
    xUnc <- 1
  }
  K <- ncol(X)
  y <- data[, all.vars(as.formula(regFormula))[1]]
  yUncertainty <- data$yUncertainty
  K1 <- if (intercept) 1 else 0
  K2 <- ncol(X2)
  
  scaleCenter <- rep(0, K)
  scaleScale <- rep(1, K)
  
  if (selectVars | scale == TRUE) {
    if (intercept) {
      scaledX <- scale(X[, -1])
      X[, -1] <- scaledX
      X2[, -1] <- scale(X2[, -1])
      if(!is.null(xUncertainty)){
        xUncertaintyMatrix[, -1] <- sweep(xUncertaintyMatrix[, -1], 2, attr(scaledX, "scaled:scale"), '/')
      }
      scaleCenter <- attr(scaledX, "scaled:center")
      scaleScale <- attr(scaledX, "scaled:scale")
    } else {
      scaledX <- scale(X)
      X <- scaledX
      X2 <- scale(X2)
      if(!is.null(xUncertainty)){
        xUncertaintyMatrix <- sweep(xUncertaintyMatrix, 2, attr(scaledX, "scaled:scale"), '/')
      }
      scaleCenter <- attr(scaledX, "scaled:center")
      scaleScale <- attr(scaledX, "scaled:scale")
    }
    X[is.na(X)] <- 0
    if(!is.null(xUncertainty)){
      xUncertaintyMatrix[is.na(xUncertaintyMatrix)] <- 0
    }
  }
  
  cores <- getOption("mc.cores", if (mc) chains else 1)
  #if (!selectVars) cores <- 1
  # Compute stan model
  model <- suppressWarnings(sampling(modelCompiled,
                                     chains = chains,
                                     iter = 2000,
                                     verbose = FALSE,
                                     refresh = 0,
                                     cores = cores,
                                     control = list(adapt_delta = 0.925,
                                                    max_treedepth = 14)))
  
  new("ConstrainedLinReg",
      model,
      formula = as.formula(regFormula),
      hasIntercept = intercept,
      scaleCenter = scaleCenter,
      scaleScale = scaleScale)
}


getLoo <- function(stanfit) {
  log_lik1 <- extract_log_lik(stanfit, merge_chains = FALSE)
  rel_n_eff <- relative_eff(exp(log_lik1))
  loo(log_lik1, rel_n_eff, cores = getOption("mc.cores", 4))
}


#' Get Best Model after Models Selection
#' 
#' @param models list of models fitted by \code{\link{constrSelEst}} function
#' @param thresholdSE numeric: How much standard errors in leave-one-out
#' prediction performance can the sparse model be worse than the best model
#' @param plotModels boolean: Plot models in leave-one-out evaluation plot TRUE/FALSE 
#' @return The best sparse model concerning leave-one-out performance within a
#' threshold
#' @export
getBestModel <- function(models, thresholdSE = 1, plotModels = TRUE) {
  loos <- suppressWarnings(lapply(models, getLoo))
  
  bestSparse <- bestModel(models, loos, thresholdSE)

  if(plotModels){
    print(plotModelFit(models = models,
                       thresholdSE = thresholdSE,
                       loos = loos,
                       markBestModel = TRUE))
  }
  return(models[[bestSparse]])
}

bestModel <- function(models, loos, thresholdSE){
  best <- which.max(lapply(loos, function(x) x$estimates["elpd_loo", "Estimate"]))
  
  comparison <- data.frame(do.call(rbind,
                                   unname(lapply(loos, function(x) compare(loos[[best]], x)))),
                           diffNumParam = unlist(unname(lapply(models, function(x) {
                             ncol(extract(models[[best]])$beta) - ncol(extract(x)$beta)}))))
  
  # only sparse models
  comparison$elpd_diff[comparison$diffNumParam < 0] <- -Inf
  
  which.max(comparison$elpd_diff + comparison$se * thresholdSE)
}

getUncertaintyModelMatrix <-
  function(data, xUncertainty, regFormula, X) {
    dataUnc <- data
    indices <- na.omit(match(names(xUncertainty), names(data)))
    if (length(indices) > 0) {
      modelMatrices <- lapply(1:1000, function(x) {
        dataUnc[, match(names(xUncertainty), names(data))] <-
          data[, match(names(xUncertainty), names(data))] +
          do.call("cbind", lapply(1:ncol(xUncertainty), function(x)
            rnorm(nrow(data), mean = 0, sd = xUncertainty[, x])))
        mUncertainty <-
          model.matrix(as.formula(regFormula), data = dataUnc)
      })
      xUncertaintyModelMatrix <-
        apply(simplify2array(modelMatrices), 1:2, sd)
    } else {
      xUncertaintyModelMatrix <-
        matrix(0, nrow = nrow(data), ncol = ncol(X))
    }
    return(xUncertaintyModelMatrix)
  }
