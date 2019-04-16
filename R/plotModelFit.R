#' Plot errors of all models
#' @description This plot is automatically produced with the execution of
#' \code{\link{getBestModel}}.
#' @param models List with models of class \code{\link{ConstrainedLinReg}}
#' @param loos List with the model fit results for all models as returned by
#' \code{BMSC:::getLoo}. If not provided, they are
#' computed from the model list, which can take some time.
#' @param thresholdSE numeric: Factor multiplied with standard error to obtain
#' ends of error bars
#' @param markBestModel boolean: highlight position of the best model in the model list
#' @export
plotModelFit <- function(models,
                         thresholdSE = 1,
                         loos = NULL,
                         markBestModel = TRUE) {
  if (is.null(loos)) {
    loos <- suppressWarnings(lapply(models, getLoo))
  }
  if (markBestModel){
    bestSparse <- bestModel(models, loos, thresholdSE)
  }
  
  modelNames <- prepModelNames(models)
  
  datPlot <- prepPlotData(loos = loos,
                          thresholdSE = thresholdSE,
                          modelNames = modelNames)
  
  colours <- if (!(markBestModel)) "black" else prepColorVec(bestSparse, length(models))
  
  plotModels(datPlot, colours, thresholdSE)
}


#' Extract model names from model objects
#' @description Extracts the model formulae from a list of model objects of
#' class \code{\link{ConstrainedLinReg}}. Elements that are superfluous for
#' reading (e.g., brackets) are removed.
#' @inheritParams plotModelFit
prepModelNames <- function(models) {
  lapply(models, function(x) as.character(x@formula)[3]) %>%
    unlist %>% 
    gsub("I\\(|\\)", "", .)
}


#' Prepare data to plot model fit
#' @inheritParams plotModelFit
#' @param modelNames Names for the models in the same order as they appear in
#' \code{loos}
#' @return A data.frame with the columns \code{Estimate} (Estimate of the
#' looic), \code{SE}, \code{model}, \code{lower}, and \code{upper}
prepPlotData <- function(loos, modelNames, thresholdSE) {
  datPlot <- lapply(loos, function(x) {
    res <- x$estimates["looic", c("Estimate", "SE")]
  }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
  datPlot$model <- modelNames
  datPlot$lower <- datPlot$Estimate - thresholdSE * datPlot$SE
  datPlot$upper <- datPlot$Estimate + thresholdSE * datPlot$SE
  datPlot
}


#' Prepare colour vector
#' @inheritParams plotModelFit
#' @param posBestModel numeric: position of best Model
#' @param length numeric: Length of colour vector
#' @return Vector of length \code{length}. It contains "black" expect for the
#' position provided in \code{posBestModel}, which is "chartreuse4" (green)
prepColorVec <- function(posBestModel, length) {
  colours <- rep("black", length)
  colours[posBestModel] <- "chartreuse4"
  colours
}


#' Plot model errors with errorbars
#' @param datPlot data.frame with prepared plot data
#' @param colours character: colour(s) for the points, bars and x-axis labels
#' @inheritParams plotModelFit
plotModels <- function(datPlot, colours, thresholdSE) {
  ggplot(datPlot, aes_string(x = "model")) +
    geom_point(aes_string(y = "Estimate"), colour = colours) +
    geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), colour = colours) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = colours, size = 6.5)) +
    labs(x = NULL,
         y = "Error",
         title = paste0("Fit of all models (Error bars show the Looic +/- ",
                       thresholdSE,
                       "*SE)")) +
    scale_x_discrete(labels = function(x) parse(text = x))
}
