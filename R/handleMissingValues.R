#' Exclude rows with missing values
#' @description All rows with missing values on the variables from the model
#' formula are excluded. If all rows are excluded, an error occurs. If only some
#' of the rows are excluded, the number and percentage of excluded rows is
#' printed via a message. In addition, the corresponding positions from
#' the yUncertainty vector are excluded.
#' @param data data.frame
#' @param formula formula object
#' @param yUncertainty numeric: vector
#' @return A list with the elements "data" (data frame containing only the
#' relevant variables and complete rows) and "yUncertainty".
handleMissingData <- function(data, formula, yUncertainty) {
  
  relevantData <- data[, names(data) %in% all.vars(formula), drop = FALSE]
  completeRows <- complete.cases(relevantData)
  relevantData <- relevantData[completeRows, ]
  yUncertainty <- yUncertainty[completeRows]
  
  if (nrow(relevantData) == 0) {
    stop("There are no rows without missing values in the data.")
  }
  
  nExcludedRows <- nrow(data) - nrow(relevantData)
  
  if (nExcludedRows > 0) {
    pctExcludedRows <- round(100 * nExcludedRows / nrow(data), 2)
    message(paste0(nExcludedRows, " rows (", pctExcludedRows,
                   "% of the data) were excluded due to missing values."))
  }
  
  return(list(data = relevantData, yUncertainty = yUncertainty))
}
