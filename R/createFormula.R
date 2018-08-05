#' Create a formula with interactions and polynomials up to a desired order
#' 
#' @description Creates a formula with interactions and polynomials up to a
#' desired order. If the input \code{formula} already includes interactions,
#' exponents or other functions (e.g., \code{\link[base]{sqrt}}), they are
#' ignored.
#' 
#' @param formula formula object: formula object without exponents or
#' interactions. If \code{formula} is not of class \code{formula}, it is turned
#' into one.
#' @param maxExponent positive integer: highest exponent included in the
#' formula. Default is 1, e.g., only linear effects.
#' @param interactionDepth positive integer: maximum order of interaction.
#' Default is 1, e.g., only main effects (no interactions).
#' @param intercept logical: include intercept or not?
#' 
#' @return A formula containing the original independent variables and their
#' polynomials and interactions.
#' 
#' @examples
#' createFormula("y ~ x1 + x2", 2, 3)
#' createFormula(as.formula("y ~ x1 + x2"), interactionDepth = 2)
#' 
#' carFormula <- createFormula("mpg ~ cyl + disp + drat", 2, 3)
#' summary(lm(carFormula, mtcars))
#' @export
createFormula <- function(formula, maxExponent = 1,
                          interactionDepth = 1, intercept = TRUE) {
  
  if (maxExponent < 1 | interactionDepth < 1)
    stop("maxExponent and interactionDepth must be positive integers.")
  
  formula <- tryAsFormula(formula)
  
  if (length(formula) < 3)
    stop("Formula is not complete.")
  
  allVars <- all.vars(formula)
  
  if (length(allVars) == 1) { # Intercept only - nothing to do
    return(formula)
  } else {
    return(createFormulaInternal(formula, allVars, maxExponent,
                                 interactionDepth, intercept))
  }
}


#' Create formula with interactions and polynomials if all checks in
#' \code{\link{createFormula}} have passed
#' @param formula formula object
#' @param allVars object returned by \code{\link[base]{all.vars}}
#' @param maxExponent positive integer
#' @param interactionDepth positive integer
#' @param intercept boolean
createFormulaInternal <- function(formula, allVars, maxExponent,
                                  interactionDepth, intercept) {
  res <- makePoly(allVars[-1], maxExponent)
  if (interactionDepth > 1) res <- makeInteractions(res, interactionDepth)
  res <- paste(res, collapse = " + ")
  if (attr(terms(formula), "intercept") == 0 | intercept == FALSE) {
    res <- paste(res, "- 1") # No intercept
    }
  return(as.formula(paste(allVars[1], res, sep = " ~ ")))
}


#' Turn character vector into formula, return error if not possible
#' @param input character
#' @return Formula or error
tryAsFormula <- function(input) {
  tryCatch(as.formula(input),
           error = function(e) {
             stop(paste0('"', input, '" cannot be turned into a formula.'))
           })
}


#' Create polynomial of degree \code{maxExponent} from variable names
#' @description Remark: Since this function is to be used only within
#' \code{\link{createFormula}}, the validity of the input is not checked here
#' but in \code{\link{createFormula}}.
#' @param vars character: variable names
#' @param maxExponent integer: highest exponent
#' @return Character vector of \code{length(vars)} times \code{maxExponent}
#' @examples BMSC:::makePoly(vars = c("x1", "x2"), maxExponent = 3)
makePoly <- function(vars, maxExponent) {
  lapply(1:maxExponent, addPowToVars, vars = vars) %>% 
    unlist
}


#' Add exponent to a vector of variables
#' @description Remark: Since this function is to be used only within
#' \code{\link{createFormula}}, the validity of the input is not checked here
#' but in \code{\link{createFormula}}.
#' @param vars character: variable names
#' @param power integer: exponent
#' @return Vector of same length as \code{vars}
#' @examples BMSC:::addPowToVars(c("x1", "x2"), 2)
addPowToVars <- function(vars, power) {
  if (power == 1) vars else paste0("I(", vars, "^", power, ")")
}


#' Add all interactions up to a desired order
#' @details Interactions of variables with themselves (including polynomials
#' of themselves) are not included.
#' @param vars character: variable names (potentially including polynomial
#' expressions)
#' @param interactionDepth integer: highest interaction order
#' @return Character vector
#' @examples BMSC:::makeInteractions(vars = c("x1", "x2",
#' "I(x1^2)", "I(x2^2)"), interactionDepth = 3)
makeInteractions <- function(vars, interactionDepth) {
  lapply(1:interactionDepth, addInteractionToVars, vars = vars) %>% 
    unlist
}


#' Add interactions of a specific order to a vector of variables
#' @details Interactions of variables with themselves (including polynomials
#' of themselves) are not included.
#' @param vars character: variables
#' @param order integer: order of the interaction
#' @return Character vector
#' @examples BMSC:::addInteractionToVars(3, c("x1", "x2", "x3"))
addInteractionToVars <- function(order, vars) {
  df <- expand.grid(rep(list(vars), order), stringsAsFactors = FALSE)
  df <- df[, rev(names(df)), drop = FALSE]
  
  # Filter out all rows where at least one variable appears repeatedly
  nDistinctVars <- mutate_all(df, extractVarname) %>% apply(1, n_distinct)
  df <- filter(df, nDistinctVars == order)
  
  # Paste together variables for interactions, filter out duplicate cases
  apply(df, 1, sortAndPaste) %>% unique
}


#' Sort a vector and collapse elements together using ":"
#' @param x Vector
#' @examples BMSC:::sortAndPaste(c("var1", "var2"))
sortAndPaste <- function(x) {
  paste(sort(x), collapse = ":")
}


#' Extract variable name from polynomial expression
#' @param x Character: variables
#' @examples BMSC:::extractVarname(c("x1",
#' "I(x2^2)"))
extractVarname <- function(x) {
  x %>%
    gsub("^I\\(", "", .) %>%
    gsub("\\^[0-9]+\\)", "", .)
}
