#' @import methods
#' @import rstantools
#' @aliases BMSC
#' @useDynLib BMSC, .registration = TRUE
#' @importFrom loo compare extract_log_lik loo relative_eff
#' @importFrom Rcpp loadModule
#' @importFrom dplyr %>% all_equal filter mutate_all n_distinct
#' @importFrom ggplot2 aes_string element_text geom_errorbar geom_point ggplot
#' labs scale_x_discrete theme
#' @importFrom rstan extract stan_model sampling
#' @importFrom R.utils insert
#' @importFrom stats as.formula complete.cases model.matrix na.exclude napredict 
#' terms sd quantile na.omit rnorm
globalVariables(".")
