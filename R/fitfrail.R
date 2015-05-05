
#' Fit a Shared Frailty model
#'
#' \code{fitfrail} fits a model.
#'
#' This function ...
#' 
#' @param ... model
#' @return parameters
#'   \url{} for more details.
#' @examples
#' fitfrail()
#' 
#' @export
fitfrail <- function(data) {
  
  # The time steps will be used throughout this function
  time_steps <- sort(unique(data$time))
  
  # The number of individuals in family i that failed up to time t
  N <- tapply(data$time, data$family, function(x) {
    colSums(outer(unlist(x), time_steps, function(x,y) {as.numeric(x>=y)}))
  })
}