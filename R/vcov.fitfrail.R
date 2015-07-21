#' 
#' Compute the variance/covariance matrix for a fitted object
#' 
#' @export
vcov.fitfrail <- function(fit, boot=TRUE, B=100) {
  
  v <- rexp(fit$n_clusters)
  v <- v/mean(v)
}