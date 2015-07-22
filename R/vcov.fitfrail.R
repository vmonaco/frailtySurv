#' 
#' Compute the variance/covariance matrix for a fitted object
#' 
#' @export
vcov.fitfrail <- function(fit, boot=TRUE, B=100, 
                          include.Lambda=FALSE, # Should determine CBH variance?
                          Lambda.time=NULL # The time points of the CBH to use
) {
  
  v <- rexp(fit$n.clusters)
  v <- v/mean(v)
  
  matrix(0, 3 + length(Lambda.time), 3 + length(Lambda.time))
}