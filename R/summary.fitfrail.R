summary.fitfrail <- function(object, Lambda.times, censored=FALSE, ...) {
  fit <- object
  
  if (!inherits(fit, "fitfrail")) 
    stop("summary.fitfrail can only be used for fitfrail objects")
  
  if (!missing(Lambda.times)) {
    if (!is.numeric(Lambda.times)) 
      stop("Lambda.times must be numeric")
    Lambda.times <- sort(Lambda.times)
  }
  
  df <- data.frame()
}