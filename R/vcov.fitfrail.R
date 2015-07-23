#' 
#' Compute the variance/covariance matrix for a fitted object
#' 
#' @export
vcov.fitfrail <- function(fit, boot=TRUE, B=100,
                          Lambda.time=NULL, # The time points of the CBH to use
                          cores=0
) {
  
  if (is.null(Lambda.time)) {
    Lambda.time <- fit$Lambda.time
  }
  
  fn <- function(s) {
    set.seed(s)
    
    weights <- rexp(fit$n.clusters)
    weights <- weights/mean(weights)
    
    new.call <- fit$call
    new.call$weights <- weights
    new.fit <- eval(new.call)
    
    c(
      new.fit$beta,
      new.fit$theta,
      setNames(new.fit$Lambdafn(Lambda.time), paste("Lambda.", 1:length(Lambda.time), sep=""))
    )
  }
  
  # To make the result reproducable when parallelized, seed each run 
  seeds <- sample(1:1e7, B, replace=TRUE)
  
  use.parallel <- cores != 1
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
    # Run in parallel and make a row-observation matrix
    if (cores <= 0) {
      cores <- parallel::detectCores() + cores
    }
    hats <- parallel::mclapply(seeds, fn, mc.preschedule=FALSE, 
                              mc.set.seed=TRUE, mc.cores=cores, mc.silent=FALSE)
  } else if (use.parallel && !have.parallel) {
    # Run in serial with a warning
    warning('Cannot run vcov.fitfrail in parallel: requires the "parallel" package.')
    hats <- lapply(seeds, fn)
  } else {
    # Run in serial
    hats <- lapply(seeds, fn)
  }
  return(hats)
  hats <- t(simplify2array(hats))
  
  cov(hats)
}