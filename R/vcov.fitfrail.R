#' Compute the variance/covariance matrix for a fitfrail object
#' 
#' @param fit a fitfrail object
#' @param boot whether to use a weighted bootstrap. If boot == FALSE, a consistent
#'             estimator is used
#' @param B number of repetitions in the weighted bootstrap. Defaults to 100.
#' @param Lambda.time time points where the variance/covariance should be 
#'                    evaluated. If Lambda.time == NULL, then the points where 
#'                    the cumulative baseline hazard increases (where failures 
#'                    occur) are used.
#' @param cores number of cores to use when computing the covariance matrix in parallel
#' 
#' @return variance/covariance matrix for the fitfrail model parameters
#' 
#' @export
vcov.fitfrail <- function(fit, boot=TRUE, B=100,
                          Lambda.time=NULL, # The time points of the CBH to use
                          cores=0
) {
  
  if (is.null(Lambda.time)) {
    Lambda.time <- fit$Lambda.time
  }
  
  # If V has already been computed and matches the size we expect
  if (!is.null(fit[["V"]]) && 
      nrow(fit[["V"]]) == (length(fit$beta)+length(fit$theta)+length(Lambda.time)))
    return(fit[["V"]])
  
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
      setNames(new.fit$Lambdafn(Lambda.time), paste("Lambda.", format(Lambda.time, nsmall=2), sep=""))
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
  
  hats <- t(simplify2array(hats))
  
  V <- cov(hats)
  
  # Cache the value for quick access later
  eval.parent(substitute(fit[["V"]] <- V))
  
  V
}