genfrail <- function(beta = c(log(2)), # Covariate coefficients
                     covariates = c("normal", "uniform"),
                     
                     # TODO: parameters names like this:
                     # coefficients = c(log(2)),
                     # covariate.distr = c("normal", "uniform"),
                     # covariate.param = c(0,1),
                     
                     # Frailty distribution and parameter vector
                     frailty = c("gamma", "lognormal", "invgauss", "posstab", "pvf"), 
                     theta = c(2), # Frailty distribution parameter vector
                     
                     # Censoring distribution and parameter vector
                     censoring = c("normal"),
                     censor.params = c(130,15), # Gaussian distribution for censoring
                     censor.rate = 0.30,
                     
                     # Number of clusters and cluster sizes
                     N = 300, # N must be provided
                     # If K is int, then fixed cluster sizes
                     # If K is numeric vector, then sizes provided
                     # Otherwise draw random
                     K = 2, #c(2, "poisson", "pareto", "uniform"), 
                     # cluster size distributian params:
                     # K-truncated Poisson: c(lambda, truncated value)
                     # Discrete Pareto: c(alpha, max value inclusive)
                     # Uniform: c(lower, upper) inclusive
                     K.params = c(2, 0), 
                     
                     # Only one of these needs to be specified
                     # Order of preference is: Lambda_0_inv, Lambda_0, lambda_0
                     lambda_0 = function(t, tau=4.6, C=0.01) (tau*(C*t)^tau)/t,
                     Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau,
                     Lambda_0_inv = NULL, #function(t) (t^(1/4.6))/0.01,
                     # Round time to nearest 
                     round.base = NULL
                     ) {
  
  # Determine cluster sizes
  if (is.numeric(K) && length(K) == 1) {
    cluster_sizes <- c(rep(K, N));
  } else if (is.numeric(K) && length(K) > 1) {
    if (length(K) != N)
      stop("length(K) must equal N when providing cluster sizes")
    cluster_sizes <- K
  } else if (K == "poisson") {
    cluster_sizes <- rtpois(N, K.params[1], K.params[2])
  } else if (K == "dpareto") {
    cluster_sizes <- rdpareto(N, K.params[1], K.params[2])
  } else if (K == "uniform") {
    cluster_sizes <- round(runif(N, K.params[1], K.params[2]))
  } else {
    stop("Wrong value for K. Must be int, string, or numeric vector")
  }
  
  NK <- sum(cluster_sizes)
  
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=NK, ncol=p)
  for (j in 1:p) {
    if (covariates == "normal") {
      Z[, j] <- rnorm(NK)
    } else if (covariates == "uniform") {
      Z[, j] <- runif(NK)
    } else if (covariates == "zero") {
      Z[, j] <- rep(0, NK)
    }
  }
  
  cluster_frailty <- rep(1, N)
  # Generate the frailty for each family
  # In each case, only generate random values for non-degenerate distributions
  if (frailty=="gamma") {
    if ( theta != 0 )
      cluster_frailty <- rgamma_r(N, theta)
  } else if (frailty=="lognormal") {
    if ( theta != 0 )
      cluster_frailty <- rlognormal_r(N, theta)
  } else if (frailty=="invgauss") {
    if ( theta != 0 )
      cluster_frailty <- rinvgauss_r(N, theta)
  } else if (frailty=="posstab") {
    if ( theta != 1 )
      cluster_frailty <- rposstab_r(N, theta)
  } else if (frailty=="pvf") {
    if ( theta != 1 )
      cluster_frailty <- rpvf_r(N, theta)
  } else if (frailty == "none") {
    # do nothing, w is already 1
  } else {
    stop("Wrong value for frailty: ", frailty)
  }
  
  w = rep(cluster_frailty, cluster_sizes)
  rate <- w*exp(Z%*%beta)
  
  if (sum(c(!is.null(lambda_0), !is.null(Lambda_0), !is.null(Lambda_0_inv)) > 1))
    stop("Must specify only one of: lambda_0, Lambda_0, or Lambda_0_inv")
  
  # Generate uniform random samples for the survival times
  u <- runif(NK)
  
  # With the inverse cumulative hazard, apply Bender's method
  if (!is.null(Lambda_0_inv)) {
    fail.time <- Lambda_0_inv(-log(u)/rate)
  # Otherwise, use the method described by Crowther
  } else {
    if (is.null(Lambda_0) && !is.null(lambda_0)) {
      # numeric integral nested in the root finding
      Lambda_0 <- Vectorize(function(t) {
        if (t <= 0) {
          return(0)
        }
        integrate(lambda_0, 0, t, subdivisions = 1000L)$value
      })
    } else if (is.null(Lambda_0)) {
      stop("Need baseline hazard if cumulative baseline hazard is not specified.")
    }
    
    # TODO: first try to invert the function analytically
    
    # if inversion fails, then need to find root for each t
    # Similar to the method described by Crowther, 
    # except solve: log(S(t)) - log(u) = 0
    fail.time <- sapply(1:(NK), function(ij) {
      fn <- function(t) {
        (-Lambda_0(t)*rate[ij]) - log(u[ij])
      }
      # nleqslv(100, fn, global="dbldog", control=list(allowSingular=TRUE))$x
      uniroot(fn, lower=0, upper=100, extendInt="yes")$root
    })
  }
  
#   Surv.fn <- ecdf(fail.time)
#   Fail.fn <- function(t) 1 - Surv.fn(t)
#   mu <- uniroot(function(mu) 
#     censor.rate - integrate(function(t) 
#       Fail.fn(t)*dnorm(t, mu, censor.params[2]), 0, Inf, subdivisions=1e3)$value,
#     lower=0, upper=100, extendInt="up")$root
  
  # Non-informative right censoring
  censor.time <- rnorm(NK, censor.params[1], censor.params[2])
  obs.status <- sign(fail.time <= censor.time)
  obs.time <- pmin(fail.time, censor.time)
  
  if (is.numeric(round.base)) {
    obs.time <- round.base*floor(obs.time/round.base + 0.5)
  }
  
  colnames(Z) <- paste("Z", 1:p, sep="")
  dat <- data.frame(family=rep(1:N, cluster_sizes),
                    rep=unlist(lapply(cluster_sizes, function(m) rep(1:m))), 
                    time=obs.time, 
                    status=obs.status,
                    Z)
  return(dat)
}
