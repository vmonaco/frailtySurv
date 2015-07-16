genfrail <- function(
                     beta = c(log(2)), # Covariate coefficients
                     covar.distr = c("normal", "uniform", "zero"),
                     covar.param = c(0,1),
                     
                     # TODO: parameters names like this:
                     # coefficients = c(log(2)),
                     # covariate.distr = c("normal", "uniform"),
                     
                     # Frailty distribution and parameter vector
                     frailty = c("gamma", "lognormal", "invgauss", 
                                 "posstab", "pvf", "none"), 
                     theta = c(2), # Frailty distribution parameter vector
                     
                     # Censoring distribution and parameters vector
                     censor.distr = c("normal", "lognormal", "none"),
                     censor.rate = NULL, # If specified, overrides mu
                     censor.mu = 130, # 5
                     censor.sigma = 15, # 0.5
                     
                     # Number of clusters and cluster sizes
                     N = 300, # N must be provided
                     # If K is int, then fixed cluster sizes
                     # If K is numeric vector, then sizes provided
                     # Otherwise draw random
                     K = 2, #c(2, "poisson", "tzeta", "duniform"), 
                     # cluster size distributian params:
                     # K-truncated Poisson: c(lambda, truncated value)
                     # Discrete Pareto: c(alpha, max value inclusive)
                     # Uniform: c(lower, upper) inclusive
                     K.params = c(2, 0), 
                     
                     # Only one of these needs to be specified
                     # Order of preference is: Lambda_0_inv, Lambda_0, lambda_0
                     lambda_0 = function(t, tau=4.6, C=0.01) (tau*(C*t)^tau)/t,
                     Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau,
                     Lambda_0_inv = NULL, #function(t, tau=4.6, C=0.01) (t^(1/tau))/C,
                     # Round time to nearest 
                     round.base = NULL
                     ) {
  
  # Determine cluster sizes
  if (is.numeric(K) && length(K) == 1) {
    cluster.sizes <- c(rep(K, N));
  } else if (is.numeric(K) && length(K) > 1) {
    if (length(K) != N)
      stop("length(K) must equal N when providing cluster sizes")
    cluster.sizes <- K
  } else if (K == "poisson") {
    cluster.sizes <- rtpois(N, K.params[1], K.params[2])
  } else if (K == "tzeta") {
    cluster.sizes <- rtzeta(N, K.params[1], K.params[2], K.params[3])
  } else if (K == "duniform") {
    cluster.sizes <- round(runif(N, K.params[1], K.params[2]))
  } else {
    stop("Wrong value for K. Must be int, string, or numeric vector")
  }
  
  NK <- sum(cluster.sizes)
  
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=NK, ncol=p)
  for (j in 1:p) {
    if (covar.distr == "normal") {
      Z[, j] <- rnorm(NK)
    } else if (covar.distr == "uniform") {
      Z[, j] <- runif(NK)
    } else if (covar.distr == "zero") {
      Z[, j] <- rep(0, NK)
    }
  }
  
  cluster.frailty <- rep(1, N)
  # Generate the frailty for each family
  # In each case, only generate random values for non-degenerate distributions
  if (frailty=="gamma") {
    if ( theta != 0 )
      cluster.frailty <- rgamma_r(N, theta)
  } else if (frailty=="lognormal") {
    if ( theta != 0 )
      cluster.frailty <- rlognormal_r(N, theta)
  } else if (frailty=="invgauss") {
    if ( theta != 0 )
      cluster.frailty <- rinvgauss_r(N, theta)
  } else if (frailty=="posstab") {
    if ( theta != 1 )
      cluster.frailty <- rposstab_r(N, theta)
  } else if (frailty=="pvf") {
    if ( theta != 1 )
      cluster.frailty <- rpvf_r(N, theta)
  } else if (frailty == "none") {
    # do nothing, w is already 1
  } else {
    stop("Wrong value for frailty: ", frailty)
  }
  
  w = rep(cluster.frailty, cluster.sizes)
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
    
    # TODO: first try to invert the function analytically, maybe using Ryacas?
    
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
  
  if (censor.distr == "normal") {
    censor.density <- dnorm
    censor.random <- rnorm
    censor.quantile <- qnorm
  } else if (censor.distr == "lognormal") {
    censor.density <- dlnorm
    censor.random <- rlnorm
    censor.quantile <- qlnorm
  } else if (censor.distr == "none") {
    warning("Censoring distribution is set to none")
  } else {
    stop("Censoring distribution must be normal or lognormal")
  }
  
  if (censor.distr == "none") {
    obs.status <- rep(1, NK)
    obs.time <- fail.time
  } else {
    if (!is.null(censor.rate)) {
      if (!is.null(censor.mu))
        warning("Censoring rate will override mean")
      
      stopifnot(0 <= censor.rate && censor.rate <= 1)
      
      # Empirical survivor and hazard functions
      esurv <- ecdf(fail.time)
      efail <- function(t) 1 - esurv(t)
      
      # Split the integral up into 2 parts since it may actually be ~0 for 
      # most of the interval. Make this split at the 50th percentile
      censor.mu <- uniroot(function(mu) 
        censor.rate - 
          (integrate(function(t) efail(t)*censor.density(t, mu, censor.sigma), 
            -Inf, censor.quantile(0.5, mu, censor.sigma), subdivisions=1e3)$value +
           integrate(function(t) efail(t)*censor.density(t, mu, censor.sigma), 
            censor.quantile(0.5, mu, censor.sigma), Inf, subdivisions=1e3)$value),
        lower=0, upper=100, extendInt="up")$root
    }
    
    # Non-informative right censoring
    censor.time <- censor.random(NK, censor.mu, censor.sigma)
    obs.status <- sign(fail.time <= censor.time)
    obs.time <- pmin(fail.time, censor.time)
  }
  
  if (is.numeric(round.base)) {
    obs.time <- round.base*floor(obs.time/round.base + 0.5)
  }
  
  colnames(Z) <- paste("Z", 1:p, sep="")
  dat <- data.frame(family=rep(1:N, cluster.sizes),
                    rep=unlist(lapply(cluster.sizes, function(m) rep(1:m))), 
                    time=obs.time, 
                    status=obs.status,
                    Z)
  return(dat)
}
