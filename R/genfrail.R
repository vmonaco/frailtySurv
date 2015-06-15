genfrail <- function(beta = c(log(2)), # Covariate coefficients
                      theta = 2, # Frailty distribution parameter
                      frailty = c("gamma", "lognormal", "invgauss", "posstab", "pvf"), # Frailty distribution
                      censor.mu = 130, # Gaussian distribution for censoring
                      censor.sigma = 15,
                      N = 300, K = 2, # Number of families and family size
                      lambda_0 = NULL,
                      Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau,
                      Lambda_0_inv = NULL,
                      covariates = c("normal", "uniform")) {
  
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=N*K, ncol=p)
  for (j in 1:p) {
    if (covariates == "normal") {
      Z[, j] <- rnorm(N*K)
    } else if (covariates == "uniform") {
      Z[, j] <- runif(N*K)
    }
  }
  
  w <- 1
  # Generate the frailty for each family
  # In each case, only generate random values for non-degenerate distributions
  if (frailty=="gamma") {
    if ( theta != 0 )
      w <- rep(rgamma_r(N, theta), rep(K, N))
  } else if (frailty=="lognormal") {
    if ( theta != 0 )
      w <- rep(rlognormal_r(N, theta), rep(K, N))
  } else if (frailty=="invgauss") {
    if ( theta != 0 )
      w <- rep(rinvgauss_r(N, theta), rep(K, N))
  } else if (frailty=="posstab") {
    if ( theta != 0 )
      w <- rep(rposstab_r(N, theta), rep(K, N))
  } else if (frailty=="pvf") {
    if ( theta != 0 )
      w <- rep(rpvf_r(N, theta), rep(K, N))
  } else if (frailty == "none") {
    # do nothing, w is already 1
  } else {
    stop("Wrong value for frailty: ", frailty)
  }
  
  rate <- w*exp(Z%*%beta)
  
  if (sum(c(!is.null(lambda_0), !is.null(Lambda_0), !is.null(Lambda_0_inv)) > 1))
    stop("Must specify only one of: lambda_0, Lambda_0, or Lambda_0_inv")
  
  # Generate uniform random samples for the survival times
  u <- runif(N*K)
  
  # With the inverse cumulative hazard, apply Bender's method
  if (!is.null(Lambda_0_inv)) {
    obs.time <- Lambda_0_inv(-log(u)/rate)
  # Otherwise, use the method described by Crowther
  } else {
    
    if (!is.null(lambda_0)) {
      # numeric integral nested in the root finding
      Lambda_0 <- Vectorize(function(t) {
        integrate(lambda_0, 0, t)$value
      })
    }
    
    # TODO: first try to invert the function analytically
    
    # if inversion fails, then need to find root for each t
    # Similar to the method described by Crowther, 
    # except solve: log(S(t)) - log(u) = 0
    obs.time <- sapply(1:(N*K), function(ij) {
      fn <- function(t) {
        (-Lambda_0(t)*rate[ij]) - log(u[ij])
      }
      # init at something big enough to avoid a seemingly 0 Jacobian
      # TODO: use uniroot instead?
      nleqslv(100, fn, global="dbldog", control=list(allowSingular=TRUE))$x
    })
  }
  
  # Non-informative right censoring
  censor.time <- rnorm(N*K, censor.mu, censor.sigma)
  obs.status <- sign(obs.time <= censor.time)
  obs.time <- pmin(obs.time, censor.time)
  
  colnames(Z) <- paste("Z", 1:p, sep="")
  dat <- data.frame(family=rep(1:N, rep(K,N)),
                    rep=rep(1:K, N), 
                    time=obs.time, 
                    status=obs.status,
                    Z)
  return(dat)
}
