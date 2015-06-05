genfrail <- function(beta = c(log(2)), # Covariate coefficients
                     frailty = c("gamma", "lognormal", "posstab"), # Frailty distribution
                     censor.mu = 130, # Gaussian distribution for censoring
                     censor.sigma = 15,
                     theta = 2, # Frailty distribution parameter
                     N = 300, K = 2, # Number of families and family size
                     tau = 4.6,
                     C = 0.01,
                     covariates = c("normal", "uniform") ) {
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=N*K, ncol=p)
  for (j in 1:p) {
    if (covariates=="normal") {
      Z[, j] <- rnorm(N*K)
    } else if (covariates=="uniform") {
      Z[, j] <- runif(N*K)
    }
  }
  betaZ <- Z%*%beta
  
  # Generate the frailty for each family
  if (frailty=="gamma") {
    w <- 1
    if ( theta != 0 ) {
      w <- rep(rgamma_r(N, theta), rep(K, N))
    }
    rate <- exp(betaZ)*w 
  } else if (frailty=="lognormal") {
    w <- 1
    if ( theta != 0 ) {
      w <- rep(rnorm_r(N, theta), rep(K, N))
    }
    rate <- exp(betaZ + w)
  } else if (frailty=="posstab") {
    w <- 1
    if ( theta != 0 ) {
      w <- rep(rposstab_r(N, theta), rep(K, N))
    }
    rate <- exp(betaZ)*w 
  }
  
  # Survival time is similar to that of a late onset disease with default params
  u <- runif(N*K)
  obs.time <- ((-log(u)/rate)^(1/tau))/C
  
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