# This file contains R implementations of frailty distribution functions
# The corresponding C++ functions are much faster and should be used instead
# Any function ending with _r has a corresponding c function ending with _c
# The function here are provided for mainly for testing purposes

#'
#' Example of calculating phi_ with numerical integrals
#' 
phi_numerical <- function(k, N_dot, H_dot, density_params, density_fn) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * density_fn(w, density_params)
  }
  integrate(f, 0, Inf)
}

#' 
#' Example of calculating phi_ with the LT
#' 
phi_laplace <- function(k, N_dot, H_dot, density_params, density_LT) {
  density_LT(N_dot + k - 1, H_dot, density_params)*(-1)^(N_dot + k - 1)
}

num_vs_lt <- function(k, N_dot, H_dot, density_params, density_fn, density_LT) {
  list(phi_numerical(k, N_dot, H_dot, density_params, density_fn),
      phi_laplace(k, N_dot, H_dot, density_params, density_LT))
}


################################################################################
# Gamma distribution

#'
#' Gamma random variable generation
#'
rgamma_r <- function(n, theta) {
  rgamma(n, 1/theta, 1/theta)
  
  # Or using rlaptrans:
  # rlaptrans(n, lt_dgamma_r, p=0, theta=theta)
}

#' 
#' Gamma density
#' 
dgamma_r <- function(x, theta) {
  dgamma(x, 1/theta, 1/theta)
}

#' 
#' Deriv of gamma density wrt. theta
#' 
deriv_dgamma_r <- function(x, theta) {
  # deriv_idx ignored here since there is only one parameter
  term1 <- (x/theta)^(1/theta - 1)
  term2 <- exp(-x/theta)
  term3 <- log(theta/x) + digamma(1/theta) + x - 1
  (term1 * term2 * term3)/(gamma(1/theta)*theta^3)
}

#' 
#' Deriv of gamma density wrt. theta, evaluated numerically
#' 
deriv_dgamma_r_numeric <- function(x, theta) {
  grad(function(theta) dgamma_r(x, theta), theta)
}

#'
#' Laplace transform of the gamma density, as defined above
#' 
lt_dgamma_r <- function(p, s, theta) {
  (-1)^p * 
    (1/theta)^(1/theta) * 
    ((1/theta) + s)^(-((1/theta) + p)) * 
    gamma((1/theta) + p)/gamma((1/theta))
}

#'
#' Derivative of the Laplace transform of the gamma density, wrt theta,
#' pth derivative wrt. s
#' 
deriv_lt_dgamma_r <- function(p, s, theta) {
  tmp1 <- (1/theta)^(1/theta)
  tmp2 <- (-1)^p
  tmp3 <- gamma(p + 1/theta)
  tmp4 <- (1/theta + s)
  tmp5 <- (-1/theta - p)
  
  term1 <- tmp1 * tmp2 * tmp3 * tmp4^tmp5 * 
    ( log(tmp4)/theta^2 - tmp5/(theta^2 * tmp4) )
  
  term2 <- tmp1 * tmp2 * (-1/theta^2 - log(1/theta)/theta^2) * tmp3 * tmp4^tmp5
  
  term3 <- (1/theta)^(1/theta + 2) * tmp2 * tmp3 * digamma(p + 1/theta) * tmp4^tmp5
  
  term4 <- (1/theta)^(1/theta + 2) * tmp2 * digamma(1/theta) * tmp3 * tmp4^tmp5
  
  (term1 + term2 - term3 + term4)/gamma(1/theta)
}

#' 
#' Same as above, evaluated numerically
#' 
deriv_lt_dgamma_r_numeric <- function(p, s, theta) {
  grad(function(theta) lt_dgamma_r(p, s, theta), theta)
}

#' 
#' pth moment of the gamma distribution, evaluated by the gamma density LT
#' 
lt_dgamma_moment_r <- function(p, theta) {
  (-1)^p * lt_dgamma_r(p, 0, theta)
}

################################################################################
# Positive stable distribution

#' 
#' Positive stable density
#' TODO: does not converge when x~0 for small alpha. The point where convergance
#' fails depends on alpha
#' eg. with alpha = 0.5, converges at x = 0.05, but not at 0.01
#' 
dposstab_r <- function(x, alpha, K=100) {
  
  f <- function(x) {
#     if (x < exp(-5*(1-alpha))) {
#       1e-3 # TODO: temporary
#     } else {
      k <- 1:K
      -1 / (pi * x) * sum(
        gamma(k * alpha + 1) / factorial(k) * 
          (-x ^ (-alpha)) ^ k * sin(alpha * k * pi))
    # }
  }
  
  Vectorize(f)(x)
}

# Where is dposstab divergent?
dposstab_divergent <- function(alpha, K=100, resolution=1e-3, eps=1e-5) {
  Vectorize(function(alpha) {
    k <- 1:K
    f <- Vectorize(function(x) {
        gamma(k * alpha + 1) / factorial(k) * 
          (-x ^ (-alpha)) ^ k * sin(alpha * k * pi)
    })
    
    x <- seq(0, 1, by=resolution)
    mindiff <- apply(f(x), 2, function(x) min(abs(diff((x)))))
    
    if (any(mindiff[!is.na(mindiff)] > eps)) {
      return(x[max(which(mindiff > eps)) + 1])
    } else {
      return(x[2])
    }
  })(alpha)
}

#' Random generation from a positive stable distribution
#'
#' \code{rposstab} generates data.
#'
#' Based on Chambers
#' 
#' @param n number of samples to generate
#' @param alpha index of the distribution, where 0 < alpha < 1
#' @return w
#' 
#'@export
rposstab_r <- function(n, alpha) {
  rlaptrans(n, lt_dposstab_r, p=0, alpha=alpha)
#   w <- rexp(n)
#   theta <- runif(n)*pi
#   
#   num1 <- sin((1 - alpha) * theta )
#   num2 <- sin(alpha * theta)^( alpha/(1 - alpha) )
#   den <- sin(theta)^( 1/(1 - alpha) ) 
#   num3 <- num1 * num2/den
#   
#   (num3/w)^( (1 - alpha)/alpha )
}

#'
#' Since there's no closed-form expression for the posstab density, the 
#' derivative is evaluated numerically. This isn't needed for estimation
#' 
deriv_dposstab_r_numeric <- function(x, alpha) {
  Vectorize(function(x) grad(function(alpha) dposstab_r(x, alpha), alpha))(x)
}

#'
#' Laplace transform of the positive stable distribution
#' 
lt_dposstab_r <- function(p, s, alpha) {
  if (p == 0) {
    return(exp(-alpha*(s^alpha)/alpha))
  }
  
  (-1)^p * lt_dposstab_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef_r(p, j, alpha) * alpha^j * s^(j*alpha - p), 0
  ))
}

#' 
#' Derivative wrt. alpha of the posstab LT
deriv_lt_dposstab_r_numeric <- function(p, s, alpha) {
  grad(function(alpha) lt_dposstab_r(p, s, alpha), alpha)
}

#'
#' Derivative wrt. alpha of the posstab LT.
#' TODO: determine coefficients for the pth derivatev wrt. s
#' currently only p = 0 (0th derivative wrt. s) has a closed-form expression
#' 
deriv_lt_dposstab_r <- function(p, s, alpha) {
  if (p == 0) {
    return(-exp(-s^alpha) * s^alpha * log(s))
  }
  
  deriv_lt_dposstab_r_numeric(p, s, alpha)
#   (-1)^p * deriv_lt_dposstab_r(0, s, alpha) * sum(vapply(1:p, function(j)
#     lt_deriv_dpvf_coef_r(p, j, alpha) * alpha^j * s^(j*alpha - p), 0
#   ))
}

#' 
#' pth moment of the posstab distribution, evaluated using the LT
#' 
lt_dposstab_moment_r <- function(alpha, p) { 
  (-1)^p * lt_dposstab_r(p, 0, alpha)
}

################################################################################
# Power variance function distribution (PVF)

#'
#' See comments from dposstab_r, same convergence issue applies here
#' 
dpvf_r <- function(x, alpha, K=100) {
  Vectorize(function(x) {
    dposstab_r(x*alpha^(1/alpha), alpha, K)*alpha^(1/alpha)*exp(-x)*exp(1/alpha)
  })(x)
}

#'
#' 
#' 
deriv_dpvf_r_numeric <- function(x, alpha) {
  Vectorize(function(x) grad(function(alpha) dpvf_r(x, alpha), alpha))(x)
}

#' 
#' 
#' 
rpvf_r <- function(n, alpha) {
  rlaptrans(n, lt_dpvf_r, p=0, alpha=alpha)
}

#'
#' Coefficients for the PVF LT derivatives
#' 
lt_dpvf_coef_r <- function(p, j, alpha) {
  if (p == 2 && j == 1) return(1 - alpha)
  if (p == 2 && j == 2) return(1)
  if (p == 3 && j == 1) return((1 - alpha)*(2 - alpha))
  if (p == 3 && j == 2) return(3*(1 - alpha))
  if (p == 3 && j == 3) return(1)
  if (j == 1) return(gamma(p - alpha)/gamma(1 - alpha))
  if (p == j) return(1)
  
  return(lt_dpvf_coef_r(p - 1, j - 1, alpha) + lt_dpvf_coef_r(p - 1, j, alpha)*((p - 1) - j*alpha))
}

#' 
#' Coefficients for the PVF LT derivates, wrt. alpha using Hougaard's definition
#' Still need coefficients with alpha = delta
#' 
lt_deriv_dpvf_coef_r <- function(p, j, alpha) {
  if (j == 1) {
    factor1 <- gamma(p - alpha)/gamma(1 - alpha)
    factor2 <- sum(1/((1:(p-1)) - alpha))
    return(factor1*factor2)
  }
  
  if (p == j) return(0)
  
  return(lt_deriv_dpvf_coef_r(p - 1, j - 1, alpha) + 
           lt_deriv_dpvf_coef_r(p - 1, j, alpha) * ((p - 1) - j*alpha) -
           j * lt_dpvf_coef_r(p - 1, j, alpha)) # last coef is not a deriv coef
}

#' 
#' Laplace transform of the one-parameter PVF distribution defined above
lt_dpvf_r <- function(p, s, alpha) {
  if (p == 0) {
    return(exp(-((1 + s)^alpha - 1)/alpha))
  }
  
  (-1)^p * lt_dpvf_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef_r(p, j, alpha) * (1 + s)^(j*alpha - p), 0
    ))
}

#' 
#' pth moment of the one-parameter PVF distribution, as defined above
lt_dpvf_moment_r <- function(alpha, p) { 
  (-1)^p * lt_dpvf_r(p, 0, alpha)
}

#'
#' Derivative wrt. alpha of the PVF LT, evaluated numerically
#' 
deriv_lt_dpvf_r_numeric <- function(p, s, alpha) {
  grad(function(alpha) lt_dpvf_r(p, s, alpha), alpha)
}

#' Similar issue as deriv_lt_dposstab_r, only closed-form expression for the
#' 0th derivative wrt. s right now
deriv_lt_dpvf_r <- function(p, s, alpha) {
  if (p == 0) {
    factor1 <- exp((1 - (s + 1)^alpha)/alpha)
    factor2_term1 <- -(1 - (s + 1)^alpha)/alpha^2
    factor2_term2 <- -((s + 1)^alpha * log(s + 1))/alpha
    return(factor1*(factor2_term1 + factor2_term2))
  }
  
  # deriv_lt_dpvf_r_numeric(p, s, alpha)
   (-1)^p * deriv_lt_dpvf_r(0, s, alpha) * sum(vapply(1:p, function(j)
     lt_deriv_dpvf_coef_r(p, j, alpha) * (1 + s)^(j*alpha - p), 0
   ))
}


################################################################################
# Log-normal

#' 
#' Normal random variable generation
#' 
rnorm_r <- function(n, theta) {
  rnorm(n, 0, sqrt(theta))
}

rlognormal_r <- function(n, theta) {
  rlnorm(n, 0, sqrt(theta))
}

dlognormal_r <- function(x, theta) {
  dlnorm(x, 0, sqrt(theta))
}

deriv_dlognormal_r <- function(x, theta) {
  term1_numer <- log(x)^2 * exp(-(log(x)^2)/(2*theta))
  term1_denom <- 2 * sqrt(2*pi) * theta^(5/2) * x
  term2_numer <- exp(-(log(x)^2)/(2*theta))
  term2_denom <- 2 * sqrt(2*pi) * theta^(3/2) * x
  term1_numer/term1_denom - term2_numer/term2_denom
}

deriv_dlognormal_r_numeric <- function(x, theta) {
  grad(function(theta) dlognormal_r(x, theta), theta)
}

################################################################################
# TODO: Inverse Gaussian

################################################################################
# K-truncated Poisson

# Sample from a k-truncated Poisson, where support is strictly greater than k
# modifed from: https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html
rtpois <- function(n, lambda, k=0) {
  qpois(runif(n, sum(dpois(0:k, lambda)), 1), lambda)
}

# Expected value of the k truncated Poisson
etpois <- function(lambda, k=0) {
  one2k <- 1:k
  zero2k <- 0:k
  
  if (k > 0)
    numer_factor2 <- sum(lambda^one2k / factorial(one2k - 1))
  else
    numer_factor2 <- 0
  
  numer <- lambda - exp(-lambda) * numer_factor2
  denom <- 1 - exp(-lambda) * sum(lambda^zero2k / factorial(zero2k))
  numer/denom
}

################################################################################
# Discrete truncated Pareto

# Density
ddpareto <- function(x, alpha){
  exp(-alpha*log(x) - log(zeta(alpha)))
}

# The "easy" way of simulating a discrete Pareto, but not exact
rdpareto <- function(n, alpha, cuttoff=10000){
  sample(x=1:x.max, size=n, replace=TRUE,
         prob=ddpareto(alpha, x=1:cuttoff))
}

################################################################################

#' 
#' rlaptrans taken from Martin Ridout: 
#' Ridout, M.S. (2009) Generating random numbers from a distribution specified 
#' by its Laplace transform. Statistics and Computing, 19, 439-450.
#' http://www.kent.ac.uk/smsas/personal/msr/rlaptrans.html
#' 
#======================================================================================
rlaptrans <- function(n, ltpdf, ..., tol=1e-7, x0=1, xinc=2, m=11, L=1, A=19, nburn=38)
#======================================================================================
{
  
  #---------------------------------------------------------
  # Function for generating a random sample of size n from a 
  # distribution, given the Laplace transform of its p.d.f.
  #---------------------------------------------------------
  
  maxiter = 500
  
  # -----------------------------------------------------
  # Derived quantities that need only be calculated once,
  # including the binomial coefficients
  # -----------------------------------------------------
  nterms = nburn + m*L
  seqbtL = seq(nburn,nterms,L)
  y = pi * (1i) * seq(1:nterms) / L
  expy = exp(y)
  A2L = 0.5 * A / L
  expxt = exp(A2L) / L
  coef = choose(m,c(0:m)) / 2^m
  
  
  # --------------------------------------------------
  # Generate sorted uniform random numbers. xrand will
  # store the corresponding x values
  # --------------------------------------------------
  u = sort(runif(n), method="qu")
  xrand = u
  
  #------------------------------------------------------------
  # Begin by finding an x-value that can act as an upper bound
  # throughout. This will be stored in upplim. Its value is
  # based on the maximum value in u. We also use the first
  # value calculated (along with its pdf and cdf) as a starting
  # value for finding the solution to F(x) = u_min. (This is
  # used only once, so doesn't need to be a good starting value
  #------------------------------------------------------------
  t = x0/xinc
  cdf = 0   
  kount0 = 0
  set1st = FALSE
  while (kount0 < maxiter & cdf < u[n]) {
    t = xinc * t
    kount0 = kount0 + 1
    x = A2L / t
    z = x + y/t
    ltx = ltpdf(x, ...)
    ltzexpy = ltpdf(z, ...) * expy
    par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
    par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
    pdf = expxt * sum(coef * par.sum[seqbtL]) / t
    cdf = expxt * sum(coef * par.sum2[seqbtL]) / t
    if (!set1st & cdf > u[1]) {
      cdf1 = cdf
      pdf1 = pdf
      t1 = t
      set1st = TRUE
    }
  }
  if (kount0 >= maxiter) {
    stop('Cannot locate upper quantile')
  }
  upplim = t
  
  #--------------------------------
  # Now use modified Newton-Raphson
  #--------------------------------
  
  lower = 0
  t = t1
  cdf = cdf1
  pdf = pdf1
  kount = numeric(n)
  
  maxiter = 1000
  
  for (j in 1:n) {
    
    #-------------------------------
    # Initial bracketing of solution
    #-------------------------------
    upper = upplim
    
    kount[j] = 0
    while (kount[j] < maxiter & abs(u[j]-cdf) > tol) {
      kount[j] = kount[j] + 1
      
      #-----------------------------------------------
      # Update t. Try Newton-Raphson approach. If this 
      # goes outside the bounds, use midpoint instead
      #-----------------------------------------------
      t = t - (cdf-u[j])/pdf 
      if (t < lower | t > upper) {
        t = 0.5 * (lower + upper)
      }
      
      #----------------------------------------------------
      # Calculate the cdf and pdf at the updated value of t
      #----------------------------------------------------
      x = A2L / t
      z = x + y/t
      ltx = ltpdf(x, ...)
      ltzexpy = ltpdf(z, ...) * expy
      par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
      par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
      pdf = expxt * sum(coef * par.sum[seqbtL]) / t
      cdf = expxt * sum(coef * par.sum2[seqbtL]) / t
      
      #------------------
      # Update the bounds 
      #------------------
      if (cdf <= u[j]) {
        lower = t}
      else {
        upper = t}
    }
    if (kount[j] >= maxiter) {
      warning('Desired accuracy not achieved for F(x)=u')
    }
    xrand[j] = t
    lower = t
  }
  
  if (n > 1) {
    rsample <- sample(xrand) }
  else {
    rsample <- xrand} 
  rsample
}
