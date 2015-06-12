# This file contains R implementations of frailty distribution functions
# The corresponding C++ functions are much faster and should be used instead
# Any function ending with _r has a corresponding c function ending with _c
# The function here are provided for mainly for testing purposes

#'
#' Example of calculating phi_ with numerical integrals
#' 
phi_numerical <- function(k, H_dot, N_dot, density_params, density_fn) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * density_fn(w, density_params)
  }
  integrate(f, 0, Inf)
}

#' 
#' Example of calculating phi_ with the LT
#' 
phi_laplace <- function(k, H_dot, N_dot, density_params, density_LT) {
  density_LT(N_dot + k - 1, H_dot, density_params)*(-1)^(N_dot + k - 1)
}


################################################################################
# Gamma distribution

#'
#' Gamma random variable generation
#'
rgamma_r <- function(n, theta) {
  rgamma(n, 1/theta, 1/theta)
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
  tiv2iv <- (1/theta)^(1/theta)
  n12p <- (-1)^p
  gptiv <- gamma(p + 1/theta)
  ootps <- (1/theta + s)
  n1otmp <- (-1/theta - p)
  
  term1 <- tiv2iv * n12p * gptiv * ootps^n1otmp * 
    ( log(ootps)/theta^2 - n1otmp/(theta^2 * ootps) )
  
  term2 <- tiv2iv * n12p * (-1/theta^2 - log(1/theta)/theta^2) * gptiv * ootps^n1otmp
  
  term3 <- (1/theta)^(1/theta + 2) * n12p * gptiv * digamma(p + 1/theta) * ootps^n1otmp
  
  term4 <- (1/theta)^(1/theta + 2) * n12p * digamma(1/theta) * gptiv * ootps^n1otmp
  
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
  #   smallx <- x[x<thresh]
  #   largex <- x[x>=thresh]
  #   dposstab_approx(x, alpha, alpha, 0)
  f <- function(x) {
    k <- 1:K
    -1 / (pi * x) * sum(
      gamma(k * alpha + 1) / factorial(k) * 
        (-x ^ (-alpha)) ^ k * sin(alpha * k * pi))
  }
  
  Vectorize(f)(x)
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
  w <- rexp(n)
  theta <- runif(n)*pi
  
  num1 <- sin((1 - alpha) * theta )
  num2 <- sin(alpha * theta)^( alpha/(1 - alpha) )
  den <- sin(theta)^( 1/(1 - alpha) ) 
  num3 <- num1 * num2/den
  
  (num3/w)^( (1 - alpha)/alpha )
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
dpvf_r <- function(x, alpha, K=100) {
  Vectorize(function(x) {
    dposstab_r(x*alpha^(1/alpha), alpha, K)*alpha^(1/alpha)*exp(-x)*exp(1/alpha)
  })(x)
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
  
  deriv_lt_dpvf_r_numeric(p, s, alpha)
#   (-1)^p * deriv_lt_dpvf_r(0, s, alpha) * sum(vapply(1:p, function(j)
#     lt_deriv_dpvf_coef_r(p, j, alpha) * (1 + s)^(j*alpha - p), 0
#   ))
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
# Inverse Gaussian

# TODO:


