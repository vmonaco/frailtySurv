# R implementations of frailty densities and generation. 
# The corresponding C++ functions are much faster and should be used instead

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
deriv_dgamma_r <- function(x, theta, deriv_idx) {
  # deriv_idx ignored here since there is only one parameter
  term1 <- (x/theta)^(1/theta - 1)
  term2 <- exp(-x/theta)
  term3 <- log(theta/x) + digamma(1/theta) + x - 1
  (term1 * term2 * term3)/(gamma(1/theta)*theta^3)
}

#'
#' Laplace transform of the gamma density, as defined above
#' 
lt_dgamma_r <- function(p, s, theta) {
  (-1)^p * (1/theta)^(1/theta) * ((1/theta) + s)^(-((1/theta) + p)) * 
    gamma((1/theta) + p)/gamma((1/theta))
}

#' 
#' pth moment of the gamma distribution
#' 
lt_dgamma_moment_r <- function(p, theta) {
  (-1)^p * lt_dgamma_r(p, 0, theta)
}

#' 
#' Example of calculating phi_1 with the LT
#' 
lt_phi_dgamma <- function(k, H_dot, N_dot, theta) {
  lt_dgamma_r(N_dot + k - 1, H_dot, theta)*(-1)^(N_dot + k - 1)
}

#'
#' Example of calculating phi_1 with numerical integrals
#' 
phi_dgamma <- function(k, H_dot, N_dot, theta) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * dgamma_r(w, theta)
  }
  integrate(f, 0, Inf)
}

################################################################################
# Positive stable distribution

#' 
#' Positive stable density
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

lt_dposstab_r <- function(p, s, alpha) {
  if (p < 1) {
    return(exp(-alpha*(s^alpha)/alpha))
  }
  
  (-1)^p * lt_dposstab_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef(p, j, alpha) * alpha^j * s^(j*alpha - p), 0
  ))
}

lt_dposstab_moment_r <- function(alpha, p) { 
  (-1)^p * lt_dposstab_r(p, 0, alpha)
}

lt_phi_dposstab <- function(k, H_dot, N_dot, alpha) {
  lt_dposstab_r(N_dot + k - 1, H_dot, alpha)*(-1)^(N_dot + k - 1)
}

phi_dposstab <- function(k, H_dot, N_dot, alpha) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * dposstab_r(w, alpha)
  }
  integrate(f, 0.01, Inf)
}

################################################################################
# Power variance function distribution (PVF)

dpvf_r <- function(x, alpha, K=100) {
  Vectorize(function(x) {
    dposstab_r(x*alpha^(1/alpha), alpha, K)*alpha^(1/alpha)*exp(-x)*exp(1/alpha)
  })(x)
}

lt_dpvf_coef <- function(p, j, alpha) {
  if (p == 2 && j == 1) {
    return(1 - alpha)
  }
  
  if (p == 2 && j == 2) {
    return(1)
  }
  
  if (p == 3 && j == 1) {
    return((1 - alpha)*(2 - alpha))
  }
  
  if (p == 3 && j == 2) {
    return(3*(1 - alpha))
  }
  
  if (p == 3 && j == 3) {
    return(1)
  }
  
  if (j == 1) {
    return(gamma(p - alpha)/gamma(1 - alpha))
  }
  
  if (p == j) {
    return(1)
  }
  
  return(lt_dpvf_coef(p-1, j-1, alpha) + lt_dpvf_coef(p-1, j, alpha)*((p-1) - j*alpha))
}

lt_dpvf_r <- function(p, s, alpha) {
  if (p < 1) {
    return(exp(-((1 + s)^alpha - 1)/alpha))
  }
  
  (-1)^p * lt_dpvf_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef(p, j, alpha) * (1 + s)^(j*alpha - p), 0
    ))
}

lt_dpvf_moment_r <- function(alpha, p) { 
  (-1)^p * lt_dpvf_r(p, 0, alpha)
}

lt_phi_dpvf <- function(k, H_dot, N_dot, alpha) {
  lt_dpvf_r(N_dot + k - 1, H_dot, alpha)*(-1)^(N_dot + k - 1)
}

phi_dpvf <- function(k, H_dot, N_dot, alpha) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * dpvf_r(w, alpha)
  }
  integrate(f, 0.01, Inf)
}

################################################################################
# Log-normal

#' 
#' Normal random variable generation
#' 
rnorm_r <- function(n, theta) {
  rnorm(n, 0, sqrt(theta))
}

#' 
#' Deriv of log norm density wrt. theta
#' 
dlnorm_dtheta_r <- function(x, theta) {
  term1 <- log(x)^2 * exp(-(log(x)^2)/(2*theta^2)) / (x*sqrt(2*pi)*theta^4)
  term2 <- exp(-(log(x)^2)/(2*theta^2)) / (x*sqrt(2*pi)*theta^2)
  
  term1 - term2
}

################################################################################
# Inverse Gaussian



