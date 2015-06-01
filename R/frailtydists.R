
#' 
#' Gamma density
#' 
# dgamma_ <- function(x, theta) {
#   dgamma(x, 1/theta, 1/theta)
# }

#' 
#' Deriv of gamma density wrt. theta
#' 
# deriv_dgamma <- function(x, theta, deriv_idx) {
#   # deriv_idx ignored here since there is only one parameter
#   term1 <- (x/theta)^(1/theta - 1)
#   term2 <- exp(-x/theta)
#   term3 <- log(theta/x) + digamma(1/theta) + x - 1
#   (term1 * term2 * term3)/(gamma(1/theta)*theta^3)
# }

#' 
#' Positive stable density
#' 
dposstab <- function(x, alpha, K) {
  # TODO: for small x, use PVF saddle point approximation
  k <- 1:K
  -1 / (pi * x) * sum(
    gamma(k * alpha + 1) / factorial(k) * 
      (-x ^ (-alpha)) ^ k * sin(alpha * k * pi))
}

#' 
#' Positive stable density
#' 
dpvf_sp_approx <- function(x, alpha, delta, theta) {
  v <- -(2 - alpha)/(2 * (1 - alpha))
  (2 * pi * (1 - alpha))^(1/2) * delta^(1/(2*(1 - alpha))) * x^v * 
    exp(-theta * x - delta^(1/(1 - alpha)) * (1/alpha - 1) * 
          x^(-alpha/(1 - alpha)) + (delta * theta^alpha)/alpha)
}

# 
# dlnorm_dtheta <- function(x, sigma) {
#   term1 <- log(x)^2 * exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^4)
#   term2 <- exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^2)
#     
#   term1 - term2
# }
