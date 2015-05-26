dgamma_ <- function(x, theta) {
  dgamma(x, 1/theta, 1/theta)
}

dgamma_dtheta <- function(x, theta, deriv_idx) {
  # deriv_idx ignored here since there is only one parameter
  if (is.infinite(x)) return(0)
  
  term1 <- (x/theta)^(1/theta - 1)
  term2 <- exp(-x/theta)
  term3 <- log(theta/x) + digamma(1/theta) + x - 1
  
  (term1 * term2 * term3)/(gamma(1/theta)*theta^3)
}
# 
# dlnorm_dtheta <- function(x, sigma) {
#   term1 <- log(x)^2 * exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^4)
#   term2 <- exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^2)
#     
#   term1 - term2
# }
