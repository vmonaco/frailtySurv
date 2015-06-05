
plot.fitfrail <- function(fit, frailty=fit$frailty, resolution=30) {
  
  fn <- Vectorize(function(beta, theta) fit$loglikfn(beta, theta, frailty))
  beta <- seq(-5, 5, length.out=resolution)
  theta <- seq(0, 10, length.out=resolution)
  z <- outer(beta, theta, fn)
  
  image2D(z, beta, theta, clab='L', NAcol = "black", 
          xlab="beta", ylab="theta", 
          main=sprintf("Loglikelihood = %.2f with beta = %.2f and theta = %.2f", 
                       fit$loglik, fit$beta, fit$theta))
}