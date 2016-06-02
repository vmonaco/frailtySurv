# Misc. functions

################################################################################
# Measures of dependence

tau.empirical <- function(theta, frailty, N=1000) {
  dat <- genfrail(N=N, K=2, theta=theta, frailty=frailty, 
                  censor.distr="none", covar.distr="none",
                  Lambda_0_inv=function(t, c=0.01, d=4.6) (t^(1/d))/c)
  xy <- matrix(dat$time, ncol=2, byrow=TRUE)
  cor(xy[,1], xy[,2], method="kendall")
}

tau.numerical <- function(theta, frailty) {
  lt <- ltfrailty[[frailty]]
  
  4 * integrate(Vectorize(function(s) {
    s * lt(0, s, theta) * lt(2, s, theta)
  }), 0, Inf)$value - 1
}

# Numerically determine theta given Kendall's tau
theta.given.tau <- function(tau, frailty) {
  uniroot(function(theta) {
    tau.numerical(theta, frailty) - tau
  }, c(lb.frailty[[frailty]], min(ub.frailty[[frailty]], 100)))$root
} 