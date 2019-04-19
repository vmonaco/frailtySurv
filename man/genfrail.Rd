\name{genfrail}
\alias{genfrail}

\title{
Generate survival data
}

\description{
Generate clustered survival data from a shared frailty model, with hazard function given by

\deqn{S(t)=\exp [-\Lambda_0(t) \omega_i \exp (\beta Z_{ij})]}

where \eqn{\Lambda_0} is the cumulative baseline hazard, \eqn{\omega_i} is the frailty value of cluster \eqn{i}, \eqn{\beta} is the regression coefficient vector, and \eqn{Z_ij} is the covariate vector for individual \eqn{i} in cluster \eqn{j}.

The baseline hazard can be specified by the inverse cumualative baseline hazard, cumulative baseline hazard, or simply the baseline hazard. Frailty values can be sampled from gamma, power variance function (PVF), log-normal, inverse Gaussian, and positive stable distributions.

}

\usage{
genfrail(N = 300, K = 2, K.param = c(2, 0), beta = c(log(2)),
         frailty = "gamma", theta = c(2), 
         covar.distr = "normal", covar.param = c(0, 1), covar.matrix = NULL,
         censor.distr = "normal", censor.param = c(130, 15), 
         censor.rate = NULL, censor.time = NULL,
         lambda_0 = NULL, Lambda_0 = NULL, Lambda_0_inv = NULL, 
         round.base = NULL, control, ...)
}

\arguments{
  \item{N}{integer; number of clusters}
  
  \item{K}{integer, string, or vector; If an integer, the number of members in each cluster. If a string, the name of the distribution to sample the cluster sizes from. This can be one of: "poisson", "pareto", or "uniform". The \code{K.param} argument specifies the distribution parameters. If a vector, must be of length N and contains the integer size of each cluster.}
  \item{K.param}{vector of the cluster size distribution parameters if \code{K} is a string. If "possion", the vector should contain the rate and truncated value (see \code{rtpois}). If "pareto", the exponent, lower, and upper bounds (see \code{rtzeta}). If "uniform", the lower (noninclusive) and upper (inclusive) bounds.}
  
  \item{beta}{vector of regression coefficients.}

  \item{frailty}{string name of the frailty distribution. Can be one of: "gamma", "pvf", "lognormal", "invgauss", "posstab", or "none". See \code{dgamma_r},\code{dpvf_r}, \code{dlognormal_r}, \code{dinvgauss_r}, \code{posstab_r} for the respective density functions. (Also see the *_c for C implementations of the respective density functions.) }
  \item{theta}{vector the frailty distribution parameters}

  \item{covar.distr}{string distribution to sample covariates from. Can be one of: "normal", "uniform", "zero"}
  \item{covar.param}{vector covariate distribution parameters.}
  \item{covar.matrix}{matrix with dimensions \code{c(NK, length(beta))} that contains the desired covariates. If not NULL, this overrides \code{covar.distr} and \code{covar.param}.}
  
  \item{censor.distr}{string censoring distribution to use. Followup times are sampled from the censoring distribution to simulate non-informative right censorship. The censoring distribution can be one of: "normal", "lognormal", "uniform", "none".}
  \item{censor.param}{vector of censoring distribution parameters. For normal and lognormal censorship, this should be c(mu,sigma) where mu is the mean and sigma is the standard deviation (Note: this is still the mean and standard deviation for lognormal). For uniform censorship, the vector \code{c(lower, upper)} should specify the lower and upper bounds.}
  \item{censor.rate}{numeric value between 0 and 1 to specify the empirical censoring rate. The mean specified in the \code{censor.param} parameter is adjusted to achieve a desired censoring rate if \code{censor.rate} is given. Note that the standard deviation (the second parameter in \code{censor.param}) must still be specified so that the problem is identifiable. For uniform censorship, the interval given by \code{c(lower, upper)} is adjusted to achieve the desired censorship, while keeping the variance fixed (i.e., upper - lower does not change).}
  \item{censor.time}{vector of right-censorship times. This must have length N*K and specifies the right-censoring times of each observation. Note that this overrides all other censor.* params and cannot be used with variable cluster sizes.}

  \item{lambda_0}{function baseline hazard. Only one of \code{lambda_0}, \code{Lambda_0}, and \code{Lambda_0_inv} need to be specified. Passing the baseline hazard (\code{lambda_0}) is the most computationally expensive since this requires numerical integration inside a root-finding algorithm.}
  \item{Lambda_0}{function cumulative baseline hazard. This overrides \code{lambda_0}.}
  \item{Lambda_0_inv}{function inverse cumulative baseline hazard. This overrides both \code{lambda_0} and \code{Lambda_0}.}
  
  \item{round.base}{numeric if specified, round the followup times to the nearest \code{round.base}}
  
  \item{control}{control parameters in the form of a \code{\link{genfrail.control}} object}
  
  \item{...}{additional arguments will be passed to \code{\link[frailtySurv]{genfrail.control}}}
}

\value{
A \code{data.frame} with row-observations is returned.
  \item{family}{the cluster}
  \item{rep}{the member within each cluster}
  \item{time}{observed followup time}
  \item{status}{failure indicator}
  \item{Z1...}{covariates, where there are \code{length(beta)} Z columns}
}

\author{
  John V. Monaco, Malka Gorfine, and Li Hsu.
}

\seealso{
\code{\link{fitfrail}}
}

\examples{
# Generate the same dataset 3 different ways

# Using the baseline hazard (least efficient)
set.seed(1234)
dat.1 <- genfrail(N = 300, K = 2, 
                  beta = c(log(2),log(3)),
                  frailty = "gamma", theta = 2,
                  lambda_0=function(t, tau=4.6, C=0.01) (tau*(C*t)^tau)/t)

# Using the cumulative baseline hazard
set.seed(1234)
dat.2 <- genfrail(N = 300, K = 2, 
                  beta = c(log(2),log(3)),
                  frailty = "gamma", theta = 2, 
                  Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau)

# Using the inverse cumulative baseline hazard (most efficient)
set.seed(1234)
dat.3 <- genfrail(N = 300, K = 2, 
                  beta = c(log(2),log(3)),
                  frailty = "gamma", theta = 2, 
                  Lambda_0_inv=function(t, tau=4.6, C=0.01) (t^(1/tau))/C)

# Generate data with PVF frailty, truncated Poisson cluster sizes, normal
# covariates, and 0.35 censorship from a lognormal distribution
set.seed(1234)
dat.4 <- genfrail(N = 100, K = "poisson", K.param=c(5, 1), 
                  beta = c(log(2),log(3)),
                  frailty = "pvf", theta = 0.3, 
                  covar.distr = "lognormal", 
                  censor.rate = 0.35) # Use the default baseline hazard

# Cluster sizes have size >= 2, summarized by
summary(dat.4)

# An oscillating baseline hazard
set.seed(1234)
dat.5 <- genfrail(lambda_0=function(t, tau=4.6, C=0.01, A=2, f=0.1) 
                              A^sin(f*pi*t) * (tau*(C*t)^tau)/t)

# Uniform censorship with 0.25 censoring rate
set.seed(1234)
dat.6 <- genfrail(N = 300, K = 2, 
                  beta = c(log(2),log(3)),
                  frailty = "gamma", theta = 2, 
                  censor.distr = "uniform", 
                  censor.param = c(50, 150), 
                  censor.rate = 0.25,
                  Lambda_0_inv=function(t, tau=4.6, C=0.01) (t^(1/tau))/C)
}

\keyword{ survival data }
\keyword{ shared frailty }
