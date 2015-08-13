\name{genfrail}
\alias{genfrail}
\title{
Generate Survival Data
}
\description{
Generate survival data from a shared frailty model, with hazard function given by

\deqn{S(t)=exp[-\Lambda_0(t) \omega_i exp(\betaZ_ij)]}

where \eqn{\Lambda_0} is the cumulative baseline hazard, \eqn{\omega_i} is the frailty value for cluster \eqn{i}, \eqn{\beta} is the regression coefficient vector, and \eqn{Z_ij} is the covariate vector for individual \eqn{i} in cluster \eqn{j}.

The baseline hazard can be specified by the inverse cumualative baseline hazard, cumulative baseline hazard, or simply the baseline hazard. Frailty values are sampled from 

}
\usage{
genfrail <- function(N = 300, K = 2, 
                     beta = c(log(2),log(3)),
                     frailty = "gamma", theta = 2, 
                     covar.distr = "normal", 
                     censor.rate = 0.30,
                     Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{integer the number of clusters}
  
  \item{K}{integer, string, or vector. If an integer, the number of members in each cluster. If a string, the name of the distribution to sample the cluster sizes from. This can be one of: "poisson", "pareto", or "uniform". The \code{K.param} argument specifies the distribution parameters. If a vector, must be of length N and contains the integer size of each cluster.}
  \item{K.param}{vector of the cluster size distribution parameters if \code{K} is a string. If "possion", the vector should contain the rate and truncated value (see \code{rtpois}). If "pareto", the exponent, lower, and upper bounds (see \code{rtzeta}). If "uniform", the lower (noninclusive) and upper (inclusive) bounds.}
  
  \item{beta}{vector of regression coefficients.}

  \item{frailty}{string the name of the frailty distribution. Can be one of: "gamma", "pvf", "lognormal", "invgauss", "posstab", or "none". See \code{dgamma_r},\code{dpvf_r}, \code{dlognormal_r}, \code{dinvgauss_r}, \code{posstab_r} for the respective density functions. (Also see the *_c for C implementations of the respective density functions.) }
  \item{theta}{vector the frailty distribution parameters}

  \item{covar.distr}{string the distribution to sample covariates from. Can be one of: "normal", "uniform", "zero"}
  \item{covar.param}{vector the covariate distribution parameters.}
  \item{covar.matrix}{matrix with dimensions \code{c(NK, length(beta))} that contains the desired covariates. If not NULL, this overrides \code{covar.distr} and \code{covar.param}.}
  
  \item{censor.distr}{string the censoring distribution to use. Followup times are sampled from the censoring distribution to simulate non-informative right censorship. The censoring distribution can be one of: "normal", "lognormal", "none".}
  \item{censor.mu}{integer the mean of the censoring distribution (log-mean if \code{censor.distr} is "lognormal")}
  \item{censor.sigma}{integer the standard deviation of the censoring distribution (log-standard deviation if \code{censor.distr} is "lognormal"}
  \item{censor.rate}{numeric value between 0 and 1 to specify the empirical censoring rate. The \code{censor.mu} parameter is adjusted to achieve a desired censoring rate if \code{censor.rate} is specified. This only overrides \code{censor.mu}, and \code{censor.sigma} must still be specified so that the problem is identifiable.}

  \item{lambda_0}{function the baseline hazard function. Only one of \code{lambda_0}, \code{Lambda_0}, and \code{Lambda_0_inv} need to be specified. Passing the baseline hazard (\code{lambda_0}) is the most computationally expensive since this requires numeric integration inside a root-finding algorithm.}
  \item{Lambda_0}{function the cumulative baseline hazard function. This overrides \code{lambda_0}.}
  \item{lambda_0_inv}{function the inverse cumulative baseline hazard function. This overrides both \code{lambda_0} and \code{Lambda_0}.}
  
  \item{round.base}{numeric if specified, round the followup times to the nearest \code{round.base}}
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

}

\keyword{ survival }
\keyword{ shared frailty }
