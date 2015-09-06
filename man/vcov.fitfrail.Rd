\name{vcov.fitfrail}
\alias{vcov.fitfrail}
\title{Compute variance/covariance matrix for fitfrail model}

\description{
Compute the variance/covariance matrix for fitfrail estimated parameters. This can be performed by a an asymptotically-normal and consistent variance estimator or a weighted bootstrap. The resulting covariance matrix is cached in the fitted object and later retrieved if the same arguments to \code{vcov.fitfrail} are supplied.
}

\usage{
\method{vcov}{fitfrail}(object, boot=FALSE, B=100, Lambda.times=NULL, cores=0, ...)
}

\arguments{
\item{object}{a \code{fitfrail} object}

\item{boot}{logical value, whether to use a weighted bootstrap. If boot == FALSE, a consistent
estimator is used and the cumulative baseline hazard variance will not be estimated.}

\item{B}{number of repetitions in the weighted bootstrap.}

\item{Lambda.times}{time points where the variance/covariance should be evaluated. If Lambda.times == NULL, then the points where the cumulative baseline hazard increases (where failures occur) are used.}

\item{cores}{number of cores to use when computing the covariance matrix in parallel}

\item{...}{extra arguments are not used}
}

\value{
variance/covariance matrix for the fitfrail model parameters
}

\seealso{
\code{\link{fitfrail}}
}

\examples{
\dontrun{
dat <- genfrail(N=200, K=2, beta=c(log(2),log(3)), 
                frailty="gamma", theta=2,
                censor.rate=0.35,
                Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau)

fit <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                dat, frailty="gamma")

# boot=TRUE will give the weighted bootstrap variance estimates
COV <- vcov(fit, boot=FALSE)
COV
}
}