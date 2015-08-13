\name{fitfrail}
\alias{fitfrail}
\title{Fit Shared Frailty Model}
\description{
Fit a Cox proportional hazards model with hidden shared frailty and unspecified baseline hazard. See Gorfine (2006) for details.
}
\usage{
fitfrail(formula, dat, control, frailty, weights, ...)
}
\arguments{

\item{formula}{a \code{formula} object, where the lhs is the response as a \code{\link[survival]{Surv}} object and rhs contain the terms, including a \code{cluster} term for the cluster identifier}

\item{dat}{data.frame that provides context for the formula}

\item{control}{control parameters in the form of a \code{\link{fitfrail.control}} object}

\item{frailty}{string name of the shared frailty distribution}

\item{weights}{vector of cluster weights}

\item{...}{additional arguments will be passed to \code{\link[frailtyr]{fitfrail.control}}}
}
\value{
A fitfrail object representing the shared frailty model.

  \item{beta}{the estimated regression coefficients}
  \item{theta}{the estimated frailty distribution parameters}
  \item{Lambda}{a \code{data.frame} with the estimated baseline hazard at each failure time}
  \item{fun.Lambda}{a function of time that returns the estimated baseline}
  \item{loglik}{the log-likelihood}
  \item{iter}{the number of iterations to convergence}
}
\references{
  Gorfine M, Zucker DM, Hsu L (2006)
  Prospective survival analysis with a general semiparametric shared frailty model: A pseudo full likelihood approach.
  \emph{Biometrika}, \bold{93}(3), 735-741.
}
\author{
  The estimation method was developed by Malka Gorfine, Li Hsu, and David Zucker. The implementation was done by John V. Monaco.
}

\examples{
# 
# Generate synthetic survival data with regression coefficients 
# beta = c(log(2),log(3)) and theta = 2, where the shared frailty 
# values from a gamma distribution with expectation 1 and variance theta.
# 
dat <- genfrail(N=300, K=2, beta=c(log(2),log(3)), 
                frailty="gamma", theta=2,
                censor.rate=0.35,
                Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau)

# Fit a shared frailty model
fit <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                data=dat, frailty="gamma")
fit

# Fit the DRS data, clustered on patient
data(drs)
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id), 
                    drs, frailty="gamma")
fit.drs

# Variance estimates
vcov(fit.drs)
}