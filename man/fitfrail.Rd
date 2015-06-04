\name{fitfrail}
\alias{fitfrail}
\title{Fit a Shared Frailty model}
\usage{
fitfrail(formula, data, control, 
         frailty = c("gamma", "lognormal", "posstab"), ...)
}
\arguments{
\item{formula}{a \code{formula} object, where the lhs is the response as a
\code{\link[survival]{Surv}} object and the terms are on the right.}

\item{data}{a data.frame that provides context for the formula}

\item{control}{control parameters in the form of a \code{\link{fitfrail.control}} object}

\item{frailty}{name of the shared frailty distribution}

\item{...}{additional arguments will be passed to fitfrail.control}
}
\value{
a fitfrail object representing the shared frailty model
}
\description{
\code{fitfrail} fits a model.
}
\details{
This function ...
}
\examples{
% \dontrun{
# 
# Generate synthetic survival data with regression coefficients 
# beta = c(log(2),log(3)) and theta = 2, where the shared frailty 
# values from a gamma distribution with expectation 1 and variance theta.
# 
simdata <- genfrail(beta=c(log(2),log(3)), N=30, K=2, 
                    theta=2, frailty="gamma", 
                    censor.mu=130, covariates="uniform")

# Fit a shared frailty model
fit.frailtyr <- fitfrail(Surv(time, status) ~ Z1 + cluster(family), 
                         data=simdata, frailty="gamma")
                         
% } % dontrun
} % examples