\name{fitfrail}
\alias{fitfrail}
\title{Fit Shared Frailty Model}
\description{
\code{fitfrail} fits a Cox proportional hazards model with hidden shared frailty.
}
\usage{
fitfrail(formula, dat, control)
}
\arguments{
\item{formula}{a \code{formula} object, where the lhs is the response as a
\code{\link[survival]{Surv}} object and the terms are on the right.}

\item{dat}{a data.frame that provides context for the formula}

\item{control}{control parameters in the form of a \code{\link{fitfrail.control}} object}

\item{frailty}{name of the shared frailty distribution}

\item{...}{additional arguments will be passed to fitfrail.control}
}
\value{
a fitfrail object representing the shared frailty model
}
\references{
  Chambers, J. M. (1992)
  \emph{Linear models.}
  Chapter 4 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  Wilkinson, G. N. and Rogers, C. E. (1973)
  Symbolic descriptions of factorial models for analysis of variance.
  \emph{Applied Statistics}, \bold{22}, 392--9.
}
\author{
  The design was inspired by the S function of the same name described
  in Chambers (1992).  The implementation of model formula by Ross Ihaka
  was based on Wilkinson & Rogers (1973).
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