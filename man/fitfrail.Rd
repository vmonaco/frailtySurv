\name{fitfrail}
\alias{fitfrail}
\title{Fit a shared frailty model}

\description{
Fit an extended Cox proportional hazards model with unobserved shared frailty variate and unspecified baseline hazard function, using a semiparametric estimation technique. See Gorfine et al.~(2006) and Zucker et al.~(2008) for details.
}

\usage{
fitfrail(formula, dat, control, frailty, weights = NULL, se = FALSE, ...)
}

\arguments{
\item{formula}{a \code{formula} object, where the lhs is the response as a \code{\link[survival]{Surv}} object and rhs contain the terms, including a \code{cluster} term for the cluster identifier}

\item{dat}{data.frame that provides context for the formula}

\item{control}{control parameters in the form of a \code{\link{fitfrail.control}} object}

\item{frailty}{string name of the shared frailty distribution}

\item{weights}{vector of cluster weights}

\item{se}{logical value, whether the standard errors of the regression coefficient and frailty distribution parameter estimates should be calculated. These are obtained using the \code{vcov.fitfrail} function}

\item{...}{additional arguments will be passed to \code{\link[frailtySurv]{fitfrail.control}}}
}

\value{
A fitfrail object representing the shared frailty model.

  \item{beta}{the estimated regression coefficients}
  \item{theta}{the estimated frailty distribution parameters}
  \item{Lambda}{a \code{data.frame} with the estimated baseline hazard at each failure time}
  \item{Lambda}{a \code{data.frame} with the estimated baseline hazard at all observed times}
  \item{Lambda.fun}{a function of time that returns the estimated baseline}
  \item{loglik}{the log-likelihood}
  \item{iter}{the number of iterations performed}
  \item{trace}{the parameter trace during estimation}
}

\section{Convergence}{
The initial values of the regression coefficients are provided by \code{coxph}. Convergence is reached when either the relative reduction or absolute reduction in loglikelihood or score equations (depending on the fitmethod used) are below a threshold. If the maxit iterations are performed before convergence, then the algorithm terminates with a warning.
}

\author{
  The estimation method was developed by Malka Gorfine, Li Hsu, and David Zucker; implemented by John V. Monaco.
}

\references{
  Gorfine M, Zucker DM, Hsu L (2006)
  Prospective survival analysis with a general semiparametric shared frailty model: A pseudo full likelihood approach.
  \emph{Biometrika}, \bold{93}(3), 735-741.
  
  Monaco JV, Gorfine M, Hsu L (2018)
  General Semiparametric Shared Frailty Model: Estimation and Simulation with {frailtySurv}
  \emph{Journal of Statistical Software}, \bold{86}(4), 1-42
  
  Zucker DM, Gorfine M, Hsu L (2008)
  Pseudo-full likelihood estimation for prospective survival analysis with a general semiparametric shared frailty model: Asymptotic theory.
  \emph{Journal of Statistical Planning and Inference}, \bold{138}(7), 1998-2016.
}

\seealso{
\code{\link{vcov.fitfrail}}, \code{\link{genfrail}}, \code{\link{simfrail}},
\code{\link[survival]{survfit}}, \code{\link[survival]{coxph}}
}

\examples{
\dontrun{
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
                dat, frailty="gamma")
fit

# The Lambda.fun function can give the estimated cumulative baseline hazard at
# any time
fit$Lambda.fun(seq(0, 100, by=10))

# Fit the DRS data, clustered on patient
data(drs)
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id), 
                    drs, frailty="gamma")
fit.drs
}

# 
# A small example with c(log(2),log(3)) coefficients, Gamma(2) frailty, and
# 0.10 censorship.
# 
dat <- genfrail(N=30, K=2, beta=c(log(2),log(3)), 
                frailty="gamma", theta=2,
                censor.rate=0.10,
                Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau)

# Fit a shared frailty model
fit <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                dat, frailty="gamma", se=TRUE)
fit

# Summarize the survival curve
head(summary(fit))
}