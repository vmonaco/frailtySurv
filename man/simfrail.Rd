\name{simfrail}
\alias{simfrail}

\title{
Simulate survival data and fit models
}
\description{
Generates simulated clustered survival data by repeatedly generating data, using a shared frailty model, and fitting the models. Respective arguments are passed to genfrail and fitfrail, and the resulting parameter estimates are aggregated and summarized.
}
\usage{
simfrail(reps, genfrail.args, fitfrail.args, Lambda.times, 
         vcov.args = list(), cores = 0, skip.SE = FALSE)
}

\arguments{
  \item{reps}{number of times to repeat the simulation}
  \item{genfrail.args}{list of arguments to pass to \code{genfrail}}
  \item{fitfrail.args}{list of arguments to pass to \code{fitfrail}}
  \item{Lambda.times}{vector of time points to obtain baseline hazard estimates at}
  \item{vcov.args}{list of arguments to pass to vcov.fitfrail for variance estimates. This is mainly used to specify whether bootstrap or estimated variances should be obtained.}
  \item{cores}{integer; if > 0, the number of cores to use; if < 0, the number of cores not to use; if 0, use all available cores}
  \item{skip.SE}{logical value, whether to skip the standard error estimates (saves time)}
}

\value{
A \code{simfrail} object that is essentially a \code{data.frame} of the resulting parameter estimates. Each row is a single run, and columns are as follows.

  \item{seed}{the seed used for the run}
  \item{runtime}{the time it took to fit the model}
  
  \item{N}{number of clusters}
  \item{mean.K}{average cluster size}
  \item{cens}{empirical censorship}
  
  \item{beta}{true regression coefficients}
  \item{hat.beta}{estimated regression coefficients}
  \item{se.beta}{standard error of each regression coefficient}
  
  \item{theta}{true frailty distribution parameters}
  \item{hat.theta}{estimated frailty distribution parameters}
  \item{se.theta}{standard error of each frailty distribution parameter}
  
  \item{Lambda}{true cumulative baseline hazard at each Lambda.times point}
  \item{hat.Lambda}{estimated cumulative baseline hazard at each Lambda.times point}
  \item{se.Lambda}{standard error at each Lambda.times point}
  
}

\author{
John. V Monaco, Malka Gorfine, Li Hsu
}

\seealso{
\code{\link{genfrail}}, \code{\link{fitfrail}}
}

\examples{
\dontrun{
sim <- simfrail(reps=100,
                genfrail.args=alist(
                  N=50, K=2,
                  beta=c(log(2),log(3)), 
                  frailty="gamma", theta=2,
                  Lambda_0 = function(t, tau=4.6, C=0.01) (C*t)^tau), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                  frailty="gamma"),
                Lambda.times=1:120, cores = 0)

# Summarize the results
summary(sim)

# Plot the residuals
plot(sim, "residuals")
}
}

\keyword{ survival simulation }
