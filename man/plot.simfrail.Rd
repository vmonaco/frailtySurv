\name{plot.simfrail}
\alias{plot.simfrail}

\title{
Plot method for \code{simfrail} objects
}
\description{
Plot the estimated parameter residuals or the mean estimated cumulative baseline hazard.
}

\usage{
\method{plot}{simfrail}(x, type = c("residuals", "cumhaz"), ...)
}

\arguments{
  \item{x}{{a \code{fitfrail} object}}
  
  \item{type}{string, the type of plot. Can be either "residuals" or "cumhaz". If \code{type="residuals"}, a boxplot of the estimated parameter residuals is created. If \code{type="cumhaz"}, the mean estimated and true cumulative baseline hazard are plotted.}
  
  \item{\dots}{extra arguments include: 
  
    \code{CI} for \code{type="cumhaz"}, the confidence interval for the empirical cumulative baseline hazard.
  
    \code{n.Lambda} for \code{type="residuals"}, the number of time points to show the cumulative baseline hazard residuals for.}
}

\value{
The plot object.
}

\author{
John. V Monaco, Malka Gorfine, Li Hsu
}

\seealso{
\code{\link{simfrail}}
}

\examples{
\dontrun{
set.seed(2015)
sim <- simfrail(1000,
    genfrail.args=alist(beta=c(log(2),log(3)), frailty="gamma", 
                        censor.rate=0.30, N=300, K=2, theta=2,
                        covar.distr="uniform", covar.param=c(0, 1),
                        Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau), 
    fitfrail.args=alist(formula=Surv(time, status) ~ Z1 + Z2 
                                                   + cluster(family), 
                        frailty="gamma"),
    Lambda.times=1:120)

# Make a boxplot of residuals
plot(sim, type="residuals")

# Plot the mean estimated cumulative baseline hazard and empirical 0.95 CI
plot(sim, type="cumhaz")
}
}