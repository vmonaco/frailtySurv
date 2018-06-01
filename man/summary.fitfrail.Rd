\name{summary.fitfrail}
\alias{summary.fitfrail}

\title{
Summary of the survival curve
}
\description{
Returns a \code{data.frame} summarizing the survival curve of the fitted model. If specified, this function uses a weighted bootstrap procedure to calculate SE of the survival curve estimates. Subsequente calls with the same arguments will use the cached SE and avoid performing the weighted bootstrap again.
}
\usage{
\method{summary}{fitfrail}(object, type = "survival", Lambda.times = NULL, 
                           censored = FALSE, se = FALSE, CI = 0.95, ...)
}

\arguments{
  \item{object}{{a \code{fitfrail} object}}
  
  \item{type}{string indicating the type of summary: either "survival" for a summary of the survival curve, or "cumhaz" for a summary of the cumulative baseline hazard.}
  
  \item{Lambda.times}{vector of times where the curve should be evaluated. The resulting \code{data.frame} will have 1 row for each time. If NULL and censored=TRUE, all observed times are used by default. If NULL and censored=FALSE, only the failure times are including in the results.
}

  \item{censored}{logical value, whether the survival curve should contain the censored times. Ignored if \code{Lambda.times} is not NULL.}
  
  \item{se}{logical value, whether the survival SE should be included with the results. If se=TRUE, a weighted bootstrap procedure is used to determine estimated survival SE.
}
  \item{CI}{numeric, the confidence interval to evaluate upper and lower limits for the survival estimate at each time point}
  
  \item{\dots}{extra arguments will be passed to \code{vcov.fitfrail}}
}

\value{
A \code{data.frame} summarizing the survival curve with the following columns.

  \item{time}{the time points}
  \item{surv/cumhaz}{survival/cumulative hazard estimate at time t+}
  \item{n.risk}{number of subjects at risk at time t-}
  \item{n.event}{the number of failures that occured from the last time point to time t+}
  \item{std.err}{the SE of the survival estimate}
  \item{lower.ci}{lower bound on the specified confidence interval}
  \item{upper.ci}{upper bound on the specified confidence interval}
}

\note{
Similar to \code{summary.survfit} function in the \code{survival} package. 
}

\seealso{
\code{\link{fitfrail}}, \code{\link{vcov.fitfrail}}
}
\examples{
\dontrun{
dat <- genfrail(N=200, K=2, beta=c(log(2),log(3)), 
                frailty="gamma", theta=2,
                censor.rate=0.35,
                Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau)

fit <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                dat, frailty="gamma")

surv <- summary(fitfrail, B=50, se=TRUE, CI=0.95)
head(surv)
}
}
