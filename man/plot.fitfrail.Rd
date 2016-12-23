\name{plot.fitfrail}
\alias{plot.fitfrail}

\title{
Plot method for \code{fitfrail} objects
}
\description{
Plot the cumulative baseline hazard estimates or the parameter trace from model estimation.
}
\usage{
\method{plot}{fitfrail}(x, type = c("cumhaz", "trace"), ...)
}

\arguments{
  \item{x}{{a \code{fitfrail} object}}
  
  \item{type}{string, the type of plot. Can be either "cumhaz" to plot the mean estimated cumulative hazard or "trace" to plot the paramater and log-likelihood trace.}
  
  \item{\dots}{extra arguments include:
  
    \code{CI} for \code{type="cumhaz"}, numeric confidence interval between 0 and 1. If \code{CI=0}, no confidence interval is displayed. Otherwise, the bootstrapped confidence interval is calculated and displayed.
  
    \code{end} for \code{type="cumhaz"}, numeric x-axis limit (plot up to time \code{end})
  
    \code{show.loglik} for \code{type="trace"}, logical whether to show the log-likelihood trace.
  }
  
}

\value{
The plot object.
}

\author{
John. V Monaco, Malka Gorfine, Li Hsu
}

\seealso{
\code{\link{fitfrail}}
}

\examples{
\dontrun{
data(drs)
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id), 
                    drs, frailty="gamma")

# Plot the parameter and log-likelihood trace
plot(fit.drs, type="trace") 

# This may take a while to run. 
# Use parameter B to specify the number of repetitions in the weighted bootstrap
plot(fit.drs, type="cumhaz", CI=0.95)
}
}