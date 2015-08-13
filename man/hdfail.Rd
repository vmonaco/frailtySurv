\name{hdfail}
\alias{hdfail}
\docType{data}
\title{
Hard drive failure dataset
}
\description{
This dataset contains the observed follow-up times and SMART statistics of 52k unique hard drives.

Daily snapshots of a large backup storage provider over 2 years were made publicly available. On each day, the Self-Monitoring, Analysis, and Reporting Technology (SMART) statistics of operational drives are recorded. When a hard drive is no longer operational, it is marked as a failure and removed from the subsequent daily snapshots. New hard drives are also continuously added to the population. In total, there are over 52k unique hard drives over approximately 2 years and 2885 (5.5\%) failures.

}
\usage{data("hdfail")}
\format{
  A data frame with 52422 observations on the following 8 variables.
  \describe{
    \item{\code{serial}}{unique serial number of the hard drive}
    \item{\code{model}}{hard drive model}
    \item{\code{time}}{the observed followup time}
    \item{\code{status}}{failure indicator}
    \item{\code{temp}}{temperature in Celsius}
    \item{\code{rsc}}{binary covariate, where 1 indicates sectors that encountered read, write, or verification errors}
    \item{\code{rer}}{binary covariate, where 1 indicates a non-zero rate of errors that occur in hardware when reading from data from disk.}
    \item{\code{psc}}{binary covariate, where 1 indicates there were sectors waiting to be remapped due to an unrecoverable error.}
  }
}

\source{
\url{https://www.backblaze.com/hard-drive-test-data.html}
}

\examples{
\dontrun{
data(hdfail)

# Select only Western Digital hard drives
dat <- subset(hdfail, grepl("WDC", model))

fit.hd <- fitfrail(Surv(time, status) ~ temp + rer + rsc 
                                      + psc + cluster(model), 
                   dat, frailty="gamma", fitmethod="score")

fit.hd
}
}
\keyword{datasets}
