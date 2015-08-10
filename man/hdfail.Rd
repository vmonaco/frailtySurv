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
    \item{\code{serial}}{Unique serial number of the hard drive}
    \item{\code{model}}{a factor with levels}
    \item{\code{time}}{a numeric vector}
    \item{\code{status}}{a numeric vector}
    \item{\code{temp}}{a numeric vector}
    \item{\code{rsc}}{a numeric vector}
    \item{\code{rer}}{a numeric vector}
    \item{\code{psc}}{a numeric vector}
  }
}
\source{
\url{https://www.backblaze.com/hard-drive-test-data.html}
}
\examples{
data(hdfail)
## maybe str(hdfail) ; plot(hdfail) ...
}
\keyword{datasets}
