\name{hdfailure}
\alias{hdfailure}
\docType{data}
\title{
Hard drive failure dataset
}
\description{
This dataset contains the observed follow-up times and SMART statistics of 52k unique hard drives.

Daily snapshots of a large backup storage provider over 2 years were made publicly available. On each day, the Self-Monitoring, Analysis, and Reporting Technology (SMART) statistics of operational drives are recorded. When a hard drive is no longer operational, it is marked as a failure and removed from the subsequent daily snapshots. New hard drives are also continuously added to the population. In total, there are over 52k unique hard drives over approximately 2 years and 2885 (5.5\%) failures.
}
\usage{data("hdfailure")}
\format{
  A data frame with 52422 observations on the following 16 variables.
  \describe{
    \item{\code{serial_number}}{a factor with levels}
    \item{\code{day}}{a numeric vector}
    \item{\code{month}}{a numeric vector}
    \item{\code{dayofyear}}{a numeric vector}
    \item{\code{dayofweek}}{a numeric vector}
    \item{\code{weekofmonth}}{a numeric vector}
    \item{\code{temperature}}{a numeric vector}
    \item{\code{pending_sector_count}}{a numeric vector}
    \item{\code{read_error_rate}}{a numeric vector}
    \item{\code{reallocated_sector_count}}{a numeric vector}
    \item{\code{power_on_hours}}{a numeric vector}
  }
}
\source{
\url{https://www.backblaze.com/hard-drive-test-data.html}
}
\examples{
data(hdfailure)
}
\keyword{datasets}
