\name{genfrail.control}
\alias{genfrail.control}

\title{
Control parameters for genfrail
}
\description{
This function creates a list of control parameters needed by genfrail. 
}

\usage{
genfrail.control(censor.reltol = .Machine$double.eps^0.25,
                 censor.subdivisions = 1000L,
                 crowther.reltol = .Machine$double.eps^0.25,
                 crowther.subdivisions = 1000L)
}

\arguments{
  \item{censor.reltol}{numeric relative tolerence for convergence of the censorship numerical integration. Default is .Machine$double.eps^0.25 (about 1.2e-4 on most machines).
  }
  \item{censor.subdivisions}{integer, maximum number of censorship numerical integration subdivisions. Default is 1000.
  }
  \item{crowther.reltol}{numeric relative tolerence for convergence of the numerical integration in Crowther's formula. Default is .Machine$double.eps^0.25 (about 1.2e-4 on most machines).
  }
  \item{crowther.subdivisions}{integer, maximum number of numerical integration subdivisions in Crowther's formula. Default is 1000.
  }
}

\value{
A list of control parameters.
}

\author{
John V. Monaco, Malka Gorfine, Li Hsu
}

\seealso{\code{\link{genfrail}}}