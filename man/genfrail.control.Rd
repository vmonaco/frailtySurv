\name{genfrail.control}
\alias{genfrail.control}

\title{
Control parameters for genfrail
}
\description{
This function creates a list of control parameters needed by genfrail. 
}

\usage{
genfrail.control(censor.reltol = 1e-4,
                 censor.subdivisions = 1000L,
                 crowther.reltol = 1e-4,
                 crowther.subdivisions = 1000L)
}

\arguments{
  \item{censor.reltol}{numeric relative tolerence for convergence of the censorship numerical integration. Default is 0.001.
  }
  \item{censor.subdivisions}{integer, maximum number of censorship numerical integration subdivisions. Default is 1000.
  }
  \item{crowther.reltol}{numeric relative tolerence for convergence of the numerical integration in Crowther's formula. Default is 0.001.
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