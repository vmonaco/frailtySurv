\docType{package}
\name{frailtyr}
\alias{frailtyr}
\alias{frailtyr-package}
\title{General semiparametric shared frailty model}
\description{
The R package frailtyr is for simulating and fitting semiparametric shared frailty models. It implements a relatively simple estimator with desirable asymptotic properties under several permissible technical conditions. The frailtyr package provides functions
\itemize{
\item{\code{genfrail} for data generation}
\item{\code{fitfrail} for model fitting}
\item{\code{simfrail} for survival simulation}
}
}

\details{
  \tabular{ll}{
     Package: \tab frailtyr\cr
     Type: \tab Package\cr
     Version: \tab 1.0.0\cr
     Date: \tab August 2015\cr
     License: \tab LGPL-2\cr
     LazyLoad: \tab Yes\cr
  }
}

\author{
  John V. Monaco, Malka Gorfine, and Li Hsu.
}

\references{
  Gorfine M, Zucker DM, Hsu L (2006)
  Prospective survival analysis with a general semiparametric shared frailty model: A pseudo full likelihood approach.
  \emph{Biometrika}, \bold{93}(3), 735-741.
  
  Zucker DM, Gorfine M, Hsu L (2008)
  Pseudo-full likelihood estimation for prospective survival analysis with a general semiparametric shared frailty model: Asymptotic theory.
  \emph{Journal of Statistical Planning and Inference}, \bold{138}(7), 1998-2016.
}

\seealso{
\code{\link{genfrail}}, \code{\link{fitfrail}}, \code{\link{simfrail}}
}