\docType{package}
\name{frailtySurv}
\alias{frailtySurv}
\alias{frailtySurv-package}
\title{General semiparametric shared frailty model}
\description{
\pkg{frailtySurv} provides a suite of functions for generating clustered survival data, fitting multivariate shared frailty models under a wide range of frailty distributions, and visualizing the output. The semi-parametric estimators have better asymptotic properties than most existing implementations, including consistent and asymptotically-normal estimators. Moreover, this is the first package that implements semi-parametric estimators with inverse Gaussian and PVF frailty models.

The \pkg{frailtySurv} package provides functions
\itemize{
\item{\code{\link{genfrail}} for data generation}
\item{\code{\link{fitfrail}} for model fitting}
\item{\code{\link{simfrail}} for survival simulation}
}
}

\details{
  \tabular{ll}{
     Package: \tab frailtySurv\cr
     Type: \tab Package\cr
     Version: \tab 1.2.0\cr
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