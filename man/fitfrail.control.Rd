\name{fitfrail.control}
\alias{fitfrail.control}

\title{
Control parameters for fitfrail
}
\description{
This function creates a list of control parameters needed by fitfrail. 
}

\usage{
fitfrail.control(fitmethod = "loglik", abstol = 1e-08, reltol = 1e-06, 
                 maxit = 100, verbose = FALSE)
}

\arguments{
  \item{fitmethod}{string indicating which fit method should be used: "loglik" or "score". If "loglik", then the loglikelihood is maximized directily using optim (L-BFGS-B algorithm). If "score", then the system of normalized score equations is solved using nleqslv (Newton algorithm).
}
  \item{abstol}{numeric absolute tolerance for convergence. If fitmethod is "loglik", convergence is reached when the absolute reduction in loglikelihood is less than abstol. If fitmethod is "score", convergence is reached when the absolute value of each normalized score equation is less then abstol.
}
  \item{reltol}{numeric relative tolerance for convergence. If fitmethod is "loglik", convergence is reached when the relative reduction in loglikelihood is less than reltol. If fitmethod is "score", 
}
  \item{maxit}{integer, the maximum number of iterations to perform}
  
  \item{verbose}{logical value, whether to print a trace of the parameter estimation}
}

\value{
A list of control parameters.
}

\author{
John V. Monaco, Malka Gorfine, Li Hsu
}

\seealso{\code{\link{fitfrail}}}