\name{fitfrail.control}
\alias{fitfrail.control}

\title{
Control parameters for fitfrail
}
\description{
This function creates a list of control parameters needed by fitfrail. 
}

\usage{
fitfrail.control(fitmethod = "loglik", 
                 abstol = 0, reltol = 1e-6, maxit = 100,
                 int.abstol = 0, int.reltol = 1, int.maxit = 1000, 
                 init.beta="coxph", init.theta=NULL, verbose = FALSE)
}

\arguments{
  \item{fitmethod}{string indicating which fit method should be used: "loglik" or "score". If "loglik", then the loglikelihood is maximized directily using optim (L-BFGS-B algorithm). If "score", then the system of normalized score equations is solved using nleqslv (Newton algorithm).
  }
  \item{abstol}{numeric absolute tolerance for convergence. If fitmethod is "loglik", convergence is reached when the absolute reduction in loglikelihood is less than abstol. If fitmethod is "score", convergence is reached when the absolute value of each normalized score equation is less then abstol.
  }
  \item{reltol}{numeric relative tolerance for convergence. If fitmethod is "loglik", convergence is reached when the relative reduction in loglikelihood is less than reltol. If fitmethod is "score", convergence is reached when the relative reduction of each estimated parameter is less than reltol.
  }
  \item{maxit}{integer, maximum number of iterations to perform
  }
  \item{int.abstol}{numeric absolute tolerence for convergence of the numeric integration. Only applicable for frailty distributions that require numerical integration. If 0, then ignore. Default is 0.
  }
  \item{int.reltol}{numeric relative tolerence for convergence of the numeric integration. Only applicable for frailty distributions that require numerical integration. If 0, then ignore. Default is 1.
  }
  \item{int.maxit}{integer, maximum number of numeric integration function evaluations. Only applicable for frailty distributions that require numerical integration. If 0, then no limit. Default is 100.
  }
  \item{init.beta}{initial regression coefficients paramater estimates, can be numeric, string, or NULL. If NULL, then a zero vector is used. If "coxph" is passed, then regression coefficients are initialized to the estimates given by coxph (with or without frailty, depending on the init.theta parameter, see below). If numeric, then the initial regression coefficients are specified explicitly.
  }
  \item{init.theta}{initial frailty distribution parameter estimates, can be numeric, string, or NULL. If NULL, then a sensible starting value is chosen for the specified frailty distribution. If "coxph", then hat.theta is initialized such that the initial frailty distribution has the same rank correlation as the estimated gamma frailty distribution. This is achieved by assuming gamma frailty and estimating theta using coxph. This estimate is then transferred to the appropriate frailty distribution through Kendall's tau. This method only works well for weak cluster dependence, i.e., small Kendall's tau. If numeric, then the initial frailty distribution parameters are specified explicitly.
  }
  \item{verbose}{logical value, whether to print a trace of the parameter estimation.
  }
}

\value{
A list of control parameters.
}

\author{
John V. Monaco, Malka Gorfine, Li Hsu
}

\seealso{\code{\link{fitfrail}}}