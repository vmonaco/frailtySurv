fitfrail.control <- function(fitmethod="loglik",
                             abstol=0,
                             reltol=1e-6,
                             maxit=100,
                             int.abstol=0,
                             int.reltol=1, 
                             int.maxit=1000, 
                             init.beta="coxph",
                             init.theta=NULL,
                             verbose=FALSE
                             ) {
  
  if (!fitmethod %in% c("loglik", "score")) 
    stop("fitmethod must be one of: loglik, score")
  if (maxit < 0) 
    stop("Invalid value for iterations")
  if (abstol < 0) 
    stop ("Invalid absolute tolerance")
  if (reltol < 0) 
    stop ("Invalid relative tolerence")
  if (int.maxit < 0) 
    stop("Invalid value for integration iterations")
  if (int.abstol < 0) 
    stop ("Invalid absolute integration tolerance")
  if (int.reltol < 0) 
    stop ("Invalid relative integration tolerence")
  
  list(fitmethod=fitmethod,
       abstol=abstol,
       reltol=reltol,
       maxit=as.integer(maxit),
       int.abstol=int.abstol,
       int.reltol=int.reltol,
       int.maxit=as.integer(int.maxit),
       init.beta=init.beta,
       init.theta=init.theta,
       verbose=as.logical(verbose)
       )
}
