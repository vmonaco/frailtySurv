fitfrail.control <- function(fitmethod="loglik",
                             abstol=1e-8,
                             reltol=1e-6,
                             maxit=100,
                             verbose=FALSE) {
  
  if (!fitmethod %in% c("loglik", "score")) 
    stop("fitmethod must be one of: score, loglik")
  if (maxit < 0) 
    stop("Invalid value for iterations")
  if (abstol <= 0) 
    stop ("Invalid absolute tolerance")
  if (reltol <= 0) 
    stop ("Invalid relative tolerence")
  
  list(fitmethod=fitmethod,
       abstol=abstol,
       reltol=reltol,
       maxit=as.integer(maxit),
       verbose=as.logical(verbose))
}
