#' 
#' Initialize control parameters
#' 
#' @export
fitfrail.control <- function(fitmethod='score',
                             eps=1e-9, 
                             toler.chol=.Machine$double.eps^.75, 
                             iter.max=20,
                             toler.inf=sqrt(eps),
                             verbose=TRUE) {
  if (!fitmethod %in% c("score", "like")) stop("fitmethod must be one of: score, like")
  if (iter.max < 0) stop("Invalid value for iterations")
  if (eps <= 0) stop ("Invalid convergence criteria")
  list(fitmethod=fitmethod,
       eps=eps,
       iter.max=as.integer(iter.max),
       verbose=as.logical(verbose))
}
