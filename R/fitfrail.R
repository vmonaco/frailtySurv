#' Fit a Shared Frailty model
#'
#' \code{fitfrail} fits a model.
#'
#' This function ...
#' 
#' @param ... model
#' @return parameters
#'   \url{} for more details.
#' @examples
#' fitfrail(genfrail())
#' 
#' @export
fitfrail <- function(formula, data, control, frailty=c("gamma","lognormal"), ...) {
  
  Call <- match.call()
  
  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "data"),
                names(Call), nomatch=0) 
  if (indx[1] ==0) stop("A formula argument is required")
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  special <- c("cluster")
  temp$formula <- if(missing(data)) terms(formula, special) else terms(formula, special, data=data)
  
  mf <- eval(temp, parent.frame())
  if (nrow(mf) ==0) stop("No (non-missing) observations")
  Terms <- terms(mf)
  
  ## Pass any ... args to nleqslv, but only valid args
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(nleqslv)) #legal arg names
    indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx==0L]),
           domain = NA)
  }
  if (missing(control)) control <- fitfrail.control(...)
  
  if (frailty== "gamma") {
    dfrailty <- dgamma_
    deriv_dfrailty <- deriv_dgamma
  } else {
    stop("frailty distribution", frailty, "not supported yet.")
  }
  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  cluster<- attr(Terms, "specials")$cluster
  if (length(cluster)) {
    tempc <- untangle.specials(Terms, 'cluster', 1:10)
    ord <- attr(Terms, 'order')[tempc$terms]
    cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
    dropterms <- tempc$terms  #we won't want this in the X matrix
    # Save away xlevels after removing cluster (we don't want to save upteen
    #  levels of that variable, which we will never need).
    xlevels <- .getXlevels(Terms[-tempc$terms], mf)
  } else {
    stop("Missing cluster, use coxme if there is no hidden frailty")
  }
  
  contrast.arg <- NULL  #due to shared code with model.matrix.coxph
  attr(Terms, "intercept") <- 1
  adrop <- 0  #levels of "assign" to be dropped; 0= intercept
  
  if (length(dropterms)) {
    temppred <- attr(terms, "predvars")
    Terms2 <- Terms[ -dropterms]
    if (!is.null(temppred)) {
      # subscripting a Terms object currently drops predvars, in error
      attr(Terms2, "predvars") <- temppred[-(1+dropterms)] # "Call" object
    }
    X <- model.matrix(Terms2, mf, constrasts=contrast.arg)
    # we want to number the terms wrt the original model matrix
    # Do not forget the intercept, which will be a zero
    renumber <- match(colnames(attr(Terms2, "factors")), 
                      colnames(attr(Terms,  "factors")))
    attr(X, "assign") <- c(0, renumber)[1+attr(X, "assign")]
  }
  else {
    X <- model.matrix(Terms, mf, contrasts=contrast.arg)
  } 
  
  Xatt <- attributes(X) 
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  
  beta_init <- coxph.fit(X, Y, strata=NULL, 
                          offset=NULL, init=NULL, 
                          control=coxph.control(), weights=NULL, 
                          method="efron", row.names(mf))$coefficients
  theta_init = c(0.1)
  fit <- sharedfrailty.fit(X, Y, cluster, 
                           beta_init, theta_init, 
                           dfrailty, deriv_dfrailty,
                           control, row.names(mf))
  class(fit) <- 'fitfrail'
  fit$call <- Call
  fit
}
