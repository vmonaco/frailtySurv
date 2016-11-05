fitfrail <- function(formula, dat, control, frailty, weights=NULL, se=FALSE, ...) {
  Call <- match.call()
  
  # Remove missing observations with a warning
  if (any(is.na(dat))) {
    prev.nrow <- nrow(dat)
    dat <- dat[complete.cases(dat),]
    new.nrow <- nrow(dat)
    warning("Removed ", prev.nrow - new.nrow, " rows with missing observations")
  }
  
  # create a call to model.frame() that contains the formula (required)
  # and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "dat"), names(Call), nomatch=0) 
  
  if (indx[1] == 0) stop("A formula argument is required")
  
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  special <- c("cluster")
  temp$formula <- if(missing(dat)) terms(formula, special) else terms(formula, special, dat=dat)
  
  mf <- eval(temp, parent.frame())
  if (nrow(mf) == 0) stop("No (non-missing) observations")
  Terms <- terms(mf)
  
  ## Match any ... args to fitfrail.control
  extraargs <- list(...)
  if (length(extraargs)) {
    controlargs <- names(formals(fitfrail.control)) #legal arg names
    indx <- pmatch(names(extraargs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(extraargs)[indx==0L]),
           domain = NA)
  }
  
  # Default to the fitfrail.control defaults
  if (missing(control)) control <- fitfrail.control(...)
  
  if (!match(frailty,c("gamma","lognormal","invgauss","pvf"), nomatch=0))
    stop("Unsupported frailty distribution:", frailty)
  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  cluster <- attr(Terms, "specials")$cluster
  if (!length(cluster))
    stop("Missing cluster, use coxph or coxme if there is no hidden frailty")
  
  # Drop cluster from response
  tempc <- untangle.specials(Terms, 'cluster', 1:10)
  ord <- attr(Terms, 'order')[tempc$terms]
  cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
  dropterms <- tempc$terms  # don't want cluster in the X matrix
  xlevels <- .getXlevels(Terms[-tempc$terms], mf)
  
  attr(Terms, "intercept") <- 1
  adrop <- 0  # levels of "assign" to be dropped; 0 = intercept
  
  temppred <- attr(terms, "predvars")
  Terms2 <- Terms[ -dropterms]
  if (!is.null(temppred)) {
    # subscripting a Terms object currently drops predvars, in error
    attr(Terms2, "predvars") <- temppred[-(1 + dropterms)] # "Call" object
  }
  
  X <- model.matrix(Terms2, mf, constrasts=NULL)
  # number the terms wrt the original model matrix
  # Do not forget the intercept, which will be a zero
  renumber <- match(colnames(attr(Terms2, "factors")), 
                    colnames(attr(Terms,  "factors")))
  attr(X, "assign") <- c(0, renumber)[1+attr(X, "assign")]
  
  Xatt <- attributes(X) 
  xdrop <- Xatt$assign %in% adrop  # columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  
  run.coxph <- FALSE
  if (is.null(control$init.beta)) {
    init.beta <- rep(0, ncol(X))
  } else if (is.numeric(control$init.beta)) {
    stopifnot(length(control$init.beta) == ncol(X))
    init.beta <- control$init.beta
  } else if (control$init.beta == "coxph") {
    run.coxph <- TRUE
  } else {
    stop("Invalid init.beta=", control$init.beta)
  }
  
  if (is.null(control$init.theta)) {
    init.theta <- init.frailty[[frailty]]
  } else if (is.numeric(control$init.theta)) {
    stopifnot(length(control$init.theta) == length(init.frailty[[frailty]]))
    if (any(control$init.theta < lb.frailty[[frailty]]) ||
        any(control$init.theta > ub.frailty[[frailty]])) {
      stop(sprintf("Frailty distribution initial parameters out of range. Parameters for %s frailty must be in: (%.2f, %.2f)", 
           frailty, lb.frailty[[frailty]], ub.frailty[[frailty]]))
    }
    init.theta <- control$init.theta
  } else if (control$init.theta == "coxph") {
    run.coxph <- TRUE
  } else {
    stop("Invalid init.theta=", control$init.theta)
  }
  
  # Initialize beta and theta from coxph estimates.
  # theta is initialized to the frailty distribution with same rank correlation 
  # as the gamma frailty estimated by coxph
  if (run.coxph) {
    # Drop the fitfrail cluster term and maybe add the coxph cluster term below
    coxph.formula <- drop.terms(terms(formula), tempc$terms, keep.response=TRUE)
    
    if (!is.null(control$init.theta) && control$init.theta == "coxph") {
      # Replace the cluster term with frailty.gamma
      # TODO: maybe a better way to specify frailty.gamma(cluster.term)
      cluster.term <- regmatches(tempc$vars, 
                                 gregexpr("(?<=\\().*?(?=\\))", tempc$vars, perl=T))[[1]]
      coxph.formula <- update(coxph.formula,
                              as.formula(sprintf("~ . + frailty.gamma(%s)", cluster.term)))
    } else {
      coxph.formula <- update(coxph.formula, as.formula("~ ."))
    }
    
    suppressWarnings(fit.coxph <- coxph(coxph.formula, data=dat))
    
    if (!is.null(control$init.beta) && control$init.beta == "coxph") {
      init.beta <- unname(fit.coxph$coefficients)
    }
    
    if (!is.null(control$init.theta) && control$init.theta == "coxph") {
      init.theta <- unname(theta.given.tau(tau.numerical(fit.coxph$history[[1]]$theta, "gamma"), frailty))
    }
  }
  
  names(init.beta) <- colnames(X)
  names(init.theta) <- 1

  fit <- fitfrail.fit(X, Y, cluster, 
                      init.beta, init.theta, 
                      frailty,
                      control, row.names(mf),
                      weights)
  class(fit) <- 'fitfrail'
  fit$call <- Call
  
  attr(fit, "description") <- paste("fitfrail: ", 
                                    fit$VARS$n.clusters, " clusters (avg. size ", 
                                    format(mean(fit$VARS$cluster_sizes), nsmall=2), "), ",
                                    toString(frailty), " frailty", sep="")
  
  if (se) {
    fit$vcov <- vcov(fit)
    fit$se <- diag(sqrt(fit$vcov))
    
    fit$se.beta <- fit$se[1:length(fit$beta)]
    fit$se.theta <- fit$se[(length(fit$beta)+1):(length(fit$beta)+length(fit$theta))]
  }
  
  fit
}