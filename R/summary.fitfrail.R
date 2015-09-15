summary.fitfrail <- function(object, Lambda.times=NULL, censored=FALSE, se=FALSE, CI=0.95, ...) {
  fit <- object
  
  if (!inherits(fit, "fitfrail")) 
    stop("summary.fitfrail can only be used for fitfrail objects")
  
  if (!is.null(Lambda.times)) {
    if (!is.numeric(Lambda.times)) 
      stop("Lambda.times must be numeric")
  } else if (censored) {
    Lambda.times <- fit$VARS$time
  } else {
    Lambda.times <- fit$VARS$time[fit$VARS$status > 0]
  }
  Lambda.times <- sort(unique(Lambda.times))
  
  status <- as.integer(fit$VARS$status > 0)
  n.risk.total <- length(status)
  
  n.risk <- n.risk.total - vapply(Lambda.times, function(t) {
    sum(fit$VARS$time < t)
  }, 0) # num failures at time t-
  
  n.event <- c(sum(status[fit$VARS$time <= Lambda.times[1]]), 
               diff(vapply(Lambda.times, 
                 function(t) {
                   sum(status[fit$VARS$time <= t])
                 }, 0))) # num failures at time t-
  
  surv <- exp(-fit$Lambda.fun(Lambda.times))
  
  result <- data.frame(time=Lambda.times,
                       surv=surv,
                       n.risk=n.risk,
                       n.event=n.event)
  rownames(result) <- NULL
  
  if (se) {
    cbh.se <- diag(vcov(fit, Lambda.times=Lambda.times, boot=TRUE, ...))
    cbh.se <- cbh.se[grepl("^Lambda", names(cbh.se))]
    std.err <- exp(-surv + cbh.se^2/2)*sqrt(exp(cbh.se^2) - 1)
    result$std.err <- std.err
    
    if (CI > 0) {
      zval <- qnorm(1- (1-CI)/2, 0,1)
      lower <- pmax(surv - zval* std.err, 0)
      upper <- pmin(surv + zval* std.err, 1)
      result$lower.ci <- lower
      result$upper.ci <- upper
    }
  }
  
  result
}