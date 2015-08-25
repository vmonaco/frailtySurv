print.fitfrail <- function(x, digits=max(options()$digits - 4, 3), ...) {
  fit <- x
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Call: ")
  print(fit$call)
  
  cat("\n")
  tmp <- data.frame(Covariate=names(fit$beta), Coefficient=fit$beta)
  print(format(tmp, width=11, justify="left", nsmall=digits), print.gap=5, row.names=FALSE, right=FALSE)
  
  cat("\n")
  cat("Frailty distribution   ", 
      fit$frailty, 
      "(", toString(format(fit$theta, nsmall=digits)), "), ",
      "VAR = ", format(fit$frailty.variance, nsmall=digits), 
      "\n",
      sep="")
  cat("Log-likelihood         ", 
      format(fit$loglik, nsmall=digits), 
      "\n", 
      sep="")
  cat("Converged (method)     ",
      fit$iter,
      " iterations, ",
      format(fit$fit.time, nsmall=digits),
      " seconds ",
      ifelse(fit$fitmethod == "score", "(solved score equations)", "(maximized log-likelihood)"),
      "\n",
      sep="")
  
  invisible()
}