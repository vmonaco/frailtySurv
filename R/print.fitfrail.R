print.fitfrail <- function(fit, digits=max(options()$digits - 3, 3), ...) {
  
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Call: ")
  dput(fit$call)
  
  cat("\nCoefficients")
  tmp <- cbind(fit$beta)
  dimnames(tmp) <- list(names(fit$beta), "")
  print(tmp, right=TRUE, print.gap = 32)
  
  cat("\n")
  cat("Frailty distribution parameters   ", 
      toString(format(round(fit$theta, digits), digits=digits)), 
      " (", fit$frailty, ")\n",
      sep="")
  cat("Frailty distribution variance     ", 
      format(round(fit$frailty.variance, digits), digits=digits), 
      "\n",
      sep="")
  cat("Log-likelihood                    ", 
      format(round(fit$loglik, digits) , digits=digits), 
      "\n", 
      sep="")
  cat("Converged (method)                ",
      fit$iter,
      " iterations ",
      ifelse(fit$fitmethod == "score", "(solved score equations)", "(maximized log-likelihood)"),
      "\n",
      sep="")
  
  invisible()
}