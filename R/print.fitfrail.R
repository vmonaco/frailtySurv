print.fitfrail <- function(fit, digits=max(options()$digits - 4, 3), ...) {
  
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Call:\n")
  dput(fit$call)
  cat("\n")
  
  cat("Coefficents:\n")
  cat(signif(fit$beta, digits))
  cat("\n")
  
  cat("Frailty distribution parameters:\n")
  cat(signif(fit$theta, digits))
  cat("\n")
  
  cat("Log likelihood:\n")
  cat(format(round(fit$loglik, digits), nsmall=2))
  cat("\n")
  
  invisible()
}