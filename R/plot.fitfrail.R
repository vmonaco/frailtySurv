plot.fitfrail <- function(fit, type="hazard", ...) {
  if (type == "hazard") {
    plot.fitfrail.hazard(fit, ...)
  } else if (type == "trace") {
    plot.fitfrail.fitter(fit, ...)
  }
}

plot.fitfrail.fitter <- function(fit, show.loglik=TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, reshape2 packages")
  }
  require(ggplot2)
  require(reshape2)
  
  trace <- fit$VARS$trace
  loglik <- trace[,c("Iteration", "loglik")]
  trace$loglik <- NULL
  
  melttrace <- reshape2::melt(trace, "Iteration", variable.name="Parameter")
  
  p <- ggplot(melttrace, aes(x=Iteration, y=value, color=Parameter)) +
    geom_line() +
    xlab("Iteration") + 
    ylab("Estimate") + 
    scale_x_continuous(breaks = 1:nrow(trace)) +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
    ggtitle("Parameter estimate trace")
  
  if (show.loglik) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Plotting requires the gridExtra package")
    }
    require(gridExtra)
    p2 <- ggplot(loglik, aes(x=Iteration, y=loglik)) +
      geom_line() + 
      xlab("Iteration") + 
      ylab("Logliklihood") + 
      scale_x_continuous(breaks = 1:nrow(trace)) + 
      ggtitle("Loglikelihood trace")
    
    gridExtra::grid.arrange(p, p2, ncol=2)
  } else {
    p
  }
  
}

plot.fitfrail.hazard <- function(fit, boot=TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package")
  }
  require(ggplot2)
  
  Lambda <- fit$Lambda
  end <- Lambda$time[nrow(Lambda)] + mean(diff(fit$Lambda$time))
  Lambda[nrow(Lambda)+1,] <- c(end, fit$Lambdafn(end))
  
  p <- qplot(time, Lambda, data = Lambda, geom="step") +
    geom_rug(sides="b", size=0.5) +
    xlab("Time") + 
    ylab("Cumulative baseline hazard") + 
    theme(legend.position="none") +
    ggtitle(attr(fit, "description"))
  
  if (boot) {
    COV <- vcov(fit, boot=TRUE, ...)
    se.Lambda <- sqrt(diag(COV))[(fit$VARS$n.gamma+1):nrow(COV)]
    
    LB <- fit$Lambda$Lambda - 1.96*se.Lambda
    UB <- fit$Lambda$Lambda + 1.96*se.Lambda
    
    time <- c(fit$Lambda$time[2:length(fit$Lambda$time)], 
              fit$Lambda$time[length(fit$Lambda$time)] + mean(diff(fit$Lambda$time)))
    
    LB <- data.frame(time=rep(time, rep(2, length(time)))[1:(2*length(time)-1)],
                     Lambda=rep(LB, rep(2, length(LB)))[2:(2*length(LB))])
    
    UB <- data.frame(time=rep(time, rep(2, length(time)))[1:(2*length(time)-1)],
                     Lambda=rep(UB, rep(2, length(UB)))[2:(2*length(UB))])
    
    df <- data.frame(px=c(LB$time, rev(UB$time)),
                     py=c(LB$Lambda, rev(UB$Lambda)))
    
    p <- p + geom_polygon(aes(x=px, y=py), data=df, fill="black", alpha=0.1)
  }
  
  p
}