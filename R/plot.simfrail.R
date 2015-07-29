plot.simfrail <- function(sim, type=c("residuals","hazard"), ...) {
  if (type == "residuals") {
    plot.simfrail.residuals(sim, ...)
  } else if (type == "hazard") {
    plot.simfrail.hazard(sim, ...)
  }
}

plot.simfrail.residuals <- function(sim, n.Lambda=3, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, and reshape2 packages")
  }
  require(reshape2)
  
  Lambda.cols <- names(sim)[grepl("^Lambda", names(sim))]
  
  # All BUT n.Lambda
  if (n.Lambda < 0) {
    n.Lambda <- max(length(Lambda.cols) + n.Lambda, 0)
  } 
  
  if (n.Lambda == 0) {
    # No Lambda
    Lambda.cols <- NULL
  } else if (length(Lambda.cols) > n.Lambda) {
    # Evenly spaced n.Lambda
    idx <- round(seq(0, length(Lambda.cols), 
                     length.out=(n.Lambda+2))[2:(n.Lambda+1)])
    
    Lambda.cols <- Lambda.cols[idx]
  }
  hat.cols <- c(names(sim)[grepl("^hat.beta|^hat.theta", names(sim))], 
                paste("hat.", Lambda.cols, sep=""))
  value.cols <- c(names(sim)[grepl("^beta|^theta", names(sim))], Lambda.cols)
  
  residuals <- data.frame(cbind(N=sim$N,
                     vapply(value.cols,
    function(col) sim[[paste("hat.", col, sep="")]] - sim[[col]], rep(0, nrow(sim)))))
  
  # Select N and residuals columns, start with res
  n.vars <- length(value.cols)
  
  residuals.melt <- reshape2::melt(residuals, id = c("N"))
  cases <- c(t(unique(residuals.melt["N"])))
  
  # bp_ylim <- 3 #0.25*mean_res
  boxplot(value ~ N + variable, data=residuals.melt, notch=TRUE, 
          col=rep(1:n.vars, each=length(cases)), 
          names=rep(cases, n.vars), 
          xlab="N", ylab="Residual",
          main=attr(summary(sim), "description"))
  abline(h=0)
  legend("topright", legend=value.cols, fill=1:n.vars, ncol=n.vars)
}

plot.simfrail.hazard <- function(sim, title=NULL, CI=0.95, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, reshape2 packages")
  }
  require(ggplot2)
  require(reshape2)
  
  hats <- sim[,grepl("^hat.Lambda", names(sim))]
  se <- sim[,grepl("^se.Lambda", names(sim))]
  values <- colMeans(sim[,grepl("^Lambda", names(sim))])
  
  Lambda.time <- sapply(names(values), function(w) gsub(".+\\.", "", w), USE.NAMES=FALSE)
  Lambda.time <- as.numeric(Lambda.time)
  names(hats) <- Lambda.time
  melthats <- reshape2::melt(t(hats))
  names(melthats) <- c("Time","instance", "value")
  melthats$type <- "Empirical 95% CI"
  
  values <- data.frame(x=Lambda.time, y=values)
  values$type <- "Actual"
  
  Z.score <- qnorm((1-CI)/2)
  mean.se <- colMeans(se)
  se <- data.frame(x=Lambda.time, lower=values$y-Z.score*mean.se, upper=values$y+Z.score*mean.se)
  se$type <- "Estimated 95% CI"
  
  ggplot(melthats, aes(x=Time,y=value,color=type)) +
    stat_summary(fun.data="mean_cl_boot", geom="smooth") +
    geom_line(aes(x=x, y=y, color=type), values) +
    geom_line(aes(x=x, y=lower, color=type), se) +
    geom_line(aes(x=x, y=upper, color=type), se) +
    scale_colour_manual("Legend", values=c("black","blue","brown")) +
    theme(legend.position=c(0,1),
          legend.justification=c(0,1)) +
    ylab("Cumulative baseline hazard") + 
    ggtitle(attr(sim, "description"))
}