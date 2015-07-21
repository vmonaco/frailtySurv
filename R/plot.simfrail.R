plot.simfrail <- function(sim, type=c("residuals","hazard"), ...) {
  if (type == "residuals") {
    plot.simfrail.residuals(sim, ...)
  } else if (type == "hazard") {
    plot.simfrail.hazard(sim, ...)
  }
}

plot.simfrail.residuals <- function(results, true.values, title="", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting residuals requires the ggplot2, and reshape2 packages")
  }
  require(reshape2)
  residuals <- cbind(N=results["N"],
                     t(apply(results[,names(true.values)], 1, 
                             function(x) x - true.values)))
  
  # Select N and residuals columns, start with res
  var.names <- names(true.values)
  n.vars <- length(true.values)
  
  residuals.melt <- reshape2::melt(residuals, id = c("N"))
  cases <- c(t(unique(residuals.melt["N"])))
  
  # bp_ylim <- 3 #0.25*mean_res
  boxplot(value ~ N + variable, data=residuals.melt, notch=TRUE, 
          col=rep(1:n.vars, each=length(cases)), 
          names=rep(cases, n.vars), 
          xlab="N", ylab="Residual", 
          # ylim=c(-bp_ylim, bp_ylim),
          main=title)
  abline(h=0)
  legend("topright", legend=names(true.values), fill=1:n.vars, ncol=n.vars)
}

plot.simfrail.hazard <- function(results, true.cbh=NULL, title=NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE) || 
      !requireNamespace("gridExtra", quietly = TRUE) ||
      !requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Plotting the hazard requires the ggplot2, gridExtra, reshape2, and Hmisc packages")
  }
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  plotter <- function(haz, ylabel, true.values=NULL) {
    base.time <- sapply(names(haz), function(w) gsub(".+\\.", "", w), USE.NAMES=FALSE)
    base.time <- as.numeric(base.time)
    names(haz) <- base.time
    melthaz <- reshape2::melt(t(haz))
    names(melthaz) <- c("Time","instance", "value")
    melthaz$type <- "Mean (95% CI)"
    
    if (!is.null(true.values)) {
      if (is.function(true.values)) {
        true.values <- true.values(base.time)
      } else if (length(base.time) != length(true.values)) {
        stop("Must provide a function or same number of BH points as evaluated in results.")
      }
      truebh <- data.frame(x=base.time, y=true.values)
      truebh$type <- "Actual"
    }
    
    p <- ggplot(melthaz, aes(x=Time,y=value,color=type)) +
      stat_summary(fun.data="mean_cl_boot", geom="smooth") +
      theme(legend.position=c(0,1),
            legend.justification=c(0,1)) +
      ylab(ylabel) + 
      ggtitle(title)
    
    if (!is.null(true.values)) {
      p <- p + 
        geom_line(aes(x=x, y=y, color=type), truebh) +
        scale_colour_manual("", values=c("black","blue"))
    } else {
      p <- p + scale_colour_manual("", values=c("blue"))
    }
    
    p
  }
  
  plotter(results[,grepl("Lambda", names(results))], 
          "Cumulative baseline hazard",
          true.cbh)
}