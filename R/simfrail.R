# Generate data and fit a model many times
simfrail <- function(reps, 
                     genfrail.args, 
                     fitfrail.args,
                     base.time, # where to evaluate the baseline hazard function
                     mc.cores=0 # 0 to use all available cores, -1 to use all but 1, etc
                     ) { 
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    data <- do.call("genfrail", as.list(genfrail.args))
    fitfrail.args[["data"]] <- data
    
    # Time the fitter
    runtime <- system.time(fit <- do.call("fitfrail", as.list(fitfrail.args)))["elapsed"]
    names(runtime) <- "runtime"
    
    cens <- 1 - sum(data$status)/length(data$status)
    names(cens) <- "cens"
    
    beta <- fit$beta
    names(beta) <- paste("beta", 1:length(fit$beta))
    
    theta=fit$theta
    names(theta) <- paste("theta", 1:length(fit$theta))
    
    lambda <- fit$lambdafn(base.time)
    names(lambda) <- paste("lambda.", base.time, sep="")
    
    Lambda <- fit$Lambdafn(base.time)
    names(Lambda) <- paste("Lambda.", base.time, sep="")
    
    c(runtime,
      cens,
      beta,
      theta,
      lambda,
      Lambda)
  }
  
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (mc.cores != 1 && have.parallel) {
    # Run in parallel and make a row-observation matrix
    if (mc.cores <= 0) {
      mc.cores <- parallel::detectCores() + mc.cores
    }
    results <- parallel::mclapply(seeds, fn, mc.preschedule=FALSE, 
                        mc.set.seed=TRUE, mc.cores=mc.cores, mc.silent=FALSE)
  } else if (use.parallel && !have.parallel) {
    # Run in serial with a warning
    warning('Cannot run simfrail in parallel: requires the "parallel" package.')
    results <- lapply(seeds, fn)
  } else {
    # Run in serial
    results <- lapply(seeds, fn)
  }
  
  results <- t(simplify2array(results))
  data.frame(results)
}

# Perform simfrail for multiple param values to passed genfrail
simfrail.multigen <- function(reps, seed, genfrail.args, fitfrail.args, base.time,
                              param.name, param.values) {
  results <- lapply(param.values, function(pvalue) {
    genfrail.args[[param.name]] <- pvalue
    set.seed(seed) # seed before each run
    partial.results <- simfrail(reps, genfrail.args, fitfrail.args, base.time)
    # Name the first column as the parameter name
    setNames(cbind(param.name=toString(pvalue), partial.results), 
             c(param.name, names(partial.results)))
  })
  
  do.call("rbind", results)
}

plot.residuals <- function(results, true.values, title="") {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("plot.residuals requires the reshape2 package")
  }
  
  residuals <- cbind(N=results[,"N"], results[,names(true.values)] - true.values)
  
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

plot.baseline_hazard <- function(results, lambdafn=NULL, Lambdafn=NULL) {
  if (!is.null(lambdafn)) {
    bh <- results[,grepl("lambda.", names(results))]
    time_steps <- sapply(names(bh), 
                        function(w) gsub("lambda\\.", "", w),
                        USE.NAMES=FALSE)
    time_steps <- as.numeric(time_steps)
    names(bh) <- time_steps
    meltbh <- reshape2::melt(t(bh))
    names(meltbh) <- c("Time","instance", "value")
    bhfn <- lambdafn
  } else if (!is.null(Lambdafn)) {
    bh <- results[,grepl("Lambda.", names(results))]
    time_steps <- sapply(names(bh), 
                         function(w) gsub("Lambda\\.", "", w),
                         USE.NAMES=FALSE)
    time_steps <- as.numeric(time_steps)
    names(bh) <- time_steps
    meltbh <- reshape2::melt(t(bh))
    names(meltbh) <- c("Time","instance", "value")
    bhfn <- Lambdafn
  } else {
    stop("Must provide BH function lambdafn or cumulative BH function Lambdafn.")
  }
  
  truebh <- data.frame(x=time_steps, y=bhfn(time_steps))
  meltbh$type <- "Estimated"
  truebh$type <- "True"
  
  ggplot(meltbh,aes(x=Time,y=value,color=type)) +
    stat_summary(fun.data = "mean_cl_boot", geom = "smooth") +
    geom_line(aes(x=x, y=y, color=type), truebh) +
    scale_colour_manual("Cumulative BH",
                        values=c("True"="black","Estimated"="blue")) +
    ylab("Cumulative baseline hazard") +
    ggtitle("todo")
}

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrailcoxph <- function(reps, genfrail.args, coxph.args, seed=2015) 
{ 
  # Seed the beginning of each simulation so that each run is independent
  set.seed(seed)
  
  results <- replicate(reps, (function(){
    # Generate new random data, give this to coxph
    data <- do.call("genfrail", as.list(genfrail.args))
    coxph.args[["data"]] = data
    coxph.args[["method"]] = "breslow"
    
    # Only interested in total elapsed time
    runtime <- system.time(m <- do.call("coxph", as.list(coxph.args)))["elapsed"]
    names(runtime) <- "runtime"
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    c(m$coefficients,
      theta=m$history[[1]]$theta,
      censoring=1-sum(data$status)/length(data$status),
      runtime)
  })())
  
  # Summarize everything
  results <- apply(results, 1, function(x) c(mean=mean(x), sd=sd(x)))
  
  results
}

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simcoxph <- function(reps, genfrail.args, coxph.args, lambda_points,
                     mc.cores=parallel::detectCores(), seed=2015) { 
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
  set.seed(seed)
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    data <- do.call("genfrail", as.list(genfrail.args))
    # data$time <- data$time + rnorm(length(data$time), 0, 0.01)
    coxph.args[["data"]] = data
    coxph.args[["method"]] = "breslow"
    
    # Only interested in total elapsed time
    runtime <- system.time(fit <- do.call("coxph", as.list(coxph.args)))["elapsed"]
    
    names(runtime) <- "runtime"
    cens <- 1 - sum(data$status)/length(data$status)
    
    beta <- fit$coefficients
    names(beta) <- paste("beta", 1:length(beta))
    
    theta <- fit$history[[1]]$theta
    names(theta) <- paste("theta", 1:length(theta))
    
    bh <- basehaz(fit, centered=FALSE)
    Lambda_hat <- c(0, bh$hazard)
    time_steps <- c(0, bh$time)
    
    Lambdafn <- Vectorize(function(t) {
      if (t <= 0) {
        return(0);
      }
      Lambda_hat[sum(t > time_steps)]
    })
    
    Lambda <- Lambdafn(lambda_points)
    names(Lambda) <- paste("Lambda.", lambda_points, sep="")
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    c(runtime,
      cens=cens,
      beta,
      # theta,
      Lambda)
  }
  
  # Run in parallel and make a row-observation matrix
  results <- parallel::mclapply(seeds, fn, mc.preschedule=FALSE, 
                      mc.set.seed=TRUE, mc.cores=mc.cores, mc.silent=FALSE)
  # results <- lapply(seeds, fn)
  results <- t(simplify2array(results))
  
  # TODO: Summarize everything somewhere else:
  # results <- apply(results, 2, function(x) c(mean=mean(x), sd=sd(x)))
  
  data.frame(results)
}