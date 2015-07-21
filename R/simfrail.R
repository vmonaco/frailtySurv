# Generate data and fit a model many times
# TODO: output progress
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
    
    Lambda <- fit$Lambdafn(base.time)
    names(Lambda) <- paste("Lambda.", base.time, sep="")
    
    COV <- fit$COV
    var.beta <- diag(COV)[1:length(beta)]
    names(var.beta) <- paste("var.beta", 1:length(fit$beta))
    var.theta <- diag(COV)[(length(beta)+1):(length(beta)+length(theta))]
    names(var.theta) <- paste("var.theta", 1:length(theta))
    
    sd.beta <- fit$SD[1:length(beta)]
    names(sd.beta) <- paste("sd.beta", 1:length(fit$beta))
    sd.theta <- fit$SD[(length(beta)+1):(length(beta)+length(theta))]
    names(sd.theta) <- paste("sd.theta", 1:length(theta))
    
    sd.Lambda <- rep(0, length(base.time))
    names(sd.Lambda) <- paste("sd.Lambda", 1:length(base.time))
    
    c(runtime,
      cens,
      beta,
      theta,
      sd.beta,
      sd.theta,
      sd.Lambda,
      var.beta,
      var.theta,
      Lambda
    )
  }
  
  use.parallel <- mc.cores != 1
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
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
  
  results <- data.frame(t(simplify2array(results)))
  class(results) <- "simfrail"
  
  results <- append(results, lapply(genfrail.args, eval))
  results <- append(results, lapply(fitfrail.args, eval))
  results$Lambda <- results$Lambda_0(base.time)
  
  results$reps <- reps
  
  results
}

# Perform simfrail for multiple param values to passed genfrail
simfrail.multigen <- function(reps, seed, genfrail.args, fitfrail.args, base.time,
                              param.name, param.values, ...) {
  results <- lapply(param.values, function(pvalue) {
    genfrail.args[[param.name]] <- pvalue
    set.seed(seed) # seed before each run
    partial.results <- simfrail(reps, genfrail.args, fitfrail.args, base.time, ...)
    # Name the first column as the parameter name
    setNames(cbind(param.name=toString(pvalue), partial.results), 
             c(param.name, names(partial.results)))
  })
  
  do.call("rbind", results)
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
simcoxph <- function(reps, 
                     genfrail.args, 
                     coxph.args, 
                     base.time,
                     mc.cores=0
) { 
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
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
    names(cens) <- "cens"
    
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
    
    Lambda <- Lambdafn(base.time)
    names(Lambda) <- paste("Lambda.", base.time, sep="")
    
    COV <- vcov(fit)
    var.beta <- diag(COV)[1:length(beta)]
    names(var.beta) <- paste("var.beta", 1:length(beta))
    #     var.theta <- diag(COV)[(length(beta)+1):(length(beta)+length(theta))]
    #     names(var.theta) <- paste("var.theta", 1:length(theta))
    
    SD <- sqrt(diag(COV))
    sd.beta <- SD[1:length(beta)]
    names(sd.beta) <- paste("sd.beta", 1:length(beta))
    #     sd.theta <- fit$SD[(length(beta)+1):(length(beta)+length(theta))]
    #     names(sd.theta) <- paste("sd.theta", 1:length(theta))
    
    c(runtime,
      cens,
      beta,
      theta,
      sd.beta,
      # sd.theta,
      var.beta,
      # var.theta,
      Lambda)
  }
  
  use.parallel <- mc.cores != 1
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
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