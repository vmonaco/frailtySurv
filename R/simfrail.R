#' Repeatedly simulate survival data and fit a model
#' 
#' Generate data and fit a model many times. 
#' Reproduceable results can be obtained by simply setting the seed before calling simfrail.
#' 
#' @param reps number of times to repeat the simulation
#' @param genfrail.args 
#' @param fitfrail.args
#' @param Lambda.time
#' @param cores
#' 
#' @return sim
#' 
#' @examples 
#' TODO
#' 
#' @export
simfrail <- function(reps, 
                     genfrail.args, 
                     fitfrail.args,
                     Lambda.time, # where to evaluate the baseline hazard function
                     vcov.args=list(),
                     cores=0, # 0 to use all available cores, -1 to use all but 1, etc
                     verbose=FALSE,
                     skip.SE=FALSE
) { 
  Call <- match.call()
  
  # To make the result reproducable when parallelized, seed each run 
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    dat <- do.call("genfrail", as.list(genfrail.args))
    fitfrail.args[["dat"]] <- dat
    
    # Obtain the true values for each parameter
    beta <- attr(dat, "beta")
    theta <- attr(dat, "theta")
    Lambda <- setNames(attr(dat, "Lambda_0")(Lambda.time),
                       paste("Lambda.", 1:length(Lambda.time), sep=""))
    
    # Time the fitter
    runtime <- system.time(fit <- do.call("fitfrail", as.list(fitfrail.args)))[["elapsed"]]
    
    # Empirical censoring rate
    cens <- 1 - sum(dat$status)/length(dat$status)
    
    # Gather the parameter estimates
    hat.beta <- setNames(fit$beta, 
                         paste("hat.beta.", 1:length(beta), sep=""))
    hat.theta <- setNames(fit$theta, 
                          paste("hat.theta.", 1:length(theta), sep=""))
    hat.Lambda <- setNames(fit$Lambdafn(Lambda.time), 
                           paste("hat.Lambda.", 1:length(Lambda.time), sep=""))
    
    if (!skip.SE) {
      # Gather the estimated standard errors
      V <- do.call(vcov, c(list(fit=fit, Lambda.time=Lambda.time, cores=1), vcov.args))
      
      SE <- sqrt(diag(V))
      
      se.beta <- setNames(SE[1:length(beta)],
                          paste("se.beta.", 1:length(beta), sep=""))
      se.theta <- setNames(SE[(length(beta)+1):(length(beta)+length(theta))],
                           paste("se.theta.", 1:length(theta), sep=""))
      se.Lambda <- setNames(SE[(length(beta)+length(theta)+1):
                                 ((length(beta)+length(theta))+length(Lambda))],
                           paste("se.Lambda.", 1:length(Lambda.time), sep=""))
    } else {
      se.beta <- NULL
      se.theta <- NULL
      se.Lambda <- NULL
    }
    
    c(seed=s,
      runtime=runtime,
      N=length(unique(dat$family)),
      mean.K=mean(table(dat$family)),
      cens=cens,
      beta,
      hat.beta,
      se.beta,
      theta,
      hat.theta,
      se.theta,
      Lambda,
      hat.Lambda,
      se.Lambda
    )
  }
  
  # TODO: output progress?
  use.parallel <- cores != 1
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
    # Run in parallel and make a row-observation matrix
    if (cores <= 0) {
      cores <- parallel::detectCores() + cores
    }
    sim <- parallel::mclapply(seeds, fn, mc.preschedule=FALSE, 
                              mc.set.seed=TRUE, mc.cores=cores, mc.silent=FALSE)
  } else if (use.parallel && !have.parallel) {
    # Run in serial with a warning
    warning('Cannot run simfrail in parallel: requires the "parallel" package.')
    sim <- lapply(seeds, fn)
  } else {
    # Run in serial
    sim <- lapply(seeds, fn)
  }
  
  sim <- data.frame(t(simplify2array(sim)))
  class(sim) <- append("simfrail", class(sim))
  
  attributes(sim) <- append(attributes(sim), list(
    call=Call,
    reps=reps,
    frailty=unique(c(genfrail.args[["frailty"]], fitfrail.args[["frailty"]]))
  ))
  
  sim
}

# Perform simfrail for multiple param values to passed genfrail
simfrail.multigen <- function(reps, seed, genfrail.args, fitfrail.args, Lambda.time,
                              param.name, param.values, ...) {
  results <- lapply(param.values, function(pvalue) {
    genfrail.args[[param.name]] <- pvalue
    set.seed(seed) # seed before each run
    partial.results <- simfrail(reps, genfrail.args, fitfrail.args, Lambda.time, ...)
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
                     Lambda.time,
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
    
    Lambda <- Lambdafn(Lambda.time)
    names(Lambda) <- paste("Lambda.", Lambda.time, sep="")
    
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