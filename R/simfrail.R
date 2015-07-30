#' Repeatedly simulate survival data and fit a model
#' 
#' Generate data and fit a model many times. 
#' Reproduceable sim can be obtained by simply setting the seed before calling simfrail.
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
                       paste("Lambda.", Lambda.time, sep=""))
    
    # Elapsed time of the fitter
    runtime <- system.time(fit <- do.call("fitfrail", as.list(fitfrail.args)))[["elapsed"]]
    
    # Empirical censoring rate
    cens <- 1 - sum(dat$status)/length(dat$status)
    
    # Gather the parameter estimates
    hat.beta <- setNames(fit$beta, 
                         paste("hat.beta.", 1:length(beta), sep=""))
    hat.theta <- setNames(fit$theta, 
                          paste("hat.theta.", 1:length(theta), sep=""))
    hat.Lambda <- setNames(fit$fun.Lambda(Lambda.time), 
                           paste("hat.Lambda.", Lambda.time, sep=""))
    
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
                           paste("se.Lambda.", Lambda.time, sep=""))
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
  
  sim <- data.frame(t(simplify2array(plapply(reps, fn, cores))))
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
  sim <- lapply(param.values, function(pvalue) {
    genfrail.args[[param.name]] <- pvalue
    set.seed(seed) # seed before each run
    partial.sim <- simfrail(reps, genfrail.args, fitfrail.args, Lambda.time, ...)
    # Name the first column as the parameter name
    setNames(cbind(param.name=toString(pvalue), partial.sim), 
             c(param.name, names(partial.sim)))
  })
  
  do.call("rbind", sim)
  class(sim) <- c("simfrail", "data.frame")
  
  sim
}

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrailcoxph <- function(reps, genfrail.args, coxph.args, seed=2015) 
{ 
  # Seed the beginning of each simulation so that each run is independent
  set.seed(seed)
  
  sim <- replicate(reps, (function(){
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
  sim <- apply(sim, 1, function(x) c(mean=mean(x), sd=sd(x)))
  
  sim
}

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simcoxph <- function(reps, 
                     genfrail.args, 
                     coxph.args,
                     Lambda.time, # where to evaluate the baseline hazard function
                     vcov.args=list(),
                     cores=0 # 0 to use all available cores, -1 to use all but 1, etc
) { 
  Call <- match.call()
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    # Generate new random data, give this to coxph
    dat <- do.call("genfrail", as.list(genfrail.args))
    coxph.args[["data"]] = dat
    
    # Obtain the true values for each parameter
    beta <- attr(dat, "beta")
    theta <- attr(dat, "theta")
    Lambda <- setNames(attr(dat, "Lambda_0")(Lambda.time),
                       paste("Lambda.", Lambda.time, sep=""))
    
    # Elapsed time of the fitter
    runtime <- system.time(fit <- do.call("coxph", as.list(coxph.args)))[["elapsed"]]
    
    cens <- 1 - sum(dat$status)/length(dat$status)
    
    hat.beta <- setNames(fit$coefficients, paste("hat.beta.", 1:length(beta), sep=""))
    
    hat.theta <- setNames(fit$history[[1]]$theta, paste("hat.theta.", 1:length(theta), sep=""))
    
    bh <- basehaz(fit, centered=FALSE)
    Lambda_hat <- c(0, bh$hazard)
    time_steps <- c(0, bh$time)
    
    fun.Lambda <- Vectorize(function(t) {
      if (t <= 0) {
        return(0);
      }
      Lambda_hat[sum(t > time_steps)]
    })
    
    hat.Lambda <- setNames(fun.Lambda(Lambda.time),
                           paste("hat.Lambda.", Lambda.time, sep=""))
    
    SE <- sqrt(diag(vcov(fit)))
    se.beta <- setNames(SE[1:length(beta)],
                        paste("se.beta.", 1:length(hat.beta), sep=""))
    
    se.theta <- setNames(rep(NA, length(hat.theta)),
                         paste("se.theta.", 1:length(hat.theta), sep=""))
    
    se.Lambda <- setNames(rep(NA, length(Lambda.time)),
                          paste("se.Lambda.", Lambda.time, sep=""))
    
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
      se.Lambda)
  }
  
  sim <- data.frame(t(simplify2array(plapply(reps, fn, cores))))
  class(sim) <- append("simfrail", class(sim))
  
  attributes(sim) <- append(attributes(sim), list(
    call=Call,
    reps=reps,
    frailty=unique(c(genfrail.args[["frailty"]]))
  ))
  
  sim
}