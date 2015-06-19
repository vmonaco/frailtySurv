# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrail <- function(reps, genfrail.args, fitfrail.args, lambda_points,
                     mc.cores=detectCores(), seed=2015, use.parallel=TRUE) { 
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
  set.seed(seed)
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    data <- do.call("genfrail", as.list(genfrail.args))
    # data$time <- data$time + rnorm(length(data$time), 0, 0.01)
    fitfrail.args[["data"]] = data
    
    # Only interested in total elapsed time
    runtime <- system.time(fit <- do.call("fitfrail", as.list(fitfrail.args)))["elapsed"]
    names(runtime) <- "runtime"
    
    cens <- 1 - sum(data$status)/length(data$status)
    
    beta <- fit$beta
    names(beta) <- paste("beta", 1:length(fit$beta))
    theta=fit$theta
    names(theta) <- paste("theta", 1:length(fit$theta))
    
    lambda <- fit$lambdafn(lambda_points)
    names(lambda) <- paste("lambda.", lambda_points, sep="")
    
    Lambda <- fit$Lambdafn(lambda_points)
    names(Lambda) <- paste("Lambda.", lambda_points, sep="")
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    c(runtime,
      cens=cens,
      beta,
      theta,
      lambda,
      Lambda)
  }
  
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
    # Run in parallel and make a row-observation matrix
    results <- mclapply(seeds, fn, mc.preschedule=FALSE, 
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
  
  # TODO: Summarize everything somewhere else:
  # results <- apply(results, 2, function(x) c(mean=mean(x), sd=sd(x)))
  
  data.frame(results)
}

fitfrail_residuals <- function(reps, beta, theta, formula, frailty, 
                               lambda_points=0:150, Lambdafn=function(t, tau=4.6, C=0.01) (C*t)^tau,
                               N=c(50,100)){#},500,1000)) {
  
  res <- lapply(N, function(Ni) {
    r <- simfrail(reps, list(beta=beta, censor.params=c(130,15), frailty=frailty, N=Ni, K=2, 
                         theta=theta, covariates="uniform"), 
             list(formula=formula, frailty=frailty, verbose=FALSE), lambda_points=lambda_points)
    
    cbind(N=Ni, r, 
               res.theta=r[,"theta.1"] - theta,
               res.beta=unname(data.matrix(r[, grepl("beta.", names(r))])) - 
                          matrix(rep(beta, each=reps), nrow=reps),
          res.Lambda.20=r$Lambda.20 - Lambdafn(20),
          res.Lambda.60=r$Lambda.60 - Lambdafn(60),
          res.Lambda.80=r$Lambda.80 - Lambdafn(80),
          res.Lambda.100=r$Lambda.100 - Lambdafn(100)
          ) })
  
  do.call("rbind", res)
}

plot_residuals <- function(res, title, lambda_points, lambdafn) {
  # Select N and residuals columns, start with res
  var_names <- sapply(names(res)[grepl("res", names(res))], 
                      function(w) gsub("res\\.", "", w),
                      USE.NAMES=FALSE)
  n_vars <- length(var_names)
  res.stacked <- melt(res[,grepl("N|res", names(res))], id = c("N"))
  cases <- c(t(unique(res.stacked["N"])))
  mean_res <- mean(c(t(abs(res.stacked["value"]))))
  bp_ylim <- 3 #0.25*mean_res
  boxplot(value~N+variable, data=res.stacked, notch=TRUE, 
          col=rep(1:n_vars, each=length(cases)), 
          names=rep(cases, n_vars), 
          xlab="N", ylab="Residual", 
          ylim=c(-bp_ylim, bp_ylim),
          main=title)
  abline(h=0)
  legend(x=0,y=bp_ylim, legend=var_names, fill=1:n_vars, ncol=n_vars)
}

plot_baseline_hazard <- function(results, lambdafn=NULL, Lambdafn=NULL) {
  if (!is.null(lambdafn)) {
    bh <- results[,grepl("lambda.", names(results))]
    time_steps <- sapply(names(bh), 
                        function(w) gsub("lambda\\.", "", w),
                        USE.NAMES=FALSE)
    time_steps <- as.numeric(time_steps)
    names(bh) <- time_steps
    meltbh <- melt(t(bh))
    names(meltbh) <- c("Time","instance", "value")
    bhfn <- lambdafn
  } else if (!is.null(Lambdafn)) {
    bh <- results[,grepl("Lambda.", names(results))]
    time_steps <- sapply(names(bh), 
                         function(w) gsub("Lambda\\.", "", w),
                         USE.NAMES=FALSE)
    time_steps <- as.numeric(time_steps)
    names(bh) <- time_steps
    meltbh <- melt(t(bh))
    names(meltbh) <- c("Time","instance", "value")
    bhfn <- Lambdafn
  } else {
    stop("Must provide BH function lambdafn or cumulative BH function Lambdafn.")
  }
  
  ggplot(meltbh,aes(x=Time,y=value)) +
    stat_summary(fun.data = "mean_cl_boot", geom = "smooth") +
    geom_line(aes(x=x, y=y), data.frame(x=time_steps, y=bhfn(time_steps)))
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
                     mc.cores=detectCores(), seed=2015) { 
  
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
    
    ss <- survfit(fit)
    
    Lambda_hat <- c(0, -log(ss$surv))
    time_steps <- c(0, ss$time)
    
    Lambdafn <- Vectorize(function(t) {
      if (t <= 0) {
        return(0);
      }
      Lambda_hat[sum(t > time_steps)]
      # Lambda_hat[which.min(abs(time_steps - t))]
    })
    
    Lambda <- Lambdafn(lambda_points)
    names(Lambda) <- paste("Lambda.", lambda_points, sep="")
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    c(runtime,
      cens=cens,
      beta,
      theta,
      Lambda)
  }
  
  # Run in parallel and make a row-observation matrix
  results <- mclapply(seeds, fn, mc.preschedule=FALSE, 
                      mc.set.seed=TRUE, mc.cores=mc.cores, mc.silent=FALSE)
  # results <- lapply(seeds, fn)
  results <- t(simplify2array(results))
  
  # TODO: Summarize everything somewhere else:
  # results <- apply(results, 2, function(x) c(mean=mean(x), sd=sd(x)))
  
  data.frame(results)
}