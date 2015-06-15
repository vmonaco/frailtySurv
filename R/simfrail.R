# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrail <- function(reps, genfrail.args, fitfrail.args, 
                     mc.cores=detectCores(), seed=2015) { 
  
  # To make the result reproducable when parallelized, seed each run from 
  # a list of random seeds sampled from a meta seed
  # set.seed(seed)
  seeds <- sample(1:1e7, reps, replace=TRUE)
  
  fn <- function(s) {
    set.seed(s)
    
    # Generate new random data, give this to coxph
    data <- do.call("genfrail", as.list(genfrail.args))
    fitfrail.args[["data"]] = data
    
    # Only interested in total elapsed time
    runtime <- system.time(fit <- do.call("fitfrail", as.list(fitfrail.args)))["elapsed"]
    names(runtime) <- "runtime"
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    c(beta=fit$beta,
      theta=fit$theta,
      cens=1-sum(data$status)/length(data$status),
      runtime)
  }
  
  # Run in parallel and make a row-observation matrix
  results <- mclapply(seeds, fn, mc.preschedule=FALSE, 
                      mc.set.seed=TRUE, mc.cores=mc.cores, mc.silent=FALSE)
  results <- t(simplify2array(results))
  
  # Summarize everything
  # results <- apply(results, 2, function(x) c(mean=mean(x), sd=sd(x)))
  
  data.frame(results)
}

fitfrail_residuals <- function(reps, beta, theta, formula, frailty, N=c(50,100,500,1000)) {
  
  res <- lapply(N, function(Ni) {
    r <- simfrail(reps, list(beta=beta, censor.mu=130, frailty=frailty, N=Ni, K=2, 
                         theta=theta, covariates="uniform"), 
             list(formula=formula, frailty=frailty, verbose=FALSE))
    
    cbind(N=Ni, r, 
               res.theta=r[,"theta"] - theta,
               res.beta=unname(data.matrix(r[, grepl("beta", names(r))])) - 
                          matrix(rep(beta, each=reps), nrow=reps)
          ) })
  
  do.call("rbind", res)
}

plot_residuals <- function(res, title) {
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