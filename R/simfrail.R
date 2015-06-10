# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrail <- function(reps, genfrail.args, fitfrail.args, 
                     mc.cores=detectCores()-1, seed=2015) { 
  
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

fitfrail_residuals <- function(reps, beta, theta, formula, N=c(50,60,70,80,90,100,500,1000)) {
  
  results <- lapply(N, function(Ni) {
    r <- simfrail(reps, list(beta=beta, censor.mu=130, frailty="gamma", N=Ni, K=2, 
                         theta=theta, covariates="uniform"), 
             list(formula=formula, frailty="gamma", verbose=FALSE))
    
    cbind(N=Ni, r, 
               res.theta=theta - r[,"theta"],
               res.beta=matrix(rep(beta, each=reps), nrow=reps) - 
            unname(data.matrix(r[, grepl("beta", names(r))]))) })
  
  do.call("rbind", results)
}

# r1.stacked = melt(r1[,c("N","res.theta","res.beta.1","res.beta.2")], id = c('N'))
# r2.stacked <- melt(r2[,c("N","res.theta","res.beta")], id = c('N'))
# boxplot(value~N+variable, data=r1.stacked, notch=TRUE, col=rep(1:3, each=8), names=rep(c(50,60,70,80,90,100,500,1000), 3), xlab="N", ylab="Residual", ylim=c(-3, 3))
# abline(h=0)
# legend(x=20,y=3, legend=c("theta=2","beta.Z1=log(2)", "beta.Z1=log(3)"), fill=1:3)

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