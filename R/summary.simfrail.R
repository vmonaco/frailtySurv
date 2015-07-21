summary.simfrail <- function(sim) {
  
  # Compare the empirical SD to the mean estimated SD for each parameter
  true.beta <- sim$beta
  true.theta <- sim$theta
  true.Lambda <- sim$Lambda
  
  est.beta <- do.call(rbind, sim[paste("beta.", 1:length(sim$beta), sep="")])
  est.theta <- do.call(rbind, sim[paste("theta.", 1:length(sim$theta), sep="")])
  est.Lambda <- do.call(rbind, sim[paste("Lambda.", 1:length(sim$Lambda), sep="")])
  
  sd.beta <- do.call(rbind, sim[paste("sd.beta.", 1:length(sim$beta), sep="")])
  sd.theta <- do.call(rbind, sim[paste("sd.theta.", 1:length(sim$theta), sep="")])
  sd.Lambda <- do.call(rbind, sim[paste("sd.Lambda.", 1:length(sim$Lambda), sep="")])
  
  emp.sd.beta <- apply(est.beta, 1, sd)
  emp.sd.theta <- apply(est.theta, 1, sd)
  emp.sd.Lambda <- apply(est.Lambda, 1, sd)
  
  est.sd.beta <- apply(sd.beta, 1, mean)
  est.sd.theta <- apply(sd.theta, 1, mean)
  est.sd.Lambda <- apply(sd.Lambda, 1, mean)
  
  cov.beta <- rowSums(est.beta - 2*sd.beta <= true.beta 
                           & true.beta <= est.beta + 2*sd.beta)/sim$reps
  cov.theta <- rowSums(est.theta - 2*sd.theta <= true.theta 
                           & true.theta <= est.theta + 2*sd.theta)/sim$reps
  cov.Lambda <- rowSums(est.Lambda - 2*sd.Lambda <= true.Lambda 
                           & true.Lambda <= est.Lambda + 2*sd.Lambda)/sim$reps
  
  summary <- data.frame(cov.beta=cov.beta)
  
  summary
}