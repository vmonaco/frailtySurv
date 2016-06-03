summary.simfrail <- function(object, ...) {
  sim <- object
  # Function to gather the empirical SD, mean estimated SD, coverage, etc.
  fn <- function(name) {
    
    value <- sim[[name]]
    hat <- sim[[paste("hat.", name, sep="")]]
    se <- sim[[paste("se.", name, sep="")]]
    
    cov.95CI <- sum(hat - 1.96*se <= value & value <= hat + 1.96*se)/length(hat)
    
    c(value=unique(value),
      mean.hat=mean(hat),
      sd.hat=sd(hat),
      mean.se=mean(se),
      cov.95CI=cov.95CI
    )
  }
  
  sum.runtime <- sum(sim$runtime)
  mean.runtime <- mean(sim$runtime)
  sd.runtime <- sd(sim$runtime)
  
  param.names <- names(sim)[grepl("^beta|^theta|^Lambda", names(sim))]
  sum.sim <- vapply(param.names, fn, rep(0, 5))
  class(sum.sim) <- append("summary.simfrail", class(sum.sim))
  attributes(sum.sim) <- append(attributes(sum.sim), list(
    reps=attr(sim, "reps"),
    N=unique(sim$N),
    mean.K=mean(sim$mean.K),
    frailty=attr(sim, "frailty"),
    description=paste("Simulation: ", 
                    attr(sim, "reps"), " reps, ",
                    toString(unique(sim$N)), " clusters (avg. size ", format(mean(sim$mean.K), digits=4), "), ",
                    toString(attr(sim, "frailty")), " frailty", "\n",
                    sprintf("Serial runtime (s): %.2f (%.2f +/- %.2f per rep)", sum.runtime, mean.runtime, sd.runtime),
                    sep="")
  ))
  
  sum.sim
}