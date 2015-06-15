# Code for GSoC 2015 proposal: Multivariate Survival Analysis
# Author: Vinnie Monaco

library(survival)
library(frailtypack)
library(cranlogs)

SEED <- 2015 # Used as a default param to the simulation fn below
N <- 300 # Number of families
K <- 2   # Family size
theta <- 2 # Strength of family clustering
reps <- 500 # Number of repetitions for each simulation

# Root dir for output
args <- commandArgs(trailingOnly = TRUE)
PROJECTDIR <- as.character(args[1])

# Generate survival data according to Gorfine et. al, Biometrika 2006
survivalgen <- function(beta = c(log(2)), # Covariate coefficients
                        frailty = c("gamma", "log.normal"), # Frailty distribution
                        censor.mu = 130, # Gaussian distribution for censoring
                        censor.sigma = 15,
                        theta = 2, # Frailty distribution parameter
                        N = 300, K = 2, # Number of families and family size
                        tau = 4.6,
                        C = 0.01,
                        covariates = c("normal", "uniform")
  )
{
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=N*K, ncol=p)
  for (j in 1:p) {
    if (covariates=="normal") {
      Z[, j] <- rnorm(N*K)
    } else if (covariates=="uniform") {
      Z[, j] <- runif(N*K)
    }
  }
  betaZ <- Z%*%beta
  
  # Generate the frailty for each family
  if (frailty=="gamma") {
    w <- 1
    if ( theta != 0 ) {
      w <- rep(rgamma(N, 1/theta, 1/theta), rep(K, N))
    }
    rate <- exp(betaZ)*w 
  } else if (frailty=="log.normal") {
    w <- 1
    if ( theta != 0 ) {
      w <- rep(rnorm(N, 0, sqrt(theta)), rep(K, N))
    }
    rate <- exp(betaZ + w)
  }
  
  # Survival time is similar to that of a late onset disease with default params
  u <- runif(N*K)
  obs.time <- ((-log(u)/rate)^(1/tau))/C
  
  # Non-informative right censoring
  censor.time <- rnorm(N*K, censor.mu, censor.sigma)
  obs.status <- sign(obs.time <= censor.time)
  obs.time <- pmin(obs.time, censor.time)
  
  colnames(Z) <- paste("Z", 1:p, sep="")
  dat <- data.frame(family=rep(1:N, rep(K,N)),
                    rep=rep(1:K, N), 
                    time=obs.time, 
                    status=obs.status,
                    Z)
  return(dat)
}

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
survivalsim <- function(reps, survivalgen.args, coxph.args, output=NULL, outdir=file.path(PROJECTDIR, "results/"), seed=SEED) 
{ 
  # Seed the beginning of each simulation so that each run is independent
  set.seed(seed)
  
  results <- replicate(reps, (function(){
    # Generate new random data, give this to coxph
    data <- do.call("survivalgen", as.list(survivalgen.args))
    coxph.args[["data"]] = data
    coxph.args[["method"]] = "breslow"
    
    # Only interested in total elapsed time
    runtime <- system.time(m <- do.call("coxph", as.list(coxph.args)))["elapsed"]
    names(runtime) <- "runtime"
    
    # Output is: beta1, beta2, ..., [theta], censoring, runtime
    return(c(m$coefficients,
             theta=m$history[[1]]$theta,
             censoring=1-sum(data$status)/length(data$status),
             runtime
             ))
  })())
  
  # Summarize everything
  results <- apply(results, 1, function(x) c(mean=mean(x), sd=sd(x)))
  
  if (is.character(output)) {
    write.csv(results, file.path(outdir, output))
  }
  
  return(results)
}

plot_weekly_package_downloads <- function(packages, dataset, filename) {   
  require("ggplot2")
  require("plyr")
  
  agg1 <- ddply(dataset[dataset$"package" %in% packages,],
                .(time=strftime(as.POSIXlt(date),format="%Y-%W"), package), 
                function(xx) {c(count = sum(xx$count))})
  
  o <- ggplot(agg1, aes(x=time, y=count, color=package, group=package)) +
        geom_line() + ylab("Downloads") + theme_bw() + 
        scale_x_discrete("Year/Week",breaks=agg1$time[seq(1, nrow(agg1), length=10)]) +
        scale_y_discrete("Downloads",breaks=seq(0, round(max(agg1$count), -2) + 1000, 1000)) +
        theme(axis.text.x  = element_text(angle=90, size=12, vjust=0.5),
              panel.grid.major = element_line(colour="gray", size=0.5),
              aspect.ratio = 1/2)
  ggsave(filename)
  return(invisible(agg1))
}

plot_gamma_vs_lognormal <- function(theta, filename=NULL) {   
  require("ggplot2")
  require("plyr")
  
  x <- seq(0, 10, length=100)
  
  labels <- c("Gamma", "Log-normal")
  
  if (is.character(filename)){
    png(filename=filename)
  }
  
  plot(x, dgamma(x, 1/theta, 1/theta), type="l", lty=1, xlab="x value", ylab="Density")
  lines(x, dlnorm(x, 0, sqrt(theta)), lty=2)
  
  legend("topright", inset=.05, title="Distributions", labels, lwd=2, lty=c(1,2))
  
  if (is.character(filename)) {
    dev.off()
  }
  
  return()
}

########## Package statistics

# Generate a figure for weekly downloads
packages <- c("survival","frailtypack","coxme","phmm","survBayes","MST","parfm")
logdata <- cran_downloads(package=packages, from="2012-10-01", to=as.Date(Sys.time()))
plot_weekly_package_downloads(packages, logdata, file.path(PROJECTDIR, "figures/packages_by-week.png"))

# Gamma vs log normal for the default theta
plot_gamma_vs_lognormal(theta, file.path(PROJECTDIR, "figures/gamma_vs_lognormal.png"))

########## Simulations
# Each simulation is run first without and then with frailty estimation

###### Gamma frailty
# beta = log 2, ~35% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2.censor_130.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log2.censor_130.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log2.censor_130.frailty_gamma.frailty_log.normal.csv")

# beta = log 2, ~85% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2.censor_60.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log2.censor_60.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log2.censor_60.frailty_gamma.frailty_log.normal.csv")

# beta = log 3, ~30% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log3.censor_130.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log3.censor_130.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log3.censor_130.frailty_gamma.frailty_log.normal.csv")

# beta = log 3, ~80% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log3.censor_60.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log3.censor_60.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log3.censor_60.frailty_gamma.frailty_log.normal.csv")

# beta = log 2, log 3, ~30% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2log3.censor_130.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gamma(family)),
            "beta_log2log3.censor_130.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gaussian(family)),
            "beta_log2log3.censor_130.frailty_gamma.frailty_log.normal.csv")

# beta = log 2, log 3, ~80% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2log3.censor_60.frailty_gamma.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gamma(family)),
            "beta_log2log3.censor_60.frailty_gamma.frailty_gamma.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gaussian(family)),
            "beta_log2log3.censor_60.frailty_gamma.frailty_log.normal.csv")

###### Log-normal frailty
# beta = log 2, ~35% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2.censor_130.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log2.censor_130.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log2.censor_130.frailty_log.normal.frailty_gamma.csv")

# beta = log 2, ~85% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2.censor_60.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log2.censor_60.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(2)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log2.censor_60.frailty_log.normal.frailty_gamma.csv")

# beta = log 3, ~30% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log3.censor_130.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log3.censor_130.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log3.censor_130.frailty_log.normal.frailty_gamma.csv")

# beta = log 3, ~80% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log3.censor_60.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "beta_log3.censor_60.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "beta_log3.censor_60.frailty_log.normal.frailty_gamma.csv")

# beta = log 2, log 3, ~30% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2log3.censor_130.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gaussian(family)),
            "beta_log2log3.censor_130.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gamma(family)),
            "beta_log2log3.censor_130.frailty_log.normal.frailty_gamma.csv")

# beta = log 2, log 3, ~80% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1),
            "beta_log2log3.censor_60.frailty_log.normal.frailty_none.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gaussian(family)),
            "beta_log2log3.censor_60.frailty_log.normal.frailty_log.normal.csv")
survivalsim(reps, alist(beta=c(log(2), log(3)), censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + Z2 + frailty.gamma(family)),
            "beta_log2log3.censor_60.frailty_log.normal.frailty_gamma.csv")

#### Uniform covariates

# beta = log 2, ~35% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2)), covariates="uniform", censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "uniform.beta_log2.censor_130.frailty_gamma.frailty_gamma.csv")

# beta = log 2, ~85% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(2)), covariates="uniform", censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "uniform.beta_log2.censor_60.frailty_gamma.frailty_gamma.csv")

# beta = log 3, ~30% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(3)), covariates="uniform", censor.mu=130, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "uniform.beta_log3.censor_130.frailty_gamma.frailty_gamma.csv")

# beta = log 3, ~80% censoring, gamma frailty
survivalsim(reps, alist(beta=c(log(3)), covariates="uniform", censor.mu=60, frailty="gamma", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gamma(family)),
            "uniform.beta_log3.censor_60.frailty_gamma.frailty_gamma.csv")

# beta = log 2, ~35% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2)), covariates="uniform", censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "uniform.beta_log2.censor_130.frailty_log.normal.frailty_log.normal.csv")

# beta = log 2, ~85% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(2)), covariates="uniform", censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "uniform.beta_log2.censor_60.frailty_log.normal.frailty_log.normal.csv")

# beta = log 3, ~30% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(3)), covariates="uniform", censor.mu=130, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "uniform.beta_log3.censor_130.frailty_log.normal.frailty_log.normal.csv")

# beta = log 3, ~80% censoring, log.normal frailty
survivalsim(reps, alist(beta=c(log(3)), covariates="uniform", censor.mu=60, frailty="log.normal", N=N, K=K, theta=theta), 
            alist(Surv(time, status) ~ Z1 + frailty.gaussian(family)),
            "uniform.beta_log3.censor_60.frailty_log.normal.frailty_log.normal.csv")

