# frailtyr simulations

library(frailtyr)

SEED <- 2015 # Meta-seed for each simulation run
N <- c(50,100,500,1000) # Number of clusters
REPS <- 100 # Number of repetitions for each simulation
BASE.TIME <- 0:200 # where to evaluate the baseline hazard

# defaults when not specified
K <- 2 # default cluster sizes
BETA <- log(2)
THETA <- 2
FRAILTY <- "gamma" 
CENSOR.PARAMS <- c(130,15)

# "Nice" baseline hazard
lambda_0 = function(t, tau=4.6, C=0.01) 
  (tau*(C*t)^tau)/t

Lambda_0 = function(t, tau=4.6, C=0.01) 
  (C*t)^tau

# Increasing oscillating baseline hazard
lambda_0.osc = function(t, tau=4.6, C=0.01, A=2, f=0.1) 
  A^sin(f*pi*t) * (tau*(C*t)^tau)/t

Lambda_0.osc = function(t, tau=4.6, C=0.01, A=2, f=0.1) 
  integrate(lambda_0.osc, 0, t, subdivisions = 1000L)$value

# dir for output
args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- as.character(args[1])

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrail.multigen.output <- function(name, ...) { 
  results <- simfrail.multigen(...)
  
  if (is.character(output)) {
    write.csv(results, file.path(outdir, paste(output, ".csv", sep="")))
  }
  
  return(results)
}


########## Simulations
# Each simulation is run first without and then with frailty estimation

###### Gamma frailty
# beta = log 2, ~35% censoring, gamma frailty
results <- simfrail.multigen(10, 2015, alist(censor.params=c(130,15), frailty="gamma", N=50, K=2, theta=2, covariates="uniform"), alist(formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), frailty="gamma", verbose=FALSE), c(1,20,50,100), "beta", list(c(2,log(2)),c(log(2),log(3))))