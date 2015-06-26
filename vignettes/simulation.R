# frailtyr simulations

library(frailtyr)

SEED <- 2015 # Meta-seed for each simulation run
N <- c(50,100,200,400)#,500,1000) # Number of clusters
REPS <- 100 # Number of repetitions for each simulation
BASE.TIME <- 0:150 # where to evaluate the baseline hazard
BASE.TIME.RES <- c(25,50,75,100) # where to evaluate residulas of the baseline hazard

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

Lambda_0.osc = Vectorize(function(t, tau=4.6, C=0.01, A=2, f=0.1) 
  ifelse(t > 0, integrate(lambda_0.osc, 0, t, subdivisions = 1000L)$value, 0))

# dir for output
args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- as.character(args[1])
# OUTDIR <- "/home/vinnie/tmp/frailtyr"

# A wrapper for coxph, gathers the coeffs, theta, censoring, runtime
simfrail.multigen.output <- function(name, title=NULL, outdir=OUTDIR, ...) {
  dots <- list(...)
  
  results <- simfrail.multigen(...)
  
  write.csv(results, file.path(outdir, paste(name, ".csv", sep="")))
  
  # Gather true values
  beta <- eval(dots$genfrail.args$beta)
  names(beta) <- paste("beta.", 1:length(beta), sep="")
  
  theta <- eval(dots$genfrail.args$theta)
  names(theta) <- paste("theta.", 1:length(theta), sep="")
  
  Lambda_0.points <- eval(dots$genfrail.args$Lambda_0)(BASE.TIME.RES)
  names(Lambda_0.points) <- paste("Lambda.", BASE.TIME.RES, sep="")
  
  true.values <- c(beta, theta, Lambda_0.points)
  
  # Plot residuals
  pdf(file.path(outdir, paste(name, "_residuals.pdf", sep="")), width=18, height=6)
  plot.simfrail.residuals(results, true.values=true.values, title=title)
  dev.off()
  
  lambdafn <- eval(dots$genfrail.args$lambda_0)
  Lambdafn <- eval(dots$genfrail.args$Lambda_0)
  
  # Plot the BH for each N
  pdf(file.path(outdir, paste(name, "_hazard.pdf", sep="")), width=12, height=6)
  lapply(unique(as.vector(t(results["N"]))), function(Ni) {
    plot.simfrail.hazard(results[results["N"]==Ni,], 
                         true.cbh=Lambdafn, true.bh=lambdafn, 
                         title=paste(title, "N = ", Ni))
  })
  dev.off()
  
  return(NULL)
}

########## Simulations
# Completed simulations are commented out

# simfrail.multigen.output("gamma-frailty_2-covars", "Gamma frailty, 2 covariates",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="gamma", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="gamma", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("gamma-frailty_poisson-clusters", "Gamma frailty, K~TP(2)",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="gamma", 
#                            censor.params=CENSOR.PARAMS,
#                            K="poisson", theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="gamma", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("gamma-frailty_oscillating-bh", "Gamma frailty, oscillating BH",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="gamma", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0.osc,
#                            Lambda_0=Lambda_0.osc), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="gamma", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("gamma-frailty_rounded10-obs-time", "Gamma frailty, obs time rounded to 10",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="gamma", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0,
#                            round=10), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="gamma", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("gamma-frailty_rounded1-obs-time", "Gamma frailty, obs time rounded to 10",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="gamma", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0,
#                            round=1), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="gamma", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("lognormal-frailty_2-covars", "Lognormal frailty, 2 covariates",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="lognormal", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="lognormal", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("pvf-frailty_2-covars", "PVF frailty, 2 covariates",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="pvf", 
#                            censor.params=c(110,15),
#                            K=K, theta=0.7, 
#                            covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="pvf", verbose=FALSE), 
#                          base.time=0:130, 
#                          param.name="N", param.values=N)
# 
# simfrail.multigen.output("invgauss-frailty_2-covars", "Inv gauss frailty, 2 covariates",
#                          reps=REPS, seed=SEED, 
#                          genfrail.args=alist(
#                            beta=c(log(2),log(3)),
#                            frailty="invgauss", 
#                            censor.params=CENSOR.PARAMS,
#                            K=K, theta=THETA, covariates="uniform",
#                            lambda_0=lambda_0,
#                            Lambda_0=Lambda_0), 
#                          fitfrail.args=alist(
#                            formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
#                            frailty="invgauss", verbose=FALSE), 
#                          base.time=BASE.TIME, 
#                          param.name="N", param.values=N)

simfrail.multigen.output("lognormal-frailty_2-covars_cens30", "Lognormal frailty, 2 covariates",
                         reps=REPS, seed=SEED, 
                         genfrail.args=alist(
                           beta=c(log(2),log(3)),
                           frailty="lognormal", 
                           censor.params=c(110,15),
                           K=K, theta=THETA, covariates="uniform",
                           lambda_0=lambda_0,
                           Lambda_0=Lambda_0), 
                         fitfrail.args=alist(
                           formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                           frailty="lognormal", verbose=FALSE), 
                         base.time=0:130, 
                         param.name="N", param.values=N)