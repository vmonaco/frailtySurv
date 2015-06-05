
sum_by_cluster <- function(A_) {
  A_dot <- sapply(names(A_), function(i) {
    rep(0, length(A_[[i]][1,]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  for (i in names(A_)) {
    A_dot[[i]] <- as.double(colSums(A_[[i]]))
  }
  
  A_dot
}

# Helper fn to get the rank, where duplicates have the same rank
increasing_rank <- function(d) {
  j <- unique(sort(d));
  return(sapply(d,function(dd) which(dd==j)));
}

fitfrail.fit <- function(x, y, cluster, beta_init, theta_init, frailty, control, rownames) {
  
  # Num observations
  n <- nrow(y)
  
  # Num frailty distribution params
  n_theta <- length(theta_init)
  
  # Num coefs
  if (is.matrix(x)) {
    n_beta <- ncol(x)
  } else {
    if (length(x)==0) n_beta <-0
    else n_beta <- 1
  }
  
  time <- y[,1]
  status <- y[,2]
  
  # Initialize theta and beta to zero vectors if initial params are missing
  if (!missing(theta_init) && length(theta_init) > 0) {
    if (length(theta_init) != n_theta) stop("Wrong length for frality distribution inital values")
  } else {
    theta_init <- rep(0,n_theta)
  }
  
  if (!missing(beta_init) && length(beta_init) > 0) {
    if (length(beta_init) != n_beta) stop("Wrong length for coefficient inital values")
  } else {
    beta_init <- rep(0,n_beta)
  }
  
  # Unique time steps, including t=0. This is where the baseline hazard actually increases
  time_steps <- c(0, sort(unique(time)))
  # Index of the end of the observation period
  tau_k <- length(time_steps)
  
  # The _ variables are all lists, indexed by the cluster name
  X_ <- split.data.frame(x, cluster) # Covariates
  T_ <- split(time, cluster) # Followup time
  I_ <- split(status, cluster) # Failure indicator
  # Failure rank + 1, where time_steps[R_[i][j]] is the failure time of member j in cluster i
  R_ <- split(1 + increasing_rank(time), cluster) 
  
  # Cluster names and sizes are needed in several places, so get these now
  cluster_names <- levels(cluster)
  cluster_sizes <- table(cluster)
  
  # N_[[i]][j, k] indicates whether individual j in cluster i failed at time <= tau_k
  N_ <- sapply(cluster_names, function(i) {
    t(sapply(1:cluster_sizes[[i]], function(j) {
      as.numeric((I_[[i]][j] > 0) & (T_[[i]][j] <= time_steps))
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # N_dot[[i]][k] is the number of individuals in cluster i that failed up to time tau_k
  N_dot <- sum_by_cluster(N_)
  
  # Y_[[i]][j,k] indicates whether individual j in cluster i failed at time >= tau_k
  Y_ <- sapply(cluster_names, function(i) {
    t(sapply(1:cluster_sizes[[i]], function(j) {
      as.numeric(T_[[i]][j] >= time_steps)
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # d_[k] is the number of failures at time tau_k
  d_ <- vapply(seq_along(time_steps), function(k) {
    sum(status[time == time_steps[k]])
  }, 0)
  
  verbose = control$verbose
  if (verbose)
    cat(c("Iter", paste("beta_", 1:n_beta, sep=""), paste("theta_", 1:n_beta, sep=""), "\n"), sep="\t")
  
  iter <- 0
  # Build the system of equations where gamma is c(beta, theta)
  U <- function(gamma) {
    beta <- gamma[1:n_beta]
    theta <- gamma[(n_beta+1):(n_beta+n_theta)]
    
    estimator <- baseline_hazard_estimator(X_, R_, d_, Y_, N_dot, beta, theta, frailty)
    H_ <<- estimator$H_
    H_dot <<- estimator$H_dot
    lambda_hat <<- estimator$lambda_hat
    
    if (verbose) {
      iter <<- iter + 1
      cat(c(signif(c(iter, beta, theta), 3), "\n"), sep="\t")
    }
    
    # TODO: This could be done in parallel, not much gain though.
    c(sapply(1:n_beta,
             function(beta_idx)
               U_r(X_, R_, I_, N_dot, H_, H_dot, beta, theta, beta_idx, frailty)),
      sapply(1:n_theta, 
             function(theta_idx) 
               U_p(X_, R_, I_, N_dot, H_, H_dot, beta, theta, theta_idx, frailty)))
  }
  
  fit <- nleqslv(c(beta_init, theta_init), U, control=list(maxit=control$iter.max))
  
  if (verbose)
    cat("Converged after", iter, "iterations\n")
  
  gamma_hat <- fit$x
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  loglik <- log_likelihood(X_, R_, I_, N_dot, H_dot, lambda_hat, beta_hat, theta_hat, frailty)
  
  list(beta = beta_hat,
       theta = theta_hat,
       lambda = lambda_hat,
       loglik = loglik,
       method='fitfrail',
       frailty=frailty,
       loglikfn = function(beta, theta, frailty) {
         estimator <- baseline_hazard_estimator(X_, R_, d_, Y_, N_dot, beta, theta, frailty)
         H_ <<- estimator$H_
         H_dot <<- estimator$H_dot
         lambda_hat <<- estimator$lambda_hat
         log_likelihood(X_, R_, I_, N_dot, H_dot, lambda_hat, beta, theta, frailty)}
      )
}