
sum_by_cluster <- function(A_) {
  A_dot <- sapply(names(A_), function(i) {
    rep(0, length(A_[[i]][1,]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  for (i in names(A_)) {
    A_dot[[i]] <- as.double(colSums(A_[[i]]))
  }
  
  A_dot
}

Rank<-function(d) {
  j<-unique(sort(d));
  return(sapply(d,function(dd) which(dd==j)));
}

sharedfrailty.fit <- function(x, y, cluster, beta_init, theta_init, frailty, control, rownames) {
  
  n_theta <- 1
  
  n <- nrow(y)
  if (is.matrix(x)) {
    n_beta <- ncol(x)
  } else {
    if (length(x)==0) n_beta <-0
    else n_beta <-1
  }
  time <- y[,1]
  status <- y[,2]

  maxiter <- control$iter.max
  
  if (!missing(beta_init) && length(beta_init)>0) {
    if (length(beta_init) != n_beta) stop("Wrong length for coefficient inital values")
  } else {
    beta_init <- rep(0,n_beta)
  }

  if (!missing(theta_init) && length(theta_init)>0) {
    if (length(theta_init) != n_theta) stop("Wrong length for frality distribution inital values")
  } else {
    theta_init <- rep(0,n_theta)
  }
  
  # The order of failure (or time index) for the data, where fail_idx==1 is t==0
  fail_idx <- 1 + Rank(time)
  # Unique time steps. This is where the baseline hazard actually increases
  time_steps <- c(0, sort(unique(time)))
  # End of the observation period
  tau_k <- length(time_steps)
  
  # Split the time, status into their respective clusters
  sptime <- split(time, cluster)
  spstatus <- split(status, cluster)
  spk <- split(fail_idx, cluster)
  spx <- split.data.frame(x, cluster)
  
  # Cluster names and sizes are needed in several places, so get these now
  cluster_names <- levels(cluster)
  cluster_sizes <- table(cluster)
  
  # N_[[i]][j, k] indicates whether individual j in cluster i failed at time <= tau_k
  N_ <- sapply(cluster_names, function(i) {
    t(sapply(1:cluster_sizes[[i]], function(j) {
      as.numeric((spstatus[[i]][j] > 0) & (sptime[[i]][j] <= time_steps))
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # N_dot[[i]][k] is the number of individuals in cluster i that failed up to time tau_k
  N_dot <- sum_by_cluster(N_)
  
  # Y_[[i]][j,k] indicates whether individual j in cluster i failed at time >= tau_k
  Y_ <- sapply(cluster_names, function(i) {
    t(sapply(1:cluster_sizes[[i]], function(j) {
      as.numeric(sptime[[i]][j] >= time_steps)
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # d_[k] is the number of failures at time tau_k
  d_ <- vapply(seq_along(time_steps), function(k) {
    sum(status[time == time_steps[k]])
  }, 0)
  
  iter <- 1
  # Build the system of equations where gamma is c(beta, theta)
  U <- function(gamma) {
    beta <- gamma[1:n_beta]
    theta <- gamma[(n_beta+1):(n_beta+n_theta)]
    cat("Iteration ", iter, " : ", "beta <- ", beta, ", theta <- ", theta, "\n")
    
    estimator <- baseline_hazard_estimator(spx, spk, d_, Y_, N_dot, beta, theta, frailty)
    H_ <<- estimator$H_
    H_dot <<- estimator$H_dot
    lambda_hat <<- estimator$lambda_hat
    
    iter <<- iter + 1
    
    # TODO: This could be done in parallel, not much gain though.
    c(sapply(1:n_beta,
             function(beta_idx)
               U_r(spk, N_dot, spstatus, spx, H_, H_dot, beta, theta, beta_idx, frailty)),
      sapply(1:n_theta, 
             function(theta_idx) 
               U_p(spk, N_dot, spstatus, spx, H_, H_dot, beta, theta, theta_idx, tau_k, frailty)))
  }
  
  fit <- nleqslv(c(beta_init, theta_init), U)
  gamma_hat <- fit$x
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  loglik <- log_likelihood(beta_hat, theta_hat, lambda_hat, H_dot, spstatus, spk, spx, N_dot, frailty)
  
  list(beta = beta_hat,
       theta = theta_hat,
       loglik = loglik,
       method='fitfrail'
       
       )
}