
sum_by_cluster <- function(A_) {
  A_dot <- sapply(names(A_), function(i) {
    rep(0, length(A_[[i]][1,]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  for (i in names(A_)) {
    A_dot[[i]] <- colSums(A_[[i]])
  }
  
  A_dot
}

sharedfrailty.fit <- function(x, y, cluster, beta_init, theta_init, dfrailty, 
                              dfrailty_dtheta, control, rownames) {
  
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
  fail_idx <- 1 + rank(time)
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
  
  # Create psi_i for each cluster. Depends on theta, N_i.(t) and H_i.(t)
  psi_ <- sapply(cluster_names, function (i) {
    function(theta, N_dot_ik, H_dot_ik) {
      phi_1 <- function (w) {
        ( w ^ N_dot_ik ) * exp(-w*H_dot_ik) * dfrailty(w, theta)
      }
      phi_2 <- function (w) {
        ( w ^ (N_dot_ik + 1) ) * exp(-w*H_dot_ik) * dfrailty(w, theta)
      }
      
#       integrate(Vectorize(phi_2), 0, Inf, stop.on.error = FALSE)$value/
#         integrate(Vectorize(phi_1), 0, Inf, stop.on.error = FALSE)$value
      # Simplified for gamma, 
      (N_dot_ik + 1/theta)/(H_dot_ik + 1/theta)
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  H_ <- sapply(cluster_names, function(i) {
    matrix(0, cluster_sizes[[i]], length(time_steps))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  H_dot <- sapply(cluster_names, function(i) {
    rep(0, length(time_steps))
  }, simplify = FALSE, USE.NAMES = TRUE)

  # Step function for baseline hazard
  estimate_lambda_hat <- function(beta, theta) {
    # delta_lambda_hat[1] = 0, since at t_0 there is no chance of failure
    delta_lambda_hat = rep(0, length(time_steps))
    lambda_hat = rep(0, length(time_steps))
    
    for (k in 2:length(time_steps)) {
      # Denom for delta_lambda_hat[k]
      denom <- sum(sapply(cluster_names, function(i) {
        psi_[[i]](theta, N_dot[[i]][k-1], H_dot[[i]][k-1]) * 
          sum(vapply(1:cluster_sizes[[i]], function(j) {
            Y_[[i]][j, k] * exp(t(beta) %*% spx[[i]][j,])
          }, 0))
      }))
      
      delta_lambda_hat[k] <- d_[k]/denom
      
      # Update the cumsum since it's needed several times in the next iteration
      lambda_hat[k] <- lambda_hat[k-1] + delta_lambda_hat[k]
      
      # Determine H_ij(t) and H_i.(t)
      for (i in cluster_names) {
        for (j in 1:cluster_sizes[[i]]) {
          k_min = min(spk[[i]][j], k)
          H_[[i]][j, k] <- lambda_hat[k_min] * exp(t(beta) %*% spx[[i]][j,])
        }
        H_dot[[i]][k] <- sum(H_[[i]][,k])
      }
    }
    
    list(lambda_hat=lambda_hat, H_=H_, H_dot=H_dot)
  }
  
  U_r <- function(H_, H_dot, beta_idx, theta) {
    term1 <- mean(sapply(cluster_names, function(i) {
      sum(sapply(1:cluster_sizes[[i]], function(j) {
          spstatus[[i]][j] * spx[[i]][j,beta_idx]
          }))
    }))
    
    term2 <- mean(sapply(cluster_names, function(i) {
      sum(sapply(1:cluster_sizes[[i]], function(j) {
        H_[[i]][j, spk[[i]][j]] * spx[[i]][j,beta_idx]})) *
          psi_[[i]](theta, N_dot[[i]][tau_k], H_dot[[i]][tau_k])
    }))
    term1 - term2
  }
  
  U_p <- function(H_dot, theta_idx, theta) {
    mean(sapply(cluster_names, function(i) {
      numer = integrate(function (w) {
        ( w ^ N_dot[[i]][tau_k] ) * exp(-w*H_dot[[i]][tau_k]) * Vectorize(function(wi) dfrailty_dtheta(wi, theta, theta_idx))(w)
      }, 0, Inf, stop.on.error = FALSE)$value
      
      denom = integrate(function (w) {
        ( w ^ N_dot[[i]][tau_k] ) * exp(-w*H_dot[[i]][tau_k]) * dfrailty(w, theta)
      }, 0, Inf, stop.on.error = FALSE)$value
      
      numer/denom
    }))
  }

  # Build the system of equations where gamma is c(beta, theta)
  U <- function(gamma) {
    beta <- gamma[1:n_beta]
    theta <- gamma[(n_beta+1):(n_beta+n_theta)]
    cat('beta',beta,'theta',theta,'\n')
    
    estimator <- estimate_lambda_hat(beta, theta)
    H_ <- estimator$H_
    H_dot <- estimator$H_dot
    
    c(sapply(1:n_beta, function(beta_idx) U_r(H_, H_dot, beta_idx, theta)),
      sapply(1:n_theta, function(theta_idx) U_p(H_dot, theta_idx, theta)))
  }
  
  fit <- nleqslv(c(beta_init, theta_init), U)
  gamma_hat <- fit$x
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  
  fit
}