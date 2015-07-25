fitfrail.fit <- function(x, y, cluster, beta_init, theta_init, frailty, 
                         control, rownames, weights) {
  
  # TODO: error check for number of frailty distr params
  if (missing(theta_init)) # || length(theta_init) != n_frailty_params(frailty))
    stop("Wrong length for frailty distribution inital values")
    
  if (!missing(beta_init) && length(beta_init) != ncol(x))
    stop("Wrong length for coefficient inital values")
  
  # gamma = c(beta, theta)
  gamma_init <- c(beta_init, theta_init)
  
  # Num observations
  n <- nrow(y)
  
  # Num coefs
  n_beta <- length(beta_init)
  
  # Num frailty distribution params
  n_theta <- length(theta_init)
  
  # Total num params
  n_gamma <- n_beta + n_theta
  
  # time and failure indicator
  time <- y[,1]
  status <- y[,2]
  
  # Sorted time
  Lambda.time <- c(0, sort(unique(time[status > 0])))
  time_sorted_idx <- order(time)
  time_steps <- c(0, time[time_sorted_idx])
  k_tau <- length(time_steps) # index of the last time step
  
  # d_[k] is the number of failures at time t_k
#   d_ <- vapply(seq_along(time_steps), function(k) {
#     sum(status[time == time_steps[k]])
#   }, 0)
  d_ <- c(0, unname(status[time_sorted_idx]))
  
  # Variables indexed by [[i]][j,] where i is cluster name and j is member
  X_ <- split.data.frame(x, cluster) # Covariates
  T_ <- split(time, cluster) # Followup time
  I_ <- split(status, cluster) # Failure indicator
  
  # Helper fn to get the rank, where duplicates have the same rank
  increasing_rank <- function(d) {
    j <- unique(sort(d));
    return(sapply(d,function(dd) which(dd==j)));
  }
  
  # Index of the observed time, K_[[i]][j]
  # K_ <- split(1 + increasing_rank(time), cluster) 
  K_ <- split(1 + rank(time, ties.method="first"), cluster)
  
  # Cluster names and sizes are needed in several places, so get these now
  cluster_names <- levels(cluster)
  names(cluster_names) <- cluster_names
  n_clusters <- length(cluster_names)

  cluster_sizes <- table(cluster)
  
  if (is.null(weights)) {
    weights <- rep(1, n_clusters)
  } else {
    stopifnot(length(weights) == n_clusters)
  }
  
  # Convenience function to apply a function over the cluster and members indices
  clustapply <- function(fn, FUN.VALUE) {
    lapply(cluster_names, function(i) {
      t(vapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      }, FUN.VALUE=FUN.VALUE))})
  }
  
  clustlapply <- function(fn) {
    lapply(cluster_names, function(i) {
      lapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      })})
  }
  
  # N_[[i]][j,k] == 1 if ij failed on/before time t_k
  N_ <- clustapply(function(i, j) {
    as.numeric((I_[[i]][j] > 0) & (K_[[i]][j] <= 1:k_tau))
    }, rep(0, k_tau))
  
  # N_dot_[[i]][k] is the number of individuals in cluster i that failed at time <= t_k
  N_dot_ <- lapply(N_, colSums)
  
  # N_tilde_[[i]][j,k] == 1 if ij was observed at time <= t_k
  N_tilde_ <- clustapply(function(i, j) {
    as.numeric(K_[[i]][j] <= 1:k_tau)
    }, rep(0, k_tau))
  
  # Y_[[i]][j,k] == 1 if ij failed at time >= t_k
  Y_ <- clustapply(function(i, j) {
      as.numeric(K_[[i]][j] >= 1:k_tau)
    }, rep(0, k_tau))
  
  # VARS holds all the variables that depend on \hat{beta} and \hat{theta}
  VARS <- new.env()
  VARS$iter <- 1
  VARS$N_ <- N_
  VARS$N_dot_ <- N_dot_
  VARS$N_tilde_ <- N_tilde_
  VARS$Y_ <- Y_
  VARS$K_ <- K_
  VARS$time_steps <- time_steps
  VARS$d_ <- d_
  
  verbose <- control$verbose
  fitmethod <- control$fitmethod
  # The objective function for either score equations or loglikelihood optimization
  # gamma is c(beta, theta), several global variables are updated since these are
  # reused in other places (jacobian, covar matrix)
  fit_fn <- function(gamma) {
    
    VARS$beta <- gamma[1:n_beta]
    VARS$theta <- gamma[(n_beta+1):(n_beta+n_theta)]
    
    VARS$R_ <- clustapply(function(i, j) {
      exp(sum(VARS$beta * X_[[i]][j,])) * Y_[[i]][j,]
    }, rep(0, k_tau))
    VARS$R_dot_ <- lapply(VARS$R_, colSums)
    
    bh_ <- bh(d_, X_, K_, Y_, N_, N_dot_, VARS$beta, VARS$theta, frailty, weights)
    
    VARS$H_ <- bh_$H_
    VARS$H_dot_ <- bh_$H_dot
    VARS$lambda <- bh_$lambda
    VARS$Lambda <- bh_$Lambda
    VARS$psi_ <- bh_$psi_
    VARS$phi_1_ <- bh_$phi_1_
    VARS$phi_2_ <- bh_$phi_2_
    
    VARS$phi_prime_1_ <- lapply(1:n_theta, function(theta_idx) 
      phi_prime_k(1, theta_idx, N_dot_, VARS$H_dot_, VARS$theta, frailty))
    
    VARS$loglik_vec <- loglikelihood(X_, K_, I_, VARS$phi_1_, VARS$lambda, VARS$beta)
    VARS$loglik <- sum(weights * VARS$loglik_vec)
    
    if (verbose) {
      if (VARS$iter == 1) {
        cat(c("Iter", 
              paste("beta.", 1:n_beta, sep=""), 
              paste("theta.", 1:n_theta, sep=""), 
              "loglik",
              "\n"), sep="\t")
      }
      cat(c(VARS$iter, format(round(c(VARS$beta, VARS$theta, VARS$loglik), 
                         4), nsmall=4, trim=TRUE), "\n"), sep="\t")
    }
    # TODO: This could be done in parallel, not much gain though.
    xi_beta_ <- matrix(vapply(1:n_beta, function(r) {
      xi_beta(X_, I_, VARS$H_, VARS$psi_, r)
    }, rep(0, n_clusters)), nrow=n_clusters, ncol=n_beta)
    
    xi_theta_ <- matrix(vapply(1:n_theta, function(r) {
      xi_theta(VARS$phi_1_, VARS$phi_prime_1_[[r]], r)
    }, rep(0, n_clusters)), nrow=n_clusters, ncol=n_theta)
    
    # xi_[i,r], where i is cluster, r is param idx
    # Update globally, this is used in the jacobian and covariance matrix
    VARS$xi_ <- cbind(xi_beta_, xi_theta_)
    VARS$U_ <- colMeans(weights * VARS$xi_)
    
    VARS$iter <- VARS$iter + 1
    
    if (fitmethod == "loglik") {
      VARS$loglik
    } else if (fitmethod == "score") {
      VARS$U_
    } else if (fitmethod == "tmp") {
      colSums(VARS$xi_)
    } else if (fitmethod == "tmp1") {
      VARS$xi_
    } else if (fitmethod == "tmp2") {
      VARS$loglik_vec
    }
  }
  
  # Already computed
  loglik_jacobian <- function(gamma=NULL) {
    VARS$U_
  }
  
  score_jacobian <- function(gamma=NULL) {
    
    VARS$phi_3_ <- phi_k(3, N_dot_, VARS$H_dot_, VARS$theta, frailty)
    
    VARS$phi_prime_2_ <- lapply(1:n_theta, function(theta_idx) 
      phi_prime_k(2, theta_idx, N_dot_, VARS$H_dot_, VARS$theta, frailty))
    
    VARS$phi_prime_prime_1_ <- lapply(1:n_theta, function(theta_idx_1) {
      lapply(1:n_theta, function(theta_idx_2) {
        phi_prime_prime_k(1, theta_idx_1, theta_idx_2, 
                          N_dot_, VARS$H_dot_, VARS$theta, frailty, k_tau)})})
    
    dH_dbeta_ <- lapply(1:n_beta, function(s) {
      dH_dbeta(s, d_, X_, K_, VARS$R_, VARS$R_dot_, 
               VARS$phi_1_, VARS$phi_2_, VARS$phi_3_, 
               VARS$Lambda, VARS$lambda,
               VARS$beta, VARS$theta, frailty)
    })
    
    dH_dtheta_ <- lapply(1:n_theta, function(s) {
      dH_dtheta(d_, X_, K_, VARS$R_, VARS$R_dot_,
                VARS$phi_1_, VARS$phi_2_, VARS$phi_3_,
                VARS$phi_prime_1_[[s]], VARS$phi_prime_2_[[s]],
                VARS$Lambda, VARS$lambda, VARS$beta)
    })
    
    jacobian_ls <- function(l, s) {
      if (l <= n_beta && s <= n_beta) {
        jacobian_beta_beta(l, 
                           X_, K_, VARS$H_, 
                           VARS$phi_1_, VARS$phi_2_, VARS$phi_3_,
                           dH_dbeta_[[s]]$dH_dbeta_, 
                           dH_dbeta_[[s]]$dH_dot_dbeta_)
      } else if (l <= n_beta && s > n_beta) {
        theta_idx <- s - n_beta
        jacobian_beta_theta(l,
                            X_, K_,  VARS$H_, 
                            VARS$phi_1_, VARS$phi_2_, VARS$phi_3_,
                            VARS$phi_prime_1_[[theta_idx]], VARS$phi_prime_2_[[theta_idx]],
                            dH_dtheta_[[theta_idx]]$dH_dtheta_, 
                            dH_dtheta_[[theta_idx]]$dH_dot_dtheta_)
      } else if (l > n_beta && s <= n_beta) {
        theta_idx <- l - n_beta
        jacobian_theta_beta(VARS$phi_1_, VARS$phi_2_,
                            VARS$phi_prime_1_[[theta_idx]], VARS$phi_prime_2_[[theta_idx]],
                            dH_dbeta_[[s]]$dH_dot_dbeta_)
      } else if (l > n_beta && s > n_beta) {
        theta_idx_1 <- l - n_beta
        theta_idx_2 <- s - n_beta
        jacobian_theta_theta(VARS$phi_1_, VARS$phi_2_,
                             VARS$phi_prime_1_[[theta_idx_1]], VARS$phi_prime_2_[[theta_idx_1]], 
                             VARS$phi_prime_prime_1_[[theta_idx_1]][[theta_idx_2]], 
                             dH_dtheta_[[theta_idx_2]]$dH_dot_dtheta_)
      }
    }
    
    # Build the jacobian matrix from the respective components
    outer(1:n_gamma, 1:n_gamma, Vectorize(function(l, s) 
      jacobian_ls(l, s)))
  }
  
  # The covariance matrix is encapsulated in a function since it takes some effort
  # to compute. By default, it is only computed when accessed by the user
  covariance <- function() {
    # jacobian is computed first, sets up some vars that are needed in steps 1-3
    jac <- score_jacobian()
    
    ################################################################## Step I
    
    V <- Reduce('+', lapply(1:n_clusters, function(i) {
      VARS$xi_[i,] %*% t(VARS$xi_[i,])
      }))/n_clusters
    # V <- (t(xi_) %*% xi_)/n_clusters
    
    ################################################################## Step II
    
    R_star <- clustapply(function(i, j) {
      exp(sum(VARS$beta * X_[[i]][j,]))
    }, 0)
    
    # Split Q up into beta and theta components, then create a Q_ function
    # that accesses component by index r
    Q_beta_ <- lapply(1:n_beta, function(r) {
      Q_beta(X_, K_, VARS$H_, R_star, VARS$phi_1_, VARS$phi_2_, VARS$phi_3_, r)
    })
    
    Q_theta_ <- lapply((n_beta+1):(n_beta+n_theta), function(r) {
      Q_theta(VARS$H_, R_star, VARS$phi_1_, VARS$phi_2_,
              VARS$phi_prime_1_[[r - n_beta]], VARS$phi_prime_2_[[r - n_beta]])
    })
    
    # Q_[[r]] has the same shape as H_dot_
    Q_ <- c(Q_beta_, Q_theta_)
    
    # calligraphic Y, Ycal_[t]
    Ycal_ <- Ycal(X_, Y_, VARS$psi_, VARS$beta)
    
    # eta_[t]
    eta_ <- eta(VARS$phi_1_, VARS$phi_2_, VARS$phi_3_)
    
    # Upsilon
    Upsilon_ <- Upsilon(X_, K_, VARS$R_dot_, eta_, Ycal_, VARS$beta)
    
    # Omega_[[i]][[j, t]
    # Compute everything
    Omega_ <- Omega(X_, N_, VARS$R_dot_, eta_, Ycal_, VARS$beta)
    
    # p_hat_[t]
    p_hat_ <- p_hat(I_, Upsilon_, Omega_, N_tilde_)
    
    # pi_[[r]][s]
    pi_ <- lapply(1:n_gamma, function(r) {
      pi_r(Q_[[r]], N_tilde_, p_hat_)
    })
    
    G <- outer(1:n_gamma, 1:n_gamma, Vectorize(function(r, l) {
      G_rl(pi_[[r]], pi_[[l]], p_hat_, Ycal_, N_)
    }))
    
    ################################################################## Step III
    
    # M_hat_[[i]][j,s]
    M_hat_ <- M_hat(X_, N_, Y_, VARS$psi_, VARS$beta, VARS$Lambda)
    
    # u_star_[i,r]
    u_star_ <- u_star(pi_, p_hat_, Ycal_, M_hat_)
    
    C <- outer(1:n_gamma, 1:n_gamma, Vectorize(function(r, l) {
        sum(vapply(1:n_clusters, function(i) {
          VARS$xi_[i,r] * u_star_[i,l] + VARS$xi_[i,l] * u_star_[i,r]
        }, 0))
    }))/n_clusters
    
    ################################################################## Step IV
    if (frailty == "none") {
      D_inv <- solve(jac[1:n_beta, 1:n_beta])
      COV <- (D_inv %*% (V[1:n_beta, 1:n_beta] 
                         + G[1:n_beta, 1:n_beta] 
                         + C[1:n_beta, 1:n_beta]) %*% t(D_inv))/n_clusters
      return(COV)
    }
    D_inv <- solve(jac)
    
    ################################################################## Step V
    
    COV <- (D_inv %*% (V + G + C) %*% t(D_inv))/n_clusters
    COV
  }
  
  if (frailty == "pvf") {
    theta_lower <- 1e-5
    theta_upper <- 1
  } else {
    theta_lower <- rep(1e-1, n_theta)
    theta_upper <- rep(Inf, n_theta)
  }
  
  # The actual optimization takes place here
  if (fitmethod == "loglik") {
    fitter <- optim(gamma_init, fit_fn,
                    gr=loglik_jacobian,
                    lower=c(rep(-Inf, n_beta),theta_lower), 
                    upper=c(rep(Inf, n_beta), theta_upper), 
                    method="L-BFGS-B",
                    control=list(factr=control$reltol/.Machine$double.eps, 
                                 pgtol=0, fnscale=-1)
                    )
    gamma_hat <- fitter$par
  } else if (control$fitmethod == "score") {
    fitter <- nleqslv(gamma_init, fit_fn, 
                      control=list(maxit=control$iter.max,
                                   xtol=0, ftol=control$abstol, btol=1e-3,
                                   allowSingular=TRUE),
                      method='Newton',
                      jac=score_jacobian,
                      jacobian=TRUE
                      )
    gamma_hat <- fitter$x
  }
  
  if (verbose)
    cat("Converged after", VARS$iter-1, "iterations\n")
  
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  
  Lambdafn <- Vectorize(function(t) {
    if (t <= 0) {
      return(0);
    }
    VARS$Lambda[sum(t >= time_steps)]
  })
  
  list(beta = beta_hat, #setNames(beta_hat, paste("beta.", names(beta_hat), sep="")),
       theta = setNames(theta_hat, paste("theta.", 1:length(theta_hat), sep="")),
       time_steps = time_steps,
       d_ = d_,
       Lambda = VARS$Lambda,
       Lambda.time=Lambda.time,
       loglik = VARS$loglik,
       frailty = frailty,
       Lambdafn = Lambdafn,
       fitter = fitter,
       fitmethod = control$fitmethod,
       score_jacobian = score_jacobian,
       fit_fn=fit_fn,
       loglik_jacobian=loglik_jacobian,
       VARS=VARS,
       iter=VARS$iter,
       frailty.variance=vfrailty[[frailty]](theta_hat),
       n.clusters=n_clusters
      )
}