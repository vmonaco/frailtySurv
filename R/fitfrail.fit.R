# Helper fn to get the rank, where duplicates have the same rank
increasing_rank <- function(d) {
  j <- unique(sort(d));
  return(sapply(d,function(dd) which(dd==j)));
}

fitfrail.fit <- function(x, y, cluster, beta_init, theta_init, frailty, control, rownames) {
  
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
  # time_steps <- c(0, sort(unique(time)))
  time_sorted_idx <- order(time)
  time_steps <- c(0, time[time_sorted_idx])
  k_tau <- length(time_steps) # index of the last time step
  
  # d_[k] is the number of failures at time t_k
  d_ <- c(0, unname(status[time_sorted_idx]))
  
  # Variables indexed by [[i]][j,] where i is cluster name and j is member
  X_ <- split.data.frame(x, cluster) # Covariates
  T_ <- split(time, cluster) # Followup time
  I_ <- split(status, cluster) # Failure indicator
  
  # Index of the observed time, K_[[i]][j]
  # K_ <- split(1 + increasing_rank(time), cluster) 
  K_ <- split(1 + rank(time, ties.method="first"), cluster)
  
  # Cluster names and sizes are needed in several places, so get these now
  cluster_names <- levels(cluster)
  names(cluster_names) <- cluster_names
  n_clusters <- length(cluster_names)

  cluster_sizes <- table(cluster)
  
  # Convenience function to apply a function over the cluster and members indices
  clustapply <- function(fn, FUN.VALUE) {
    lapply(cluster_names, function(i) {
      t(vapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      }, FUN.VALUE=FUN.VALUE))})
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
  
  # Several global vars are updated in the objective function: 
  # beta, theta, H_, H_dot_, Lambda_hat, lambda_hat, psi_
  # TODO: maybe use a different name, like coef, otherwise locked binding
  beta <- NULL
  
  iter <- 1
  verbose <- control$verbose
  fitmethod <- control$fitmethod
  # The objective function, either score equations or loglikelihood
  # where gamma is c(beta, theta)
  fit_fn <- function(gamma) {
    beta <<- gamma[1:n_beta]
    theta <<- gamma[(n_beta+1):(n_beta+n_theta)]
    
    R_ <<- clustapply(function(i, j) {
      exp(sum(beta * X_[[i]][j])) * Y_[[i]][j,]
    }, rep(0, k_tau))
    R_dot_ <<- lapply(R_, colSums)
    
    bh_ <- bh(d_, X_, K_, Y_, N_dot_, beta, theta, frailty)
    H_ <<- bh_$H_
    H_dot_ <<- bh_$H_dot
    lambda <<- bh_$lambda
    Lambda <<- bh_$Lambda
    psi_ <<- bh_$psi_
    phi_1_ <<- bh_$phi_1_
    phi_2_ <<- bh_$phi_2_
    
    # Update the loglik globally
    loglik <<- loglikelihood(X_, K_, I_, phi_1_, Lambda, beta)
    
    if (verbose) {
      if (iter == 1) {
        cat(c("Iter", 
              paste("beta_", 1:n_beta, sep=""), 
              paste("theta_", 1:n_theta, sep=""), 
              "loglik",
              "\n"), sep="\t")
      }
      cat(c(signif(c(iter, beta, theta, loglik), 4), "\n"), sep="\t")
      iter <<- iter + 1
    }
    
    if (fitmethod == "loglik") {
      loglik
    } else if (fitmethod == "score") {
      # TODO: This could be done in parallel, not much gain though.
      xi_beta_ <- vapply(1:n_beta, function(r) {
        xi_beta(X_, I_, N_dot_, H_, H_dot_, beta, theta, r, frailty)
      }, rep(0, n_clusters))
      
      xi_theta_ <- vapply(1:n_theta, function(r) {
        xi_theta(X_, I_, N_dot_, H_, H_dot_, beta, theta, r, frailty)
      }, rep(0, n_clusters))
      
      # xi_[i,r], where i is cluster, r is param idx
      xi_ <- cbind(xi_beta_, xi_theta_)
      
      colMeans(xi_)
    }
  }
  
  jacobian_ij <- function(i, j) {
    if (i <= n_beta && j <= n_beta) {
      dH <- dH_dbeta(d_, K_, X_, R_, R_dot_, N_dot_, H_, H_dot_, 
                            Lambda_hat, lambda_hat, beta, theta, j, frailty)
      
      jacobian_beta_beta(K_, X_, N_dot_, H_, H_dot_, dH$dH_dbeta_, dH$dH_dot_dbeta_,
                         beta, theta, i, j, frailty)
    } else if (i <= n_beta && j > n_beta) {
      dH <- dH_dtheta(d_, K_, X_, R_, R_dot_, N_dot_, H_, H_dot_, 
                      Lambda_hat, lambda_hat, beta, theta, n_beta - j, frailty)
      
      jacobian_beta_theta(K_, X_, N_dot_, H_, H_dot_, dH$dH_dtheta_, dH$dH_dot_dtheta_,
                         beta, theta, i, n_beta - j, frailty)
    } else if (i > n_beta && j <= n_beta) {
      dH <- dH_dbeta(d_, K_, X_, R_, R_dot_, N_dot_, H_, H_dot_, 
                     Lambda_hat, lambda_hat, beta, theta, j, frailty)
      
      jacobian_theta_beta(K_, X_, N_dot_, H_, H_dot_, dH$dH_dbeta_, dH$dH_dot_dbeta_,
                         beta, theta, n_beta - i, j, frailty)
    } else if (i > n_beta && j > n_beta) {
      dH <- dH_dtheta(d_, K_, X_, R_, R_dot_, N_dot_, H_, H_dot_, 
                      Lambda_hat, lambda_hat, beta, theta, n_beta - j, frailty)
      
      jacobian_theta_theta(K_, X_, N_dot_, H_, H_dot_, dH$dH_dtheta_, dH$dH_dot_dtheta_,
                          beta, theta, n_beta - i, n_beta - j, frailty)
    }
  }
  
  jacobian <- function(gamma=NULL) {
    cat('jac', gamma, '\n')
    
    # Build the jacobian matrix from the respective components
    outer(1:n_gamma, 1:n_gamma, Vectorize(function(i, j) 
      jacobian_ij(i, j)))
  }
  
  # The covariance matrix is encapsulated in a function since it takes some effort
  # to compute. By default, it is only computed when accessed by the user
  covariance <- function() {
    # First, gather some variables that were defined before the covar steps
    
    R_ <- clapply(cluster_names, cluster_sizes, function(i, j) {
      exp(sum(beta * X_[[i]][j])) * Y_[[i]][j,]
    })
    
    R_dot_ <- lapply(R_, function(R_i) colSums(R_i))
    
    # Sum over custers and members
    N_dot_dot <- colSums(Reduce('+', N_))

    ################################################################## Step I
    xi_beta_ <- vapply(1:n_beta, function(r) {
      xi_beta(X_, I_, N_dot, H_, H_dot_, beta, theta, r, frailty)
    }, rep(0, n_clusters))
    
    xi_theta_ <- vapply(1:n_theta, function(r) {
      xi_theta(X_, I_, N_dot, H_, H_dot_, beta, theta, r, frailty)
    }, rep(0, n_clusters))
    
    # xi_[i, r], where i is cluster, r is param idx
    xi_ <- cbind(xi_beta_, xi_theta_)
    
    ################################## V
    V <- Reduce('+', lapply(1:n_clusters, function(i) {
      outer(xi_[i,], xi_[i,])}
      ))/n_clusters
    
    ################################################################## Step II
    
    R_star <- clapply(cluster_names, cluster_sizes, function(i, j) {
      exp(sum(beta * X_[[i]][j, ]))
    })
    
    # Split Q up into beta and theta components, then create a Q_ function
    # that accesses component by index r
    Q_beta_ <- lapply(1:n_beta, function(r) {
      Q_beta(X_, R_star, N_dot, H_, H_dot_, theta, r, frailty)
    })
    
    Q_theta_ <- lapply((n_beta+1):(n_beta+n_theta), function(r) {
      Q_theta(X_, R_star, N_dot, H_, H_dot_, theta, r, frailty)
    })
    
    # Q_[[r]] has the same shape as H_dot_
    Q_ <- c(Q_beta_, Q_theta_)
    
    # calligraphic Y
    Ycal_ <- Ycal(X_, Y_, N_dot, H_dot_, beta, theta, frailty)
    
    # eta
    eta_ <- eta(N_dot, H_dot_, theta, frailty)
    
    # Upsilon
    Upsilon_ <- Upsilon(X_, K_, R_dot_, eta_, Ycal_, beta, theta, frailty)
    
    # Omega[[i]][j][s, t]
    Omega <- clapply2(function(i, j) {
      outer(1:tau, 1:tau, Vectorize(function(s, t) {
        integrate(Vectorize(function(u) {
          Ycal_[u]^(-2) * R_dot_[[i]][u] * eta_[[i]][u] * 
            exp(sum(beta*X_[[i]][j,])) * N_dot_dot[u]
        }), s, t, subdivisions = 1000)$value/n_clusters^2
      }))
    })
    
    p_hat <- vapply(1:tau, function(t) {
      prod(vapply(1:t, function(s) {
        1 + sum_clapply(function(i, j) {
          (I_[[i]][j] * Upsilon_[s] + Omega[[i]][[j]][s,t]) * N_tilde[[i]][j, s]
        })
      }, 0))
    }, 0)
    
    # TODO: correct?
    QN_sum <- lapply(1:n_gamma, function(r) {
      colSums(matrix(unlist(clapply2(function(i, j) {
        Q_[[r]][[i]][j,] * N_tilde[[i]][j,]
      })), ncol=tau))
    })
    
    pi_ <- outer(1:n_gamma, 1:tau, Vectorize(function(r, s) {
      integrate(Vectorize(function(t) {
        QN_sum[[r]][t]/p_hat[t]
      }), s, tau)$value/n_clusters
      
    }))
    
    G <- outer(1:n_gamma, 1:n_gamma, Vectorize(function(r, l) {
      integrate(Vectorize(function(s) {
        pi_[r, s] * pi_[l, s] * p_hat[s - 1]^2 * N_dot_dot[s]/Ycal_[s]^2
      }), 2, tau)$value/n_clusters
    }))
    
    ################################################################## Step III
    
    # M_[i,s]
    M_ <- lapply(1:n_clusters, function(i) {
      outer(1:cluster_sizes[[i]], 1:tau, Vectorize(function(j, t) {
        if (t <= 2) {
          N_[[i]][j, t]
        } else {
          N_[[i]][j, t] - integrate(Vectorize(function(u) {
            exp(sum(beta * X_[[i]][j,])) * Y_[[i]][j, u] * psi_[[i]][u - 1] * Lambda_hat[u]
          }), 2, t)$value
        }
      }))
    })
    
    Msum_ <- do.call(rbind, lapply(M_, colSums))
    
    # u_[i,r]
    u_ <- outer(1:n_clusters, 1:n_gamma, Vectorize(function(i, r) {
      integrate(Vectorize(function (s) {
        Msum_[i, s] * pi_[r, s] * p_hat[s - 1]/Ycal_[s]
      }), 2, tau)$value
    }))
    
    C <- outer(1:n_gamma, 1:n_gamma, Vectorize(function(r, l) {
        Reduce('+', vapply(1:n_clusters, function(i) {
          xi_[i,r] * u_[i,l] + xi_[i,l] * u_[i,r]
        }, 0))/n_clusters
    }))
    
    ################################################################## Step IV
    D <- jacobian()
    
    ################################################################## Step V
    D_inv <- solve(D)
    D_inv %*% (V + G + C) %*% t(D_inv)
  }
  
  # The actual optimization takes place here
  if (fitmethod == "loglik") {
    fitter <- optim(gamma_init, fit_fn, 
                    lower=c(-4,0.05), upper=c(4,10), method="L-BFGS-B",
                    control=list(fnscale=-1, factr=1e8, pgtol=1e-8)
                    # TODO: gradiant function
                    )
    gamma_hat <- fitter$par
  } else if (control$fitmethod == "score") {
    fitter <- nleqslv(gamma_init, fit_fn, 
                      control=list(maxit=control$iter.max,xtol=1e-8,ftol=1e-8,btol=1e-3),
                      method='Newton',
                      # jac=jacobian, # TODO
                      jacobian=TRUE
                      )
    gamma_hat <- fitter$x
  }
  
  if (verbose)
    cat("Converged after", iter, "iterations\n")
  
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  
  lambda_hat_nz <- lambda_hat[d_ > 0]
  time_steps_nz <- time_steps[d_ > 0]
  
  # Group lambda by time steps 
  lambda_hat_nz <- vapply(split(lambda_hat_nz, factor(time_steps_nz)), sum, 0)
  time_steps_nz <- unique(time_steps_nz)
  
  lambdafn <- Vectorize(function(t) {
    if (t <= 0) {
      return(0);
    }
    
    idx <- sum(t > time_steps_nz)
    
    if (idx == 0) {
      return(0)
    }
    
    if (idx >= length(time_steps_nz)) {
      idx <- length(time_steps_nz) - 1
    }
    lambda_hat_nz[idx]/(time_steps_nz[idx+1] - time_steps_nz[idx])
  })
  
  Lambdafn <- Vectorize(function(t) {
    if (t <= 0) {
      return(0);
    }
    Lambda_hat[sum(t > time_steps)]
  })
  
  list(beta = beta_hat,
       theta = theta_hat,
       time_steps = time_steps,
       lambda_hat_nz=lambda_hat_nz,
       time_steps_nz=time_steps_nz,
       d_ = d_,
       lambda = lambda_hat,
       Lambda = Lambda_hat,
       loglik = loglik,
       frailty = frailty,
       lambdafn = lambdafn,
       Lambdafn = Lambdafn,
       fitter = fitter,
       fitmethod = control$fitmethod,
       method = 'fitfrail',
       jacobian = jacobian,
       covariance = covariance
       # TODO: estimated frailty values
       # TODO: use numerical jac for variance?
      )
}