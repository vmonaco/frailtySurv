
sum_by_cluster <- function(A_) {
  A_dot <- sapply(names(A_), function(i) {
    rep(0, length(A_[[i]][1,]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  for (i in names(A_)) {
    A_dot[[i]] <- as.double(colSums(A_[[i]]))
  }
  
  A_dot
}

clapply <- function(clnames, clsizes, fn) {
  sapply(clnames, function(i) {
    t(sapply(1:clsizes[[i]], function(j) {
      fn(i, j)
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
}

# Helper fn to get the rank, where duplicates have the same rank
increasing_rank <- function(d) {
  j <- unique(sort(d));
  return(sapply(d,function(dd) which(dd==j)));
}

fitfrail.fit <- function(x, y, cluster, beta_init, theta_init, frailty, control, rownames) {
  
  # Num observations
  n <- nrow(y)
  
  # Num coefs
  if (is.matrix(x)) {
    n_beta <- ncol(x)
  } else {
    if (length(x)==0) n_beta <-0
    else n_beta <- 1
  }
  
  # Num frailty distribution params
  n_theta <- length(theta_init)
  
  # total num params
  n_gamma <- n_beta + n_theta
  
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
  # time_steps <- c(0, sort(unique(time)))
  
  time_steps <- c(0, sort(time))
  time_steps_status <- c(0, unname(status[order(time)]))
  tau <- length(time_steps)
  
  # The _ variables are all lists, indexed by the cluster name
  X_ <- split.data.frame(x, cluster) # Covariates
  T_ <- split(time, cluster) # Followup time
  I_ <- split(status, cluster) # Failure indicator
  # Failure rank + 1, where time_steps[K_[i][j]] is the failure time of member j in cluster i
  # K_ <- split(1 + increasing_rank(time), cluster) 
  K_ <- split(1 + rank(time, ties.method="first"), cluster)
  
  # Cluster names and sizes are needed in several places, so get these now
  cluster_names <- levels(cluster)
  names(cluster_names) <- cluster_names
  n_clusters <- length(cluster_names)

  cluster_sizes <- table(cluster)
  
  clapply2 <- function(fn) {
    lapply(cluster_names, function(i) {
      lapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      })
    })
  }
  
  sum_clapply <- function(fn) {
    sum(unlist(clapply(cluster_names, cluster_sizes, fn)))
  }
  
  # N_[[i]][j, k] indicates whether individual j in cluster i failed at time <= tau_k
#   N_ <- sapply(cluster_names, function(i) {
#     t(sapply(1:cluster_sizes[[i]], function(j) {
#       # as.numeric((I_[[i]][j] > 0) & (T_[[i]][j] <= time_steps))
#       as.numeric((I_[[i]][j] > 0) & (K_[[i]][j] <= 1:length(time_steps)))
#     }))
#   }, simplify = FALSE, USE.NAMES = TRUE)
  N_ <- clapply(cluster_names, cluster_sizes, function(i, j) {
    as.numeric((I_[[i]][j] > 0) & (K_[[i]][j] <= 1:length(time_steps)))})
  
  # N_dot[[i]][k] is the number of individuals in cluster i that failed up to time tau_k
  N_dot <- sum_by_cluster(N_)
  
  N_tilde <- clapply(cluster_names, cluster_sizes, function(i, j) {
    as.numeric((K_[[i]][j] <= 1:length(time_steps)))})
  
  # Y_[[i]][j,k] indicates whether individual j in cluster i failed at time >= tau_k
#   Y_ <- sapply(cluster_names, function(i) {
#     t(sapply(1:cluster_sizes[[i]], function(j) {
#       # as.numeric(T_[[i]][j] >= time_steps)
#       as.numeric(K_[[i]][j] >= 1:length(time_steps))
#     }))
#   }, simplify = FALSE, USE.NAMES = TRUE)
  Y_ <- clapply(cluster_names, cluster_sizes, function(i, j) {
      as.numeric(K_[[i]][j] >= 1:length(time_steps))})
  
  d_ <- time_steps_status
  # d_[k] is the number of failures at time tau_k
#   d_ <- vapply(seq_along(time_steps), function(k) {
#     sum(status[time == time_steps[k]])
#   }, 0)
  
  verbose <- control$verbose
  if (verbose)
    cat(c("Iter", paste("beta_", 1:n_beta, sep=""), paste("theta_", 1:n_theta, sep=""), "\n"), sep="\t")
  
  # TODO: maybe use a different name, like coef, otherwise locked binding
  beta <- NULL
  
  iter <- 0
  # Build the system of equations where gamma is c(beta, theta)
  U <- function(gamma) {
    beta <<- gamma[1:n_beta]
    theta <<- gamma[(n_beta+1):(n_beta+n_theta)]
    
    estimator <- baseline_hazard_estimator(X_, K_, d_, Y_, N_dot, beta, theta, frailty)
    H_ <<- estimator$H_
    H_dot_ <<- estimator$H_dot
    lambda_hat <<- estimator$lambda_hat
    Lambda_hat <<- estimator$Lambda_hat
    
    if (verbose) {
      iter <<- iter + 1
      cat(c(signif(c(iter, beta, theta), 4), "\n"), sep="\t")
    }
    
    # TODO: This could be done in parallel, not much gain though.
    c(sapply(1:n_beta,
             function(beta_idx)
               U_r(X_, K_, I_, N_dot, H_, H_dot_, 
                   beta, theta, beta_idx, frailty)),
      sapply(1:n_theta, 
             function(theta_idx) 
               U_p(X_, K_, I_, N_dot, H_, H_dot_, 
                   beta, theta, theta_idx, frailty)))/n_clusters
  }
  
  loglikfn <- function(beta, theta, frailty) {
    if (verbose) {
      iter <<- iter + 1
      cat(c(signif(c(iter, beta, theta), 4), "\n"), sep="\t")
    }
    estimator <- baseline_hazard_estimator(X_, K_, d_, Y_, N_dot, beta, theta, frailty)
    loglikelihood(X_, K_, I_, N_dot, estimator$H_dot, estimator$Lambda_hat, beta, theta, frailty)
  }
  
  jacobian_ij <- function(i, j, beta, theta, R_, R_dot_, H_, H_dot_, Lambda, lambda) {
    if (i <= n_beta && j <= n_beta) {
      dH <- dH_dbeta(d_, K_, X_, R_, R_dot_, N_dot, H_, H_dot_, 
                            Lambda, lambda, beta, theta, j, frailty)
      
      jacobian_beta_beta(K_, X_, N_dot, H_, H_dot, dH$dH_dbeta_, dH$dH_dot_dbeta_,
                         beta, theta, i, j, frailty)
    } else if (i <= n_beta && j > n_beta) {
      dH <- dH_dtheta(d_, K_, X_, R_, R_dot_, N_dot, H_, H_dot_, 
                     Lambda, lambda, beta, theta, n_beta - j, frailty)
      
      jacobian_beta_theta(K_, X_, N_dot, H_, H_dot, dH$dH_dtheta_, dH$dH_dot_dtheta_,
                         beta, theta, i, n_beta - j, frailty)
    } else if (i > n_beta && j <= n_beta) {
      dH <- dH_dbeta(d_, K_, X_, R_, R_dot_, N_dot, H_, H_dot_, 
                     Lambda, lambda, beta, theta, j, frailty)
      
      jacobian_theta_beta(K_, X_, N_dot, H_, H_dot, dH$dH_dbeta_, dH$dH_dot_dbeta_,
                         beta, theta, n_beta - i, j, frailty)
    } else if (i > n_beta && j > n_beta) {
      dH <- dH_dtheta(d_, K_, X_, R_, R_dot_, N_dot, H_, H_dot_, 
                      Lambda, lambda, beta, theta, n_beta - j, frailty)
      
      jacobian_theta_theta(K_, X_, N_dot, H_, H_dot, dH$dH_dtheta_, dH$dH_dot_dtheta_,
                          beta, theta, n_beta - i, n_beta - j, frailty)
    }
  }
  
  jacobian <- function(gamma) {
    n_gamma <- length(gamma)
    stopifnot(length(gamma) == n_beta + n_theta)
    
    beta <- gamma[1:n_beta]
    theta <- gamma[(n_beta+1):(n_beta+n_theta)]
    
    R_ <- clapply(cluster_names, cluster_sizes, function(i, j) {
      exp(sum(beta * X_[[i]][j])) * Y_[[i]][j,]
    })
    R_dot_ <- sum_by_cluster(R_)
    
    bh <- baseline_hazard_estimator(X_, K_, d_, Y_, N_dot, beta, theta, frailty)
    Lambda <- bh$Lambda_hat
    lambda <- bh$lambda_hat
    H_ <- bh$H_
    H_dot_ <- bh$H_dot
    
    # Build the jacobian matrix from the respective components
    outer(1:n_gamma, 1:n_gamma, Vectorize(function(i, j) 
      jacobian_ij(i, j, beta, theta, R_, R_dot_, H_, H_dot_, Lambda, lambda)))/n_clusters
  }
  
  est <- function(beta, theta, frailty="gamma") {
    baseline_hazard_estimator(X_, K_, d_, Y_, N_dot, beta, theta, frailty)
  }
  
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
      }), 1, tau)$value/n_clusters
    }))
    
    ################################################################## Step III
    
    u_
    
    M_
    
    C <- outer(1:n_gamma, 1:n_gamma, Vectorize(function(r, l) {
        Reduce('+', vapply(1:n_clusters, function(i) {
          xi_[i,r] * u_[i,l] + xi_[i,l] * u_[i,r]
        }, 0))/n_clusters
    }))
    
    ################################################################## Step IV
    D <- jacobian(gamma)
    
    ################################################################## Step V
    solve(D)*(V + G + C)*solve(t(D))
  }
  
  if (control$fitmethod == 'score') {
    fitter <- nleqslv(c(beta_init, theta_init), U, 
                      control=list(maxit=control$iter.max,xtol=1e-8,ftol=1e-8,btol=1e-3),
                      # jac=jacobian,
                      jacobian=TRUE
                      )
    gamma_hat <- fitter$x
  } else if (control$fitmethod == 'like') {
    stop("TODO")
    fitter <- optim(c(0, 0.01), function(x) loglikfn(x[1], x[2], frailty), 
                 lower=c(-4,0.05), upper=c(4,10), method="L-BFGS-B",
                 control=list(fnscale=-1, factr=1e8, pgtol=1e-8))
    gamma_hat <- fitter$par
  }
  
  if (verbose)
    cat("Converged after", iter, "iterations\n")
  
  beta_hat <- gamma_hat[1:n_beta]
  theta_hat <- gamma_hat[(n_beta+1):(n_beta+n_theta)]
  loglik <- loglikfn(beta_hat, theta_hat, frailty)
  
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
       loglikfn=loglikfn,
       lambdafn = lambdafn,
       Lambdafn = Lambdafn,
       fitter = fitter,
       fitmethod = control$fitmethod,
       method = 'fitfrail',
       jacobian = jacobian,
       U = U,
       est = est,
       covariance = covariance
       # TODO: estimated frailty values
       # TODO: use numerical jac for variance?
      )
}