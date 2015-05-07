
#' Fit a Shared Frailty model
#'
#' \code{fitfrail} fits a model.
#'
#' This function ...
#' 
#' @param ... model
#' @return parameters
#'   \url{} for more details.
#' @examples
#' fitfrail(genfrail())
#' 
#' @export
fitfrail <- function(data) {
  
  # Temporary, these values will be extractd from the observed data based on
  # the survival formula
  timecol <- 'time'
  statuscol <- 'status'
  covars <- c('Z1')
  
  # Store the rank to avoid finding the min between t and T_ij
  data$rank <- order(data[[timecol]])
  
  # The time steps will be used throughout this function
  time_steps <- sort(unique(data$time))
  families <- split(data, data$family)
  family_names <- names(families)
  
  # The number of individuals in family i that failed up to time step t, N[i,t]
  # N_[[i]][k]
  N_ <- sapply(families, function(family) {
    colSums(outer(family[,c(timecol,statuscol)], time_steps, function(row, t_k) {
      (as.numeric(row[[statuscol]]) > 0) & (row[[timecol]] <= t_k)
    }))    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Y_[[i]][j,k] = I(T_ij >= t_k)
  Y_ <- sapply(families, function(family) {
    t(apply(family, 1, function(individual) {
      as.numeric(individual[timecol] >= time_steps)
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # d_[k] = The number of failures at time tau_k
  d_ <- sapply(seq_along(time_steps), function(k) {
    sum(data$status[data$time == time_steps[k]])
  })
  
  # Create psi_i for each family. Depends on theta, t, and H_i.(t)
  psi_ <- sapply(family_names, function (i) {
    function(theta, k, H_dot_ik) {
      
      phi_1 <- function (w) {
        ( w ^ N_[[i]][k] ) * exp(-w*H_dot_ik) * dgamma(w, theta, theta)
      }
      
      phi_2 <- function (w) {
        ( w ^ (N_[[i]][k] + 1) ) * exp(-w*H_dot_ik) * dgamma(w, theta, theta)
      }
      
      integrate(phi_2, 0, Inf)$value/integrate(phi_1, 0, Inf)$value
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Step function for baseline hazard
  estimate_lambda_hat <- function(beta, theta) {
    # delta_lambda_hat[1] = 0, since at t_0 there are no failures
    delta_lambda_hat = rep(0, length(time_steps))
    lambda_hat = rep(0, length(time_steps))
    
    H_ <- sapply(family_names, function(i) {
      matrix(0, nrow(families[[i]]), length(time_steps))
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    H_dot <- sapply(family_names, function(i) {
      rep(0, length(time_steps))
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    for (k in 2:length(time_steps)) {
      # Denom for delta_lambda_hat[k]
      denom <- sum(sapply(family_names, function(i) {
          psi_[[i]](theta, k, H_dot[[i]][k]) * 
          sum(sapply(1:nrow(families[[i]]), function(j) {
            Y_[[i]][j, k] * exp(t(beta) %*% families[[i]][j,covars])
          }))
      }))
      
      delta_lambda_hat[k] <- d_[k]/denom
    
      # Update the cumsum since it's needed several times in the next iteration
      lambda_hat[k] <- sum(delta_lambda_hat[1:k])
      
      # Determine H_ij(t) and H_i.(t)
      for (i in family_names) {
        for (j in 1:nrow(families[[i]])) {
          k_min = min(families[[i]][j,"rank"], k)
          H_[[i]][j, k] <- lambda_hat[k_min] * exp(t(beta) %*% families[[i]][j,covars])
        }
        H_dot[[i]][k] <- sum(H_[[i]][,k])
      }
    }
    
    list(lambda_hat=lambda_hat, H_=H_, H_dot=H_dot)
  }
  
  # TODO: Build the system of equations
  # E.g. will look something like this, without the duplicate integrations
  U_r <- function(beta, theta, covarcol) {
    estimator <- estimate_lambda_hat(beta, theta)
    lambda_hat <- estimator$lambda_hat
    H_ <- estimator$H_
    H_dot <- estimator$H_dot
    
    tau <- length(time_steps)
    term1 <- sum(data[,statuscol]*data[,covarcol])
    
    term2 <- sum(sapply(family_names, function(i) {
      psi_[[i]](theta, tau, H_dot[[i]][tau]) * 
        sum(sapply(1:nrow(families[[i]]), function(j) {
          H_[[i]][j, families[[i]][j,"rank"]] * families[[i]][j,covarcol]
        }))
    }))
    
    term1 - term2
  }
  
  U_p <- function(beta, theta) {
    estimator <- estimate_lambda_hat(beta, theta)
    lambda_hat <- estimator$lambda_hat
    H_ <- estimator$H_
    H_dot <- estimator$H_dot
    
    tau <- length(time_steps)
    
    sum(sapply(family_names, function(i) {
      integrate(function (w) {
        # TODO: df(w)/dtheta correct? Getting a divergent integral with large theta right now.
        ( w ^ N_[[i]][tau] ) * exp(-w*H_dot[[i]][tau]) * sapply(w, function (wi) {grad(function(p) {dgamma(wi, 1/p, 1/p)}, x=theta)})
      }, 0, Inf)$value /
        integrate(function (w) {
        ( w ^ N_[[i]][tau] ) * exp(-w*H_dot[[i]][tau]) * dgamma(w, 1/theta, 1/theta)
      }, 0, Inf)$value
    }))
  }
  
  U_r
}

