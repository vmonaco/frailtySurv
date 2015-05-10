
dgamma_ <- function(x, theta) {
  dgamma(x, 1/theta, 1/theta)
}

dgamma_dtheta <- function(x, theta) {
  if (is.infinite(x)) return(0)
  
  term1 <- (x/theta)^(1/theta - 1)
  term2 <- exp(-x/theta)
  term3 <- log(theta/x) + digamma(1/theta) + x - 1
  
  (term1 * term2 * term3)/(gamma(1/theta)*theta^3)
}
# 
# dlnorm_dtheta <- function(x, sigma) {
#   term1 <- log(x)^2 * exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^4)
#   term2 <- exp(-(log(x)^2)/(2*sigma^2)) / (x*sqrt(2*pi)*sigma^2)
#     
#   term1 - term2
# }

sum_by_family <- function(A_) {
  A_dot <- sapply(names(A_), function(i) {
    rep(0, length(A_[[i]][1,]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  for (i in names(A_)) {
    A_dot[[i]] <- colSums(A_[[i]])
  }
  
  A_dot
}

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
  
  # Store the time step index to avoid lookups when the rank is needed later
  data$k <- 1 + rank(data[[timecol]])
  
  # The time steps will be used throughout this function
  time_steps <- c(0, sort(unique(data$time)))
  families <- split(data, data$family)
  family_names <- names(families)
  tau_k <- length(time_steps)
  
  # The number of individuals in family i that failed up to time step t, N[i,t]
  # N_dot[[i]][k]
  N_ <- sapply(family_names, function(i) {
    t(sapply(1:nrow(families[[i]]), function(j) {
      as.numeric((families[[i]][j,statuscol] > 0) & (families[[i]][j,timecol] <= time_steps))
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)

  N_dot <- sum_by_family(N_)
  
  # Y_[[i]][j,k] = I(T_ij >= t_k)
  Y_ <- sapply(families, function(family) {
    t(apply(family, 1, function(individual) {
      as.numeric(individual[timecol] >= time_steps)
    }))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # d_[k] = The number of failures at time tau_k
  d_ <- sapply(seq_along(time_steps), function(k) {
    sum(data[[statuscol]][data[[timecol]] == time_steps[k]])
  })
  
  # Create psi_i for each family. Depends on theta, t, and H_i.(t)
  psi_ <- sapply(family_names, function (i) {
    function(theta, H_dot, k) {
      
      phi_1 <- function (w) {
        ( w ^ N_dot[[i]][k] ) * exp(-w*H_dot[[i]][k]) * dgamma_(w, theta)
      }
      
      phi_2 <- function (w) {
        ( w ^ (N_dot[[i]][k] + 1) ) * exp(-w*H_dot[[i]][k]) * dgamma_(w, theta)
      }
      # For gamma, 
      # psi_[i](theta, k, H_dot_ik) == (N_dot[[i]][k] + 1/theta)/(H_dot[[i]][k] + 1/theta)
      integrate(Vectorize(phi_2), 0, Inf, stop.on.error = FALSE)$value/
        integrate(Vectorize(phi_1), 0, Inf, stop.on.error = FALSE)$value
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Step function for baseline hazard
  estimate_lambda_hat <- function(beta, theta) {
    # delta_lambda_hat[1] = 0, since at t_0 there is no chance of failure
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
          psi_[[i]](theta, H_dot, k - 1) * 
          sum(sapply(1:nrow(families[[i]]), function(j) {
            Y_[[i]][j, k] * exp(t(beta) %*% families[[i]][j,covars])
          }))
      }))
      
      delta_lambda_hat[k] <- d_[k]/denom
    
      # Update the cumsum since it's needed several times in the next iteration
      lambda_hat[k] <- lambda_hat[k-1] + delta_lambda_hat[k]
      
      # Determine H_ij(t) and H_i.(t)
      for (i in family_names) {
        for (j in 1:nrow(families[[i]])) {
          k_min = min(families[[i]][j,"k"], k)
          H_[[i]][j, k] <- lambda_hat[k_min] * exp(t(beta) %*% families[[i]][j,covars])
        }
        H_dot[[i]][k] <- sum(H_[[i]][,k])
      }
    }
    
    list(lambda_hat=lambda_hat, H_=H_, H_dot=H_dot)
  }
  
  # TODO: Build the system of equations
  # E.g. will look something like this, without the duplicate integrations
  U <- function(beta, theta) {
    estimator <- estimate_lambda_hat(beta, theta)
    lambda_hat <- estimator$lambda_hat
    H_ <- estimator$H_
    H_dot <- estimator$H_dot
    
    U_r <- function(covarcol) {
      term1 <- mean(sapply(family_names, function(i) {
          sum(sapply(1:nrow(families[[i]]), function(j) {
            as.numeric(families[[i]][j,statuscol]) * families[[i]][j,covarcol]
          }))
      }))
      
      term2 <- mean(sapply(family_names, function(i) {
        sum(sapply(1:nrow(families[[i]]), function(j) {
          H_[[i]][j, families[[i]][j,"k"]] * families[[i]][j,covarcol]
        })) *
        psi_[[i]](theta, H_dot, tau_k)
      }))
      term1 - term2
    }
    
    U_p <- function() {
      mean(sapply(family_names, function(i) {
        numer = integrate(function (w) {
          ( w ^ N_dot[[i]][tau_k] ) * exp(-w*H_dot[[i]][tau_k]) * Vectorize(function(wi) dgamma_dtheta(wi, theta))(w)
        }, 0, Inf, stop.on.error = FALSE)$value

        denom = integrate(function (w) {
          ( w ^ N_dot[[i]][tau_k] ) * exp(-w*H_dot[[i]][tau_k]) * dgamma_(w, theta)
        }, 0, Inf, stop.on.error = FALSE)$value
        
        numer/denom
      }))
    }
  
    c(U_r("Z1"), U_p())
  }
  
  function (x) {
    beta <- x[1]
    theta <- x[2]
    U(beta, theta)
  }
}
