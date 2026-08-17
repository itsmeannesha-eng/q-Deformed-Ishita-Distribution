rm(list = ls())
set.seed(12345)
N <- 1000
sample_sizes <- c(25, 50, 75, 100, 200, 300, 600)

true_theta <- 1.00
true_q <- 1.10
C_qDI <- function(theta, q) {
  
  if (abs(q - 1) < 1e-8) {
    
    return(theta^3 / (theta^3 + 2))
    
  } else {
    
    numerator <-
      (2 - q) * (3 - 2*q) * (4 - 3*q) * theta^3
    
    denominator <-
      theta^3 * (3 - 2*q) * (4 - 3*q) + 2
    
    return(numerator / denominator)
  }
}

qexp_DI <- function(x, theta, q) {
  
  if (abs(q - 1) < 1e-8) {
    
    return(exp(-theta * x))
    
  }
  
  z <- 1 + (q - 1) * theta * x
  
  ans <- numeric(length(x))
  
  valid <- z > 0
  
  ans[valid] <- z[valid]^(-1 / (q - 1))
  
  ans[!valid] <- 0
  
  return(ans)
}

d_qDI <- function(x, theta, q) {
  
  if (theta <= 0 ||
      q >= 4/3 ||
      q >= 2) {
    return(rep(0, length(x)))
  }
  
  z <- 1 + (q - 1) * theta * x
  
  f <- numeric(length(x))
  
  valid <- x >= 0 & z > 0
  
  if (any(valid)) {
    
    f[valid] <-
      C_qDI(theta, q) *
      (theta + x[valid]^2) *
      z[valid]^(-1/(q - 1))
  }
  
  return(f)
}

logLik_qDI <- function(par, data) {
  
  theta <- par[1]
  q <- par[2]
  
  # Parameter restrictions
  if (theta <= 0 ||
      q >= 4/3 ||
      q <= -10) {
    return(-1e100)
  }
  
  z <- 1 + (q - 1) * theta * data
  
  # Invalid support
  if (any(z <= 0)) {
    return(-1e100)
  }
  
  C <- C_qDI(theta, q)
  
  if (!is.finite(C) || C <= 0) {
    return(-1e100)
  }
  
  logf <-
    log(C) +
    log(theta + data^2) -
    (1/(q - 1)) * log(z)
  
  # q = 1 case
  if (abs(q - 1) < 1e-8) {
    
    logf <-
      log(theta^3/(theta^3 + 2)) +
      log(theta + data^2) -
      theta * data
  }
  
  if (any(!is.finite(logf))) {
    return(-1e100)
  }
  
  return(sum(logf))
}

MLE_qDI <- function(data) {
  
  # Several starting values improve numerical stability
  starts <- list(
    c(0.5, 1.05),
    c(1.0, 1.05),
    c(1.5, 1.10),
    c(2.0, 0.90)
  )
  
  results <- vector("list", length(starts))
  
  for (i in seq_along(starts)) {
    
    results[[i]] <- try(
      optim(
        par = starts[[i]],
        fn = function(par)
          -logLik_qDI(par, data),
        method = "L-BFGS-B",
        lower = c(0.001, -5),
        upper = c(20, 1.32),
        control = list(
          maxit = 500,
          factr = 1e7
        )
      ),
      silent = TRUE
    )
  }
  
  #optimizations
  good <- sapply(
    results,
    function(x)
      !inherits(x, "try-error") &&
      is.finite(x$value)
  )
  
  if (!any(good)) {
    return(c(theta = NA, q = NA))
  }
  
  results_good <- results[good]
  
  #solution with largest likelihood
  values <- sapply(
    results_good,
    function(x) x$value
  )
  
  best <- results_good[[which.min(values)]]
  
  return(
    c(
      theta = best$par[1],
      q = best$par[2]
    )
  )
}

F_qDI <- function(x, theta, q) {
  
  if (abs(q - 1) < 1e-8) {
    
    # Classical Ishita CDF
    C <- theta^3 / (theta^3 + 2)
    
    return(
      C * (
        2/theta^2 -
          exp(-theta*x) *
          (theta^2 + theta*x^2 + 2*x + 2/theta)
      )
    )
  }
  
  C <- C_qDI(theta, q)
  
  a <- (q - 1) * theta
  y <- 1 + a*x
  
  # Support
  if (q < 1) {
    
    xmax <- 1/((1-q)*theta)
    
    if (x >= xmax) {
      return(1)
    }
  }
  
  if (x <= 0) {
    return(0)
  }
  
  if (y <= 0) {
    return(0)
  }
  
  A <- function(y) {
    
    term1 <-
      y^((q-2)/(q-1)) / (q-2)
    
    term2 <-
      1/((q-1)^2 * theta^3) *
      (
        y^((3*q-4)/(q-1))/(3*q-4)
        -
          2*y^((2*q-3)/(q-1))/(2*q-3)
        +
          y^((q-2)/(q-1))/(q-2)
      )
    
    return(term1 + term2)
  }
  
  if (q < 1) {
    
    upper <- 0
    
  } else {
    
    # q > 1
    upper <- 0
  }
  
  value <- C * (A(y) - A(1))
  
  # Numerical correction
  value <- max(0, min(1, value))
  
  return(value)
}
r_qDI <- function(n, theta, q) {
  
  u <- runif(n)
  x <- numeric(n)
  
  # Numerical inverse CDF
  for (i in seq_len(n)) {
    
    if (q < 1) {
      
      upper <- 1/((1-q)*theta)
      
    } else {
      
      # Large enough upper bound
      upper <- 100/theta
      
      while (
        F_qDI(upper, theta, q) < 0.999999
      ) {
        
        upper <- 2 * upper
        
        if (upper > 1e7) break
      }
    }
    
    x[i] <- uniroot(
      function(xx)
        F_qDI(xx, theta, q) - u[i],
      interval = c(0, upper),
      tol = 1e-8
    )$root
  }
  
  return(x)
}


run_MLE_simulation <- function(
    theta0,
    q0,
    n,
    N = 1000) {
  
  theta_hat <- numeric(N)
  q_hat <- numeric(N)
  
  AAD <- numeric(N)
  
  cat(
    "\nRunning: theta0 =",
    theta0,
    ", q0 =",
    q0,
    ", n =",
    n,
    "\n"
  )
  
  for (r in seq_len(N)) {
    
    # Generate sample
    data <- r_qDI(
      n = n,
      theta = theta0,
      q = q0
    )
    
    # MLE
    est <- MLE_qDI(data)
    
    theta_hat[r] <- est["theta"]
    q_hat[r] <- est["q"]
    
    # AAD
    if (is.finite(theta_hat[r]) &&
        is.finite(q_hat[r])) {
      
      F_hat <- sapply(
        data,
        function(x)
          F_qDI(
            x,
            theta_hat[r],
            q_hat[r]
          )
      )
      
      F_true <- sapply(
        data,
        function(x)
          F_qDI(
            x,
            theta0,
            q0
          )
      )
      
      AAD[r] <-
        mean(abs(F_hat - F_true))
      
    } else {
      
      AAD[r] <- NA
    }
  }
  
  # Remove failed estimates
  valid <- is.finite(theta_hat) &
    is.finite(q_hat)
  
  theta_hat <- theta_hat[valid]
  q_hat <- q_hat[valid]
  AAD <- AAD[valid]
  
  MSE_theta <-
    mean((theta_hat - theta0)^2)
  
  MRE_theta <-
    mean(abs(
      (theta_hat - theta0) / theta0
    ))
  
  Bias_theta <-
    mean(abs(theta_hat - theta0))
  
  
  MSE_q <-
    mean((q_hat - q0)^2)
  
  MRE_q <-
    mean(abs(
      (q_hat - q0) / q0
    ))
  
  Bias_q <-
    mean(abs(q_hat - q0))
  
  AAD_value <-
    mean(AAD, na.rm = TRUE)
  
  result <- data.frame(
    n = n,
    MSE_theta = MSE_theta,
    MRE_theta = MRE_theta,
    Bias_theta = Bias_theta,
    MSE_q = MSE_q,
    MRE_q = MRE_q,
    Bias_q = Bias_q,
    AAD = AAD_value
  )
  
  return(result)
}

results <- lapply(
  sample_sizes,
  function(n)
    run_MLE_simulation(
      theta0 = true_theta,
      q0 = true_q,
      n = n,
      N = N
    )
)

results_table <- do.call(
  rbind,
  results
)

print(results_table)

write.csv(
  results_table,
  "qDI_MLE_theta1_q1.10_simulation_results.csv",
  row.names = FALSE
)

