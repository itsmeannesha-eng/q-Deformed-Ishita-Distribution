rm(list = ls())
set.seed(12345)
N <- 1000
n_values <- c(25, 50, 75, 100, 200, 300, 600)
theta0 <- 1.0
q0 <- 1.1
Cq <- function(theta, q) {
  
  if (abs(q - 1) < 1e-8) {
    return(theta^3 / (theta^3 + 2))
  }
  
  I <- 1 / (theta * (2 - q)) +
    2 / (theta^3 *
           (2 - q) *
           (3 - 2*q) *
           (4 - 3*q))
  
  1 / I
}

dqDI <- function(x, theta, q) {
  
  if (theta <= 0)
    return(rep(0, length(x)))
  
  z <- 1 + (q - 1) * theta * x
  
  ans <- numeric(length(x))
  
  ok <- x >= 0 & z > 0
  
  ans[ok] <-
    Cq(theta, q) *
    (theta + x[ok]^2) *
    z[ok]^(-1/(q - 1))
  
  ans
}

r_qDI <- function(n, theta, q) {
  
  # Exponential 
  x <- rexp(n, rate = theta)
  
  # Rejection constant
  M <- 100
  
  result <- numeric(0)
  
  while (length(result) < n) {
    
    x <- rexp(n, rate = theta)
    
    u <- runif(n)
    
    ratio <-
      dqDI(x, theta, q) /
      (M * dexp(x, rate = theta))
    
    accept <- u < ratio
    
    result <- c(
      result,
      x[accept]
    )
  }
  
  result[1:n]
}

loglik <- function(par, x) {
  
  theta <- par[1]
  q <- par[2]
  
  # Parameter restrictions
  if (theta <= 0 || q <= 0)
    return(-1e100)
  
  # Avoid q = 1 singularity
  if (abs(q - 1) < 1e-5)
    return(-1e100)
  
  z <- 1 + (q - 1) * theta * x
  
  if (any(z <= 0))
    return(-1e100)
  
  C <- Cq(theta, q)
  
  logf <-
    log(C) +
    log(theta + x^2) -
    log(z)/(q - 1)
  
  if (any(!is.finite(logf)))
    return(-1e100)
  
  sum(logf)
}
MLE_qDI <- function(x, theta.start, q.start) {
  
  fit <- tryCatch({
    
    optim(
      par = c(theta.start, q.start),
      fn = function(par)
        -loglik(par, x),
      
      method = "L-BFGS-B",
      
      lower = c(0.01, 0.10),
      upper = c(10, 1.80),
      
      control = list(
        maxit = 1000
      )
    )
    
  }, error = function(e) NULL)
  
  
  if (is.null(fit))
    return(c(theta = NA, q = NA))
  
  if (!is.finite(fit$value))
    return(c(theta = NA, q = NA))
  
  c(
    theta = fit$par[1],
    q = fit$par[2]
  )
}

All_Results <- data.frame()

All_Estimates <- data.frame()

for (n in n_values) {
  
  cat("\n\n")
  cat("====================================================\n")
  cat("        SAMPLE SIZE n =", n, "\n")
  cat("====================================================\n")
  
  est <- matrix(
    NA,
    nrow = N,
    ncol = 2
  )
  
  colnames(est) <-
    c("theta_hat", "q_hat")
  
  
  cat("\nSimulation started for n =", n, "\n\n")
  
  for (i in 1:N) {
    
    # Generate sample
    x <- r_qDI(
      n = n,
      theta = theta0,
      q = q0
    )
    
    
    # Estimate theta and q
    est[i, ] <-
      MLE_qDI(
        x,
        theta.start = theta0,
        q.start = q0
      )
    
    if (i %% 100 == 0) {
      
      cat(
        "n =", n,
        ": Completed",
        i,
        "of",
        N,
        "\n"
      )
    }
  }
  
  valid <- complete.cases(est)
  
  est <- est[valid, ]
  
  R <- nrow(est)
  
  
  cat("\nSuccessful replications =", R, "\n")
  
  theta_hat <- est[, "theta_hat"]
  q_hat <- est[, "q_hat"]
  
  Bias_theta <-
    mean(theta_hat) - theta0
  
  Bias_q <-
    mean(q_hat) - q0
  
  MSE_theta <-
    mean(
      (theta_hat - theta0)^2
    )
  
  MSE_q <-
    mean(
      (q_hat - q0)^2
    )
  
  MRE_theta <-
    mean(
      abs(theta_hat - theta0) /
        theta0
    )
  
  MRE_q <-
    mean(
      abs(q_hat - q0) /
        q0
    )
  
  AAD_theta <-
    mean(
      abs(theta_hat - theta0)
    )
  
  AAD_q <-
    mean(
      abs(q_hat - q0)
    )
  
  AAD <-
    mean(
      c(
        abs(theta_hat - theta0),
        abs(q_hat - q0)
      )
    )
  
  Results_n <- data.frame(
    
    n = n,
    
    Parameter = c(
      "theta",
      "q"
    ),
    
    True = c(
      theta0,
      q0
    ),
    
    Mean = c(
      mean(theta_hat),
      mean(q_hat)
    ),
    
    Bias = c(
      Bias_theta,
      Bias_q
    ),
    
    MSE = c(
      MSE_theta,
      MSE_q
    ),
    
    MRE = c(
      MRE_theta,
      MRE_q
    ),
    
    AAD = c(
      AAD_theta,
      AAD_q
    )
  )
  
  
  # Add to overall results
  All_Results <-
    rbind(
      All_Results,
      Results_n
    )
  
  est_df <- data.frame(
    n = n,
    theta_hat = theta_hat,
    q_hat = q_hat
  )
  
  All_Estimates <-
    rbind(
      All_Estimates,
      est_df
    )
  
  cat("\n\nResults for n =", n, "\n")
  
  print(
    Results_n,
    digits = 6,
    row.names = FALSE
  )
  
  cat(
    "\nOverall AAD for n =",
    n,
    ":",
    AAD,
    "\n"
  )
}

cat("\n\n")
cat("============================================================\n")
cat("          FINAL MONTE CARLO SIMULATION RESULTS\n")
cat("============================================================\n\n")


print(
  All_Results,
  digits = 6,
  row.names = FALSE
)

Theta_Results <-
  All_Results[
    All_Results$Parameter == "theta",
  ]

cat("\n\n============================================================\n")
cat("                    THETA RESULTS\n")
cat("============================================================\n\n")

print(
  Theta_Results,
  digits = 6,
  row.names = FALSE
)

Q_Results <-
  All_Results[
    All_Results$Parameter == "q",
  ]

cat("\n\n============================================================\n")
cat("                       q RESULTS\n")
cat("============================================================\n\n")

print(
  Q_Results,
  digits = 6,
  row.names = FALSE
)

write.csv(
  All_Results,
  "qDI_MLE_All_n_Results.csv",
  row.names = FALSE
)

write.csv(
  All_Estimates,
  "qDI_MLE_All_n_Estimates.csv",
  row.names = FALSE
)

write.csv(
  Theta_Results,
  "qDI_MLE_Theta_Results.csv",
  row.names = FALSE
)

write.csv(
  Q_Results,
  "qDI_MLE_q_Results.csv",
  row.names = FALSE
)

cat("\n\n============================================================\n")
cat("Simulation completed for all sample sizes.\n")
cat("Sample sizes:", n_values, "\n")
cat("Results saved successfully.\n")
cat("============================================================\n")

