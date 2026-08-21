rm(list = ls())

qDI_kernel <- function(x, theta, q) {
  
  if (theta <= 0 || q <= 0) {
    return(rep(0, length(x)))
  }
  
  # q > 1 : support is (0, infinity)
  if (q > 1) {
    
    z <- 1 + (q - 1) * theta * x
    
    return(
      ifelse(
        z > 0,
        (theta + x^2) * z^(-1 / (q - 1)),
        0
      )
    )
  }
  
  if (q < 1) {
    
    upper <- 1 / ((1 - q) * theta)
    
    z <- 1 + (q - 1) * theta * x
    
    return(
      ifelse(
        x >= 0 & x <= upper & z > 0,
        (theta + x^2) * z^(-1 / (q - 1)),
        0
      )
    )
  }
  
  return(
    (theta + x^2) * exp(-theta * x)
  )
}


C_qDI <- function(theta, q) {
  if (q > 1) {
    
    integ <- integrate(
      function(x)
        qDI_kernel(x, theta, q),
      lower = 0,
      upper = Inf,
      subdivisions = 1000,
      rel.tol = 1e-9
    )$value
  }
  else if (q < 1) {
    
    upper <- 1 / ((1 - q) * theta)
    
    integ <- integrate(
      function(x)
        qDI_kernel(x, theta, q),
      lower = 0,
      upper = upper,
      subdivisions = 1000,
      rel.tol = 1e-9
    )$value
  }
  else {
    
    integ <- integrate(
      function(x)
        qDI_kernel(x, theta, q),
      lower = 0,
      upper = Inf,
      subdivisions = 1000,
      rel.tol = 1e-9
    )$value
  }
  
  return(1 / integ)
}

Tsallis_entropy <- function(theta, q, alpha) {
  
  if (theta <= 0) {
    return(NA)
  }
  
  if (alpha <= 0 || alpha == 1) {
    return(NA)
  }
  
  C <- C_qDI(theta, q)
  
  
  if (q > 1) {
    
    integral_value <- integrate(
      function(x) {
        
        f <- C * qDI_kernel(x, theta, q)
        
        f^alpha
      },
      lower = 0,
      upper = Inf,
      subdivisions = 1000,
      rel.tol = 1e-8
    )$value
  }
  
  
  else if (q < 1) {
    
    upper <- 1 / ((1 - q) * theta)
    
    integral_value <- integrate(
      function(x) {
        
        f <- C * qDI_kernel(x, theta, q)
        
        f^alpha
      },
      lower = 0,
      upper = upper,
      subdivisions = 1000,
      rel.tol = 1e-8
    )$value
  }
  
  else {
    
    integral_value <- integrate(
      function(x) {
        
        f <- C * qDI_kernel(x, theta, q)
        
        f^alpha
      },
      lower = 0,
      upper = Inf,
      subdivisions = 1000,
      rel.tol = 1e-8
    )$value
  }
  
  H <- (1 - integral_value) / (alpha - 1)
  
  return(H)
}

theta_values <- seq(0.1, 2.0, by = 0.1)

q <- 1.30
alpha <- 8

entropy_values <- sapply(
  theta_values,
  function(theta)
    Tsallis_entropy(
      theta = theta,
      q = q,
      alpha = alpha
    )
)

results <- data.frame(
  theta = theta_values,
  Tsallis_Entropy = entropy_values
)

print(results)
plot(
  theta_values,
  entropy_values,
  type = "o",
  pch = 19,
  lwd = 2,
  xlab = expression(theta),
  ylab = expression(H[alpha](X)),col="green"
)

grid()





