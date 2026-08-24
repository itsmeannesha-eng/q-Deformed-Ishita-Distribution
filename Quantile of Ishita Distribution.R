# q-Deformed Ishita PDF

dqDI <- function(x, theta, q) {
  
  C <- ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q) + 2)
  
  C * (theta + x^2) *
    (1 + (q-1)*theta*x)^(-1/(q-1))
}


# CDF obtained by numerical integration
pDI <- function(x, theta, q) {
  
  sapply(x, function(xx) {
    integrate(
      function(t) dqDI(t, theta, q),
      lower = 0,
      upper = xx
    )$value
  })
}


# Quantile function
qDI <- function(p, theta, q) {
  
  sapply(p, function(pp) {
    
    uniroot(
      function(x) pDI(x, theta, q) - pp,
      interval = c(0, 1000)
    )$root
    
  })
}


# Parameters
theta <- 0.5
q <- 1.05

# Required probabilities
p <- c(0.25, 0.50, 0.75, 0.90, 0.95)

# Numerical quantiles
quantiles <- qDI(p, theta, q)

# Display results
data.frame(
  p = p,
  Quantile = quantiles
)

