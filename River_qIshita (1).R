data("rivers")
force(rivers)
x <- rivers / 100  
n <- length(x)


#(a) q-Deformed Ishita (qDI)
dqDI <- function(x, theta, q) {
  if (theta <= 0 || q >= 4/3 || q < 1) return(rep(0, length(x)))
  num <- (2 - q) * (3 - 2 * q) * (4 - 3 * q) * (theta^3)
  den <- 2 + (3 - 2 * q) * (4 - 3 * q) * (theta^3)
  if (den <= 0 || num <= 0) return(rep(0, length(x)))
  
  C_val <- num / den
  pdf_vals <- numeric(length(x))
  valid <- (x > 0)
  
  if (abs(q - 1) < 1e-5) {
    pdf_vals[valid] <- (theta^3 / (theta^3 + 2)) * (theta + x[valid]^2) * exp(-theta * x[valid])
  } else {
    term <- 1 + (q - 1) * theta * x[valid]
    pdf_vals[valid] <- C_val * (theta + x[valid]^2) * (term^(-1 / (q - 1)))
  }
  return(pdf_vals)
}

pqDI <- function(x, theta, q) {
  sapply(x, function(val) {
    if (val <= 0) return(0)
    tryCatch(integrate(dqDI, lower = 0, upper = val, theta = theta, q = q)$value, error = function(e) NA)
  })
}

# (b) q-Deformed Lindley 
dqDL <- function(x, theta, q) {
  if (theta <= 0 || q >= 1.5 || q < 1) return(rep(0, length(x)))
  num <- (2 - q) * (3 - 2 * q) * (theta^2)
  den <- 1 + (3 - 2 * q) * theta
  if (den <= 0 || num <= 0) return(rep(0, length(x)))
  
  C_val <- num / den
  pdf_vals <- numeric(length(x))
  valid <- (x > 0)
  
  if (abs(q - 1) < 1e-5) {
    pdf_vals[valid] <- (theta^2 / (1 + theta)) * (1 + x[valid]) * exp(-theta * x[valid])
  } else {
    term <- 1 + (q - 1) * theta * x[valid]
    pdf_vals[valid] <- C_val * (1 + x[valid]) * (term^(-1 / (q - 1)))
  }
  return(pdf_vals)
}

pqDL <- function(x, theta, q) {
  sapply(x, function(val) {
    if (val <= 0) return(0)
    tryCatch(integrate(dqDL, lower = 0, upper = val, theta = theta, q = q)$value, error = function(e) NA)
  })
}

# (c) q-Exponential 
dqExp <- function(x, rate, q) {
  if (rate <= 0 || q >= 2 || q < 1) return(rep(0, length(x)))
  pdf_vals <- numeric(length(x))
  valid <- (x > 0)
  
  if (abs(q - 1) < 1e-5) {
    pdf_vals[valid] <- rate * exp(-rate * x[valid])
  } else {
    term <- 1 + (q - 1) * rate * x[valid]
    pdf_vals[valid] <- (2 - q) * rate * (term^(-1 / (q - 1)))
  }
  return(pdf_vals)
}


pqExp <- function(x, rate, q) {
  if (abs(q - 1) < 1e-5) {
    return(pexp(x, rate = rate))
  }
  cdf_vals <- numeric(length(x))
  valid <- (x > 0)
  
  
  term <- 1 + (q - 1) * rate * x[valid]
  cdf_vals[valid] <- 1 - term^(-(2 - q) / (q - 1))
  return(cdf_vals)
}

# (d) Classical Ishita
dIshita <- function(x, theta) {
  if (theta <= 0) return(rep(0, length(x)))
  pdf_vals <- numeric(length(x))
  valid <- (x > 0)
  pdf_vals[valid] <- (theta^3 / (theta^3 + 2)) * (theta + x[valid]^2) * exp(-theta * x[valid])
  return(pdf_vals)
}

pIshita <- function(x, theta) {
  sapply(x, function(val) {
    if (val <= 0) return(0)
    tryCatch(integrate(dIshita, lower = 0, upper = val, theta = theta)$value, error = function(e) NA)
  })
}


# 3. Negative Log-Likelihood Functions
nll_qDI <- function(par) {
  theta <- par[1]; q <- par[2]
  if (theta <= 1e-4 || q >= 1.33 || q <= 1.0) return(1e10)
  pdf_vals <- dqDI(x, theta, q)
  if (any(pdf_vals <= 0) || any(is.na(pdf_vals))) return(1e10)
  return(-sum(log(pdf_vals)))
}

nll_qDL <- function(par) {
  theta <- par[1]; q <- par[2]
  if (theta <= 1e-4 || q >= 1.49 || q <= 1.0) return(1e10)
  pdf_vals <- dqDL(x, theta, q)
  if (any(pdf_vals <= 0) || any(is.na(pdf_vals))) return(1e10)
  return(-sum(log(pdf_vals)))
}

nll_qExp <- function(par) {
  rate <- par[1]; q <- par[2]
  if (rate <= 1e-4 || q >= 1.99 || q <= 1.0) return(1e10)
  pdf_vals <- dqExp(x, rate, q)
  if (any(pdf_vals <= 0) || any(is.na(pdf_vals))) return(1e10)
  return(-sum(log(pdf_vals)))
}

nll_Ishita  <- function(par) -sum(log(dIshita(x, par[1])))
nll_Weibull <- function(par) -sum(dweibull(x, shape = par[1], scale = par[2], log = TRUE))
nll_Gamma   <- function(par) -sum(dgamma(x, shape = par[1], rate = par[2], log = TRUE))


# 4. Maximum Likelihood Estimation Optimization
fit_qDI     <- optim(c(0.5, 1.15), nll_qDI,     method = "L-BFGS-B", lower = c(0.01, 1.001), upper = c(5, 1.32))
fit_qDL     <- optim(c(0.5, 1.15), nll_qDL,     method = "L-BFGS-B", lower = c(0.01, 1.001), upper = c(5, 1.48))
fit_qExp    <- optim(c(0.3, 1.15), nll_qExp,    method = "L-BFGS-B", lower = c(0.01, 1.001), upper = c(5, 1.98))
fit_Ishita  <- optim(c(0.5),      nll_Ishita,  method = "L-BFGS-B", lower = c(0.01), upper = c(10))
fit_Weibull <- optim(c(1, 1),      nll_Weibull, method = "L-BFGS-B", lower = c(0.01, 0.01))
fit_Gamma   <- optim(c(1, 1),      nll_Gamma,   method = "L-BFGS-B", lower = c(0.01, 0.01))

# 5. Goodness-of-Fit Summary & Model Comparison
models <- list(
  "q-Deformed Ishita"  = list(fit = fit_qDI,     k = 2, pfun = function(v) pqDI(v, fit_qDI$par[1], fit_qDI$par[2])),
  "q-Deformed Lindley" = list(fit = fit_qDL,     k = 2, pfun = function(v) pqDL(v, fit_qDL$par[1], fit_qDL$par[2])),
  "q-Exponential"      = list(fit = fit_qExp,    k = 2, pfun = function(v) pqExp(v, fit_qExp$par[1], fit_qExp$par[2])),
  "Classical Ishita"   = list(fit = fit_Ishita,  k = 1, pfun = function(v) pIshita(v, fit_Ishita$par[1])),
  "Weibull"            = list(fit = fit_Weibull, k = 2, pfun = function(v) pweibull(v, fit_Weibull$par[1], fit_Weibull$par[2])),
  "Gamma"              = list(fit = fit_Gamma,   k = 2, pfun = function(v) pgamma(v, fit_Gamma$par[1], fit_Gamma$par[2]))
)

comparison_results <- data.frame(
  Model = character(),
  Params = numeric(),
  LogLik = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  KS_Stat = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(models)) {
  m <- models[[name]]
  logL <- -m$fit$value
  k <- m$k
  aic <- 2 * k - 2 * logL
  bic <- k * log(n) - 2 * logL
  
  ks_test <- ks.test(x, m$pfun)
  
  comparison_results[nrow(comparison_results) + 1, ] <- list(
    Model   = name,
    Params  = k,
    LogLik  = round(logL, 4),
    AIC     = round(aic, 4),
    BIC     = round(bic, 4),
    KS_Stat = round(ks_test$statistic, 4),
    p_value = round(ks_test$p.value, 4)
  )
}

comparison_results <- comparison_results[order(comparison_results$AIC), ]

print("--- MODEL COMPARISON TABLE ---")
print(comparison_results, row.names = FALSE)


# 6. Comparative Density and CDF Visualizations
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

x_grid <- seq(0.01, max(x) + 2, length.out = 500)

# Histogram & Density Overlay
hist(x, breaks = 25, probability = TRUE, col = "gray92", border = "gray60", main="", xlab = "Length (100 of miles)",
     ylab = "Density", ylim = c(0, 0.35))

lines(x_grid, dqDI(x_grid, fit_qDI$par[1], fit_qDI$par[2]), col = "firebrick", lwd = 2.5, lty = 1)
lines(x_grid, dqDL(x_grid, fit_qDL$par[1], fit_qDL$par[2]), col = "blue3", lwd = 2, lty = 2)
lines(x_grid, dqExp(x_grid, fit_qExp$par[1], fit_qExp$par[2]), col = "darkgreen", lwd = 2, lty = 3)
lines(x_grid, dIshita(x_grid, fit_Ishita$par[1]), col = "purple", lwd = 2, lty = 4)
lines(x_grid, dweibull(x_grid, fit_Weibull$par[1], fit_Weibull$par[2]), col = "darkorange", lwd = 2, lty = 5)

legend("topright", legend = c("q-Deformed Ishita", "q-Deformed Lindley", "q-Exponential", "Classical Ishita", "Weibull"),
       col = c("firebrick", "blue3", "darkgreen", "purple", "darkorange"),
       lty = 1:5, lwd = 2, bty = "n", cex = 0.8)

# Empirical vs Fitted CDFs
plot(ecdf(x), main = "Empirical CDF vs. Fitted CDFs", xlab = "Length (hundreds of miles)",
     ylab = "F(x)", col = "gray40", lwd = 1.5)

lines(x_grid, pqDI(x_grid, fit_qDI$par[1], fit_qDI$par[2]), col = "firebrick", lwd = 2.5, lty = 1)
lines(x_grid, pqDL(x_grid, fit_qDL$par[1], fit_qDL$par[2]), col = "blue3", lwd = 2, lty = 2)
lines(x_grid, pqExp(x_grid, fit_qExp$par[1], fit_qExp$par[2]), col = "darkgreen", lwd = 2, lty = 3)
lines(x_grid, pIshita(x_grid, fit_Ishita$par[1]), col = "purple", lwd = 2, lty = 4)
lines(x_grid, pweibull(x_grid, fit_Weibull$par[1], fit_Weibull$par[2]), col = "darkorange", lwd = 2, lty = 5)

legend("bottomright", legend = c("Emperical CDF", "q-Deformed Ishita", "q-Deformed Lindley", "q-Exponential", "Ishita", "Weibull"),
       col = c("gray40", "firebrick", "blue3", "darkgreen", "purple", "darkorange"),
       lty = c(1, 1:5), lwd = c(1.5, rep(2, 5)), bty = "n", cex = 0.75)

par(mfrow = c(1, 1))









# Box plot of River Length Data

data("rivers")

x <- rivers / 100

boxplot(
  x,
  main = "",
  ylab = "Length (hundreds of miles)",
  col = "red",
  border = "blue",
  outline = TRUE
)
