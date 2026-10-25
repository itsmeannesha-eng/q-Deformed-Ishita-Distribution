
rm(list = ls())

Cq <- function(theta, q){
  ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
}

pqDI <- function(x, theta, q){
  C <- Cq(theta, q)
  a <- (q - 1) * theta
  
  term1 <- (1 - (1 + a * x)^(-(2 - q) / (q - 1))) / (2 - q)
  
  term2 <- (1 / (theta^3 * (q - 1)^2)) * (
    -(1 + a * x)^(-(4 - 3 * q) / (q - 1)) / (4 - 3 * q)
    + 2 * (1 + a * x)^(-(3 - 2 * q) / (q - 1)) / (3 - 2 * q)
    - (1 + a * x)^(-(2 - q) / (q - 1)) / (2 - q)
    + 1 / (4 - 3 * q) - 2 / (3 - 2 * q) + 1 / (2 - q)
  )
  
  F <- C * (term1 + term2)
  
  # Handle bounded support threshold safely if q < 1
  if(q < 1){
    xmax <- 1 / ((1 - q) * theta)
    F[x > xmax] <- 1
  }
  return(F)
}

# Setup 2x2 multi-panel layout, fonts, and margins
par(mfrow = c(2, 2), family = "serif", las = 1, tck = -0.02, mar = c(4, 4, 2, 1))

cols <- c("black", "blue", "brown3", "darkgreen")
ltys <- c(1, 2, 3, 4)


# ==========================================
# Panel 1: Theta = 0.5, q > 1 (Heavy-tailed)
# ==========================================
theta1 <- 0.5
q.values1 <- c(1.05, 1.20, 1.30, 1.32)
x1 <- seq(0, 40, length = 1000)

plot(x1, pqDI(x1, theta1, q.values1[1]),
     type = "l", lwd = 2, col = cols[1], lty = ltys[1],
     ylim = c(0, 1), xlab = "x", ylab = "F(x)",
     cex.lab = 1.2, cex.axis = 1.0)

for(i in 2:length(q.values1)){
  lines(x1, pqDI(x1, theta1, q.values1[i]), col = cols[i], lty = ltys[i], lwd = 2)
}

legend("bottomright",
       legend = c(
         expression(theta==0.5~","~q==1.05),
         expression(theta==0.5~","~q==1.20),
         expression(theta==0.5~","~q==1.30),
         expression(theta==0.5~","~q==1.32)
       ),
       col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.75)


# ==========================================
# Panel 2: Theta = 1.5, q > 1 (Heavy-tailed)
# ==========================================
theta2 <- 1.5
q.values2 <- c(1.05, 1.20, 1.30, 1.32)
x2 <- seq(0, 10, length = 1000)

plot(x2, pqDI(x2, theta2, q.values2[1]),
     type = "l", lwd = 2, col = cols[1], lty = ltys[1],
     ylim = c(0, 1), xlab = "x", ylab = "F(x)",
     cex.lab = 1.2, cex.axis = 1.0)

for(i in 2:length(q.values2)){
  lines(x2, pqDI(x2, theta2, q.values2[i]), col = cols[i], lty = ltys[i], lwd = 2)
}

legend("bottomright",
       legend = c(
         expression(theta==1.5~","~q==1.05),
         expression(theta==1.5~","~q==1.20),
         expression(theta==1.5~","~q==1.30),
         expression(theta==1.5~","~q==1.32)
       ),
       col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.75)


# ==========================================
# Panel 3: Theta = 0.01, q < 1 (Bounded support)
# ==========================================
theta3 <- 0.01
q.values3 <- c(0.75, 0.50, 0.25, 0.1)
x3 <- seq(0, 120, length = 1000)

plot(x3, pqDI(x3, theta3, q.values3[1]),
     type = "l", lwd = 2, col = cols[1], lty = ltys[1],
     ylim = c(0, 1), xlab = "x", ylab = "F(x)",
     cex.lab = 1.2, cex.axis = 1.0)

for(i in 2:length(q.values3)){
  lines(x3, pqDI(x3, theta3, q.values3[i]), col = cols[i], lty = ltys[i], lwd = 2)
}

legend("topleft",
       legend = c(
         expression(theta==0.01~","~q==0.75),
         expression(theta==0.01~","~q==0.50),
         expression(theta==0.01~","~q==0.25),
         expression(theta==0.01~","~q==0.1)
       ),
       col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.75)


# ==========================================
# Panel 4: Theta = 0.025, q < 1 (Bounded support)
# ==========================================
theta4 <- 0.025
q.values4 <- c(0.75, 0.5, 0.25, 0.10)
x4 <- seq(0, 60, length = 1000)

plot(x4, pqDI(x4, theta4, q.values4[1]),
     type = "l", lwd = 2, col = cols[1], lty = ltys[1],
     ylim = c(0, 1), xlab = "x", ylab = "F(x)",
     cex.lab = 1.2, cex.axis = 1.0)

for(i in 2:length(q.values4)){
  lines(x4, pqDI(x4, theta4, q.values4[i]), col = cols[i], lty = ltys[i], lwd = 2)
}

legend("topleft",
       legend = c(
         expression(theta==0.025~","~q==0.75),
         expression(theta==0.025~","~q==0.50),
         expression(theta==0.025~","~q==0.25),
         expression(theta==0.025~","~q==0.1)
       ),
       col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.75)

# Reset layout parameters back to single panel default
par(mfrow = c(1, 1))

