rm(list = ls())

# Parameters
theta <- 0.25
q.values <- c(0.75, 0.50, 0.25, 0.10)

# q-Deformed Ishita PDF (q < 1)
dqDI <- function(x, theta, q){
  
  # Normalizing constant
  C <- ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
  
  # Upper support
  xmax <- 1/((1-q)*theta)
  
  pdf <- numeric(length(x))
  
  ind <- (x >= 0 & x <= xmax)
  
  pdf[ind] <- C*(theta + x[ind]^2)*
    (1-(1-q)*theta*x[ind])^(1/(1-q))
  
  return(pdf)
}

# x range
x <- seq(0,25,length=1500)

# Plot settings
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("goldenrod","black","brown3","darkgreen")
ltys <- c(2,1,2,4)

# Plot first curve
plot(x,
     dqDI(x,theta,q.values[1]),
     type="l",
     col=cols[1],
     lty=ltys[1],
     lwd=2,
     ylim=c(0,max(sapply(q.values,function(q)
       dqDI(x,theta,q)))),
     xlab="x",
     ylab="Density",
     cex.lab=1.3,
     cex.axis=1.2)

# Remaining curves
for(i in 2:length(q.values)){
  lines(x,
        dqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
}

# Legend
legend("topright",
       legend=c(
         expression(theta==0.25~","~q==0.75),
         expression(theta==0.25~","~q==0.50),
         expression(theta==0.25~","~q==0.25),
         expression(theta==0.25~","~q==0.10)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)



rm(list = ls())

# Parameters
theta <- 1.5
q.values <- c(0.95, 0.80, 0.70, 0.60)

# q-Deformed Ishita PDF
dqDI <- function(x, theta, q){
  
  # Normalizing constant
  C <- ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
  
  # Upper support
  xmax <- 1/((1-q)*theta)
  
  pdf <- numeric(length(x))
  
  ind <- (x >= 0 & x <= xmax)
  
  pdf[ind] <- C*(theta + x[ind]^2)*
    (1-(1-q)*theta*x[ind])^(1/(1-q))
  
  return(pdf)
}

# x values
x <- seq(0, 3, length = 1000)

# Plot style
par(
  family = "serif",
  las = 1,
  tck = -0.02
)

cols <- c("goldenrod", "black", "brown3", "darkgreen")
ltys <- c(2, 1, 2, 4)

# Plot first curve
plot(x,
     dqDI(x, theta, q.values[1]),
     type = "l",
     col = cols[1],
     lty = ltys[1],
     lwd = 2,
     ylim = c(0,
              max(sapply(q.values,
                         function(q) dqDI(x, theta, q)),
                  na.rm = TRUE)),
     xlab = "x",
     ylab = "Density",
     cex.lab = 1.3,
     cex.axis = 1.2)

# Remaining curves
for(i in 2:length(q.values)){
  lines(x,
        dqDI(x, theta, q.values[i]),
        col = cols[i],
        lty = ltys[i],
        lwd = 2)
}

# Legend
legend("topright",
       legend = c(
         expression(theta==1.5~","~q==0.95),
         expression(theta==1.5~","~q==0.80),
         expression(theta==1.5~","~q==0.70),
         expression(theta==1.5~","~q==0.60)
       ),
       col = cols,
       lty = ltys,
       lwd = 2,
       bty = "n",
       cex = 0.9)

rm(list = ls())

# Parameters
theta <- 0.01
q.values <- c(0.90, 0.50, 0.30, 0.20)

# q-Deformed Ishita PDF
dqDI <- function(x, theta, q){
  
  # Normalizing constant
  C <- ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
  
  # Upper support
  xmax <- 1/((1-q)*theta)
  
  pdf <- numeric(length(x))
  
  ind <- (x >= 0 & x <= xmax)
  
  pdf[ind] <- C*(theta + x[ind]^2)*
    (1-(1-q)*theta*x[ind])^(1/(1-q))
  
  return(pdf)
}

# x values
x <- seq(0, 50, length = 1000)

# Plot style
par(
  family = "serif",
  las = 1,
  tck = -0.02
)

cols <- c("goldenrod", "black", "brown3", "darkgreen")
ltys <- c(2, 1, 2, 4)

# Plot first curve
plot(x,
     dqDI(x, theta, q.values[1]),
     type = "l",
     col = cols[1],
     lty = ltys[1],
     lwd = 2,
     ylim = c(0,
              max(sapply(q.values,
                         function(q) dqDI(x, theta, q)),
                  na.rm = TRUE)),
     xlab = "x",
     ylab = "Density",
     cex.lab = 1.3,
     cex.axis = 1.2)

# Remaining curves
for(i in 2:length(q.values)){
  lines(x,
        dqDI(x, theta, q.values[i]),
        col = cols[i],
        lty = ltys[i],
        lwd = 2)
}

# Legend
legend("topleft",
       legend = c(
         expression(theta==0.01~","~q==0.90),
         expression(theta==0.01~","~q==0.50),
         expression(theta==0.01~","~q==0.30),
         expression(theta==0.01~","~q==0.20)
       ),
       col = cols,
       lty = ltys,
       lwd = 2,
       bty = "n",
       cex = 0.9)

