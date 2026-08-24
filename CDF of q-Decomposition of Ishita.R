rm(list = ls())
theta <- 0.5
q.values <- c(1.05, 1.20, 1.30, 1.32)
Cq <- function(theta, q){
  ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
}
pqDI <- function(x, theta, q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*(
    -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
    +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
    -(1+a*x)^(-(2-q)/(q-1))/(2-q)
    +1/(4-3*q)-2/(3-2*q)+1/(2-q)
  )
  
  F <- C*(term1+term2)
  
  return(F)
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,100,length=1000)

plot(x,
     pqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,1),
     xlab="x",
     ylab="F(x)",
     cex.lab=1.4,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  lines(x,
        pqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
}

legend("topleft",
       legend=c(
         expression(theta==0.5~","~q==1.05),
         expression(theta==0.5~","~q==1.20),
         expression(theta==0.5~","~q==1.30),
         expression(theta==0.5~","~q==1.32)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)



rm(list = ls())
theta <- 1.5
q.values <- c(1.05, 1.20, 1.30, 1.32)
Cq <- function(theta, q){
  ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
}
pqDI <- function(x, theta, q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*(
    -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
    +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
    -(1+a*x)^(-(2-q)/(q-1))/(2-q)
    +1/(4-3*q)-2/(3-2*q)+1/(2-q)
  )
  
  F <- C*(term1+term2)
  
  return(F)
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,5,length=1000)

plot(x,
     pqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,1),
     xlab="x",
     ylab="F(x)",
     cex.lab=1.4,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  lines(x,
        pqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
}

legend("topright",
       legend=c(
         expression(theta==1.5~","~q==1.05),
         expression(theta==1.5~","~q==1.20),
         expression(theta==1.5~","~q==1.30),
         expression(theta==1.5~","~q==1.32)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)





rm(list = ls())
theta <- 0.01
q.values <- c(0.75, 0.50, 0.25, 0.1)
Cq <- function(theta, q){
  ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
}
pqDI <- function(x, theta, q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*(
    -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
    +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
    -(1+a*x)^(-(2-q)/(q-1))/(2-q)
    +1/(4-3*q)-2/(3-2*q)+1/(2-q)
  )
  
  F <- C*(term1+term2)
  
  return(F)
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,100,length=1000)

plot(x,
     pqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,1),
     xlab="x",
     ylab="F(x)",
     cex.lab=1.4,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  lines(x,
        pqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
}

legend("topright",
       legend=c(
         expression(theta==0.01~","~q==0.75),
         expression(theta==0.01~","~q==0.50),
         expression(theta==0.01~","~q==0.25),
         expression(theta==0.01~","~q==0.1)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)



rm(list = ls())
theta <- 0.025
q.values <- c(0.75, 0.5, 0.25, 0.10)
Cq <- function(theta, q){
  ((2-q)*(3-2*q)*(4-3*q)*theta^3) /
    (theta^3*(3-2*q)*(4-3*q)+2)
}

pqDI <- function(x, theta, q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*(
    -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
    +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
    -(1+a*x)^(-(2-q)/(q-1))/(2-q)
    +1/(4-3*q)-2/(3-2*q)+1/(2-q)
  )
  
  F <- C*(term1+term2)
  
  return(F)
}

par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,50,length=1000)

plot(x,
     pqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,1),
     xlab="x",
     ylab="F(x)",
     cex.lab=1.4,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  lines(x,
        pqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
}

legend("topleft",
       legend=c(
         expression(theta==0.025~","~q==0.75),
         expression(theta==0.025~","~q==0.50),
         expression(theta==0.025~","~q==0.25),
         expression(theta==0.025~","~q==0.1)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)

