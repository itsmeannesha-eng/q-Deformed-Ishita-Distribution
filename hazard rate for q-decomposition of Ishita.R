rm(list = ls())
theta <- 0.5
q.values <- c(1.05,1.20,1.30,1.32)
Cq <- function(theta,q){
  
  ((2-q)*(3-2*q)*(4-3*q)*theta^3)/
    (theta^3*(3-2*q)*(4-3*q)+2)
  
}
dqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  C*(theta+x^2)*
    (1+(q-1)*theta*x)^(-1/(q-1))
  
}
pqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*
    (
      -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
      +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
      -(1+a*x)^(-(2-q)/(q-1))/(2-q)
      +1/(4-3*q)
      -2/(3-2*q)
      +1/(2-q)
    )
  
  C*(term1+term2)
  
}
hqDI <- function(x,theta,q){
  
  f <- dqDI(x,theta,q)
  
  F <- pqDI(x,theta,q)
  
  f/(1-F)
  
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,25,length=1200)

plot(x,
     hqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,max(sapply(q.values,function(q)
       hqDI(x,theta,q)),na.rm=TRUE)),
     xlab="x",
     ylab="h(x)",
     cex.lab=1.3,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  
  lines(x,
        hqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
  
}

legend("topright",
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
q.values <- c(1.05,1.20,1.30,1.32)
Cq <- function(theta,q){
  
  ((2-q)*(3-2*q)*(4-3*q)*theta^3)/
    (theta^3*(3-2*q)*(4-3*q)+2)
  
}
dqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  C*(theta+x^2)*
    (1+(q-1)*theta*x)^(-1/(q-1))
  
}
pqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*
    (
      -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
      +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
      -(1+a*x)^(-(2-q)/(q-1))/(2-q)
      +1/(4-3*q)
      -2/(3-2*q)
      +1/(2-q)
    )
  
  C*(term1+term2)
  
}
hqDI <- function(x,theta,q){
  
  f <- dqDI(x,theta,q)
  
  F <- pqDI(x,theta,q)
  
  f/(1-F)
  
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,25,length=1200)

plot(x,
     hqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,max(sapply(q.values,function(q)
       hqDI(x,theta,q)),na.rm=TRUE)),
     xlab="x",
     ylab="h(x)",
     cex.lab=1.3,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  
  lines(x,
        hqDI(x,theta,q.values[i]),
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
q.values <- c(0.75,0.50,0.25,0.10)
Cq <- function(theta,q){
  
  ((2-q)*(3-2*q)*(4-3*q)*theta^3)/
    (theta^3*(3-2*q)*(4-3*q)+2)
  
}
dqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  C*(theta+x^2)*
    (1+(q-1)*theta*x)^(-1/(q-1))
  
}
pqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*
    (
      -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
      +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
      -(1+a*x)^(-(2-q)/(q-1))/(2-q)
      +1/(4-3*q)
      -2/(3-2*q)
      +1/(2-q)
    )
  
  C*(term1+term2)
  
}
hqDI <- function(x,theta,q){
  
  f <- dqDI(x,theta,q)
  
  F <- pqDI(x,theta,q)
  
  f/(1-F)
  
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,25,length=1200)

plot(x,
     hqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,max(sapply(q.values,function(q)
       hqDI(x,theta,q)),na.rm=TRUE)),
     xlab="x",
     ylab="h(x)",
     cex.lab=1.3,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  
  lines(x,
        hqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
  
}

legend("topleft",
       legend=c(
         expression(theta==0.01~","~q==0.75),
         expression(theta==0.01~","~q==0.50),
         expression(theta==0.01~","~q==0.25),
         expression(theta==0.01~","~q==0.10)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)




rm(list = ls())

theta <- 0.025
q.values <- c(0.75,0.50,0.25,0.10)
Cq <- function(theta,q){
  
  ((2-q)*(3-2*q)*(4-3*q)*theta^3)/
    (theta^3*(3-2*q)*(4-3*q)+2)
  
}

dqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  C*(theta+x^2)*
    (1+(q-1)*theta*x)^(-1/(q-1))
  
}

pqDI <- function(x,theta,q){
  
  C <- Cq(theta,q)
  
  a <- (q-1)*theta
  
  term1 <- (1-(1+a*x)^(-(2-q)/(q-1)))/(2-q)
  
  term2 <- (1/theta^3)*
    (
      -(1+a*x)^(-(4-3*q)/(q-1))/(4-3*q)
      +2*(1+a*x)^(-(3-2*q)/(q-1))/(3-2*q)
      -(1+a*x)^(-(2-q)/(q-1))/(2-q)
      +1/(4-3*q)
      -2/(3-2*q)
      +1/(2-q)
    )
  
  C*(term1+term2)
  
}

hqDI <- function(x,theta,q){
  
  f <- dqDI(x,theta,q)
  
  F <- pqDI(x,theta,q)
  
  f/(1-F)
  
}
par(
  family="serif",
  las=1,
  tck=-0.02
)

cols <- c("black","black","brown3","darkgreen")
ltys <- c(1,2,3,4)

x <- seq(0,25,length=1200)

plot(x,
     hqDI(x,theta,q.values[1]),
     type="l",
     lwd=2,
     col=cols[1],
     lty=ltys[1],
     ylim=c(0,max(sapply(q.values,function(q)
       hqDI(x,theta,q)),na.rm=TRUE)),
     xlab="x",
     ylab="h(x)",
     cex.lab=1.3,
     cex.axis=1.2)

for(i in 2:length(q.values)){
  
  lines(x,
        hqDI(x,theta,q.values[i]),
        col=cols[i],
        lty=ltys[i],
        lwd=2)
  
}

legend("topleft",
       legend=c(
         expression(theta==0.025~","~q==0.75),
         expression(theta==0.025~","~q==0.50),
         expression(theta==0.025~","~q==0.25),
         expression(theta==0.025~","~q==0.10)
       ),
       col=cols,
       lty=ltys,
       lwd=2,
       bty="n",
       cex=0.9)

