## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# data preparation
set.seed(111)
n <- 200
Num.Cmp <- 8
pro <- rep(1/8, Num.Cmp)
multi <- sample(1:Num.Cmp, n, replace = T, prob=pro)
mu <- 3 * ((2/3)^(1:Num.Cmp) - 1)
sigma <- (2/3)^(1:Num.Cmp)
x <- NULL
for (ii in 1:Num.Cmp) {
  com_txt <- paste("com", ii, " <- rnorm(length(which(multi==", ii, ")), mean=", mu[ii], ", sd=", sigma[ii], ")",sep="")
  eval(parse(text=com_txt))
  com_txt <- paste("x <- c(x, com", ii, ")", sep="")
  eval(parse(text=com_txt))
}
# true density function, y is h, and z is v.
y <- seq(-3, 1, 0.01)
z <- rep(0, length(y))
for (ii in 1:Num.Cmp) {
  z <- z + pro[ii] * dnorm(y, mean=mu[ii], sd=sigma[ii])
}

## -----------------------------------------------------------------------------
# Rodeo local for 1-dim
rodeo.local.bw1 <- function(xx, x, h.init=1.3/log(log(n)), beta=0.9, cn=log(n)/n){
  # bandwidth selection
  # para: xx        target point at which we want to estimate f(x)
  # para: x         samples
  # para: beta      learning rate
  # value: h.init   bandwidth
  h1 <- h.init
  while(TRUE){
    Z.i <- ((xx - x)^2 - h1^2) * exp(- (xx - x)^2 / h1^2 / 2) / h1^4 / sqrt(2*pi)
    Z <- mean(Z.i) 
    s <- var(Z.i)
    lam <- sqrt(2*s*log(n*cn)) / 10
    if (abs(Z) > lam) {
      h1 <- h1 * beta
    } else {
      return(h1)
    }
  }
}
rodeo.local1 <- function(t, x){
  # estimate target points t
  # para: t   target points, a vector
  # para: x   data points, a vector
  h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
  K <- stats::dnorm
  f.hat <- unlist(base::lapply(X=1:length(t), FUN=function(ii) mean(K((t[ii] - x) / h[ii]))/ h[ii]))
  return(f.hat)
}
# plot h and density
t <- seq(-3, 1, 0.05)
h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
plot(t, h, "l", main="Bandwidth of Rodeo", xlab="x", ylab="h")
fit.rodeo <- rodeo.local1(t=t, x=x)
plot(t, fit.rodeo, "l", lty=2, ylim=c(0,2.5), xlim=c(-3, 1), main="Rodeo", xlab="x", ylab="Density")
lines(y, z, lty=1, lwd=2)
legend("topright", legend=c("True Density", "Rodeo"), col=1, lty=c(1, 2), lwd=c(2, 1))

## -----------------------------------------------------------------------------
# Required packages
library(kedd)

# Traditional density estimation method without using boundary correction techniques
set.seed(123)
dat <- rexp(1000, 1)
plot(density(dat, kernel = "epanechnikov"), ylim = c(0, 1.2), main = "", xlab = "x")
lines(seq(0, 8, by = 0.02), dexp(seq(0, 8, by = 0.02), 1), col = "blue")
abline(v = 0, col = "black", lty = 2)
legend("topright", c("true density", "estimated density"), lty = c(1,1), col = c("blue","black"))

## -----------------------------------------------------------------------------
# Generlized Jacknifing method
kernel.new <- function(x, u, h) {
  # Compute a0, a1, a2
  lb <- -1
  ub <- pmin(u / h, 1)
  a0 <- 0.75 * (ub - lb) - 0.25 * (ub^3 - lb^3)
  a1 <- 3 / 8 * (ub^2 - lb^2) - 3 / 16 * (ub^4 - lb^4)
  a2 <- 0.25 * (ub^3 - lb^3) - 0.15 * (ub^5 - lb^5)
  ((3/4)*(1-x^2)*(abs(x) <= 1)) * (a2 - a1*x) / (a0*a2 - a1^2)
}
den.est <- function(u, ui, h) {
sapply(u, function(u) ifelse(u < 0, 0, mean(kernel.new((u - ui) / h, u, h)) / h))
}
# Estimated value at boundary 0
x <- seq(0, 8, by = 0.02)
y <- den.est(x, dat, 2 * bw.bcv(dat))
y0 <- den.est(0, dat, 2 * bw.bcv(dat))
y0_True <- dexp(0,1)
cat("The estimated value at boundary 0 by using GJack is ", y0, "\nThe true value at boundary 0 is ", y0_True)

# Plot
plot(x, y, type = "l", ylim = c(0, 1.2), ylab = "Density")
lines(x, dexp(x, 1), col = "blue")
abline(v = 0, col = "black", lty = 2)
legend("topright", c("true density", "estimated density"), lty = c(1,1), col = c("blue","black"))

