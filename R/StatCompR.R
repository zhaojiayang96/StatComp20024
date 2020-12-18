#' @title A kernel density estimation bw selection method
#' @description A KDE bw selection method for one dimension
#' @importFrom stats var
#' @param xx the point to be estimated
#' @param x data points
#' @param h.init initial smooth parameter
#' @param beta learning rate, default 0.9
#' @param cn a turning parameter, default log(n)/n
#' @param n the number of trials.
#' @return h selected by this method
#' @examples
#' \dontrun{
#' set.seed(111)
#' n <- 200
#' Num.Cmp <- 8
#' pro <- rep(1/8, Num.Cmp)
#' multi <- sample(1:Num.Cmp, n, replace = T, prob=pro)
#' mu <- 3 * ((2/3)^(1:Num.Cmp) - 1)
#' sigma <- (2/3)^(1:Num.Cmp)
#' x <- NULL
#' for (ii in 1:Num.Cmp) {
#'   com_txt <- paste("com", ii, " <- rnorm(length(which(multi==", ii, ")), 
#'                mean=", mu[ii], ", sd=", sigma[ii], ")",sep="")
#'   eval(parse(text=com_txt))
#'   com_txt <- paste("x <- c(x, com", ii, ")", sep="")
#'   eval(parse(text=com_txt))
#' }
#' 
#' # true density function, y is h, and z is v.
#' y <- seq(-3, 1, 0.01)
#' z <- rep(0, length(y))
#' for (ii in 1:Num.Cmp) {
#'   z <- z + pro[ii] * dnorm(y, mean=mu[ii], sd=sigma[ii])
#' }
#' t <- seq(-3, 1, 0.05)
#' h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
#' plot(t, h, "l", main="Bandwidth of Rodeo", xlab="x", ylab="h")
#' }
#' @export
rodeo.local.bw1 <- function(xx, x, h.init=1.3/log(log(n)), beta=0.9, cn=log(n)/n){
  n <-200
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


#' @title A kernel density estimation bw selection method
#' @description A KDE bw selection method for one dimension
#' @param t target points, a vector
#' @param x data points, vector
#' @return a vector of estimated f(x)
#' @examples
#' \dontrun{
#' set.seed(111)
#' n <- 200
#' Num.Cmp <- 8
#' pro <- rep(1/8, Num.Cmp)
#' multi <- sample(1:Num.Cmp, n, replace = T, prob=pro)
#' mu <- 3 * ((2/3)^(1:Num.Cmp) - 1)
#' sigma <- (2/3)^(1:Num.Cmp)
#' x <- NULL
#' for (ii in 1:Num.Cmp) {
#'   com_txt <- paste("com", ii, " <- rnorm(length(which(multi==", ii, ")), 
#'                        mean=", mu[ii], ", sd=", sigma[ii], ")",sep="")
#'   eval(parse(text=com_txt))
#'   com_txt <- paste("x <- c(x, com", ii, ")", sep="")
#'   eval(parse(text=com_txt))
#' }
#' 
#' # true density function, y is h, and z is v.
#' y <- seq(-3, 1, 0.01)
#' z <- rep(0, length(y))
#' for (ii in 1:Num.Cmp) {
#'   z <- z + pro[ii] * dnorm(y, mean=mu[ii], sd=sigma[ii])
#' }
#' t <- seq(-3, 1, 0.05)
#' h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
#' t <- seq(-3, 1, 0.05)
#' h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
#' fit.rodeo <- rodeo.local1(t=t, x=x)
#' plot(t, fit.rodeo, "l", lty=2, ylim=c(0,2.5), xlim=c(-3, 1), main="Rodeo", xlab="x", ylab="Density")
#' lines(y, z, lty=1, lwd=2)
#' legend("topright", legend=c("True Density", "Rodeo"), col=1, lty=c(1, 2), lwd=c(2, 1))
#' }
#' @export
rodeo.local1 <- function(t, x){
  # estimate target points t
  # para: t   target points, a vector
  # para: x   data points, a vector
  h <- unlist(base::lapply(X=t, FUN=rodeo.local.bw1, x=x))
  K <- stats::dnorm
  f.hat <- unlist(base::lapply(X=1:length(t), FUN=function(ii) mean(K((t[ii] - x) / h[ii]))/ h[ii]))
  return(f.hat)
}


#' @title linear multiple kernel function for Generlized Jacknifing method of KDE
#' @description This new kernel function is special kernel function form of Generlized Jacknifing method
#' @param x the input
#' @param u the input
#' @param h the bandwidth (smoothing parameter)
#' @return function value(s) of the kernel function at certain point(s)
#' @examples
#' \dontrun{
#' set.seed(123)
#' library(stats)
#' x <- seq(0, 8, by = 0.02)
#' dat <- rexp(1000, 1)
#' h <- 2 * bw.bcv(dat)
#' kernel.new.values <- kernel.new((x-dat)/h, x, h)
#' }
#' @export
kernel.new <- function(x, u, h) {
  # Compute a0, a1, a2
  lb <- -1
  ub <- pmin(u / h, 1)
  a0 <- 0.75 * (ub - lb) - 0.25 * (ub^3 - lb^3)
  a1 <- 3 / 8 * (ub^2 - lb^2) - 3 / 16 * (ub^4 - lb^4)
  a2 <- 0.25 * (ub^3 - lb^3) - 0.15 * (ub^5 - lb^5)
  ((3/4)*(1-x^2)*(abs(x) <= 1)) * (a2 - a1*x) / (a0*a2 - a1^2)
}

#' @title Kernel density estimation results at boundary
#' @description Kernel density estimation results using linear multiple kernel of Generlized Jacknifing method
#' @param u the point to be estimated
#' @param ui input samples of the target density function
#' @param h the bandwidth (smoothing parameter)
#' @return a vector of estimated density function values at the points to be estimated
#' @examples
#' \dontrun{
#' x <- seq(0, 8, by = 0.02)
#' y <- den.est(x, dat, 2 * bw.bcv(dat))
#' y0 <- den.est(0, dat, 2 * bw.bcv(dat))
#' cat("The estimated value at boundary 0 by using GJack is ", y0)
#' }
#' @export
den.est <- function(u, ui, h) {
  sapply(u, function(u) ifelse(u < 0, 0, mean(kernel.new((u - ui) / h, u, h)) / h))
}