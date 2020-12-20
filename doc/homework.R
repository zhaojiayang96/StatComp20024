## -----------------------------------------------------------------------------
n<-1000
u<-runif(n)
x<-2/(sqrt(1-u))
hist(x,prob=TRUE,breaks=80,main=expression(f(x)==8/x^3),xlim=c(2,30))
y<-seq(0,30,.01)
lines(y,8/(y^3))

## -----------------------------------------------------------------------------
u1<-runif(1000,-1,1)
u2<-runif(1000,-1,1)
u3<-runif(1000,-1,1)
i<-1
z<-rep(0,1000)
while(i<=1000)
{
  if(abs(u3[i])>abs(u2[i]) &&  abs(u3[i])>abs(u1[i]))
     {z[i]<-u2[i]} 
  else
     {z[i]<-u3[i]} 
i=i+1
}
hist(z,prob=TRUE,main=expression(f(x)==3/4*(1-x^2)))

## -----------------------------------------------------------------------------
n <- 1000
r <- 4
beta <- 2
gammarv <- rgamma(n, shape = r, rate = beta)
x <- rexp(n, rate = gammarv)
hist(x[x<5], freq = FALSE, breaks = seq(0,5,0.1), main = "Histogram of the Pareto sample", xlab = "value")  
f <- function(x) {64/(2+x)^5}    # pdf of Pareto distribution
curve(f, 0, 5, col = 2, add = TRUE)    
legend(3, 1, "true density", col = 2, lwd = 1)    # add a legend

## -----------------------------------------------------------------------------
set.seed(1)
m <- 1e5
x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi / 3

## -----------------------------------------------------------------------------
c(theta.hat, 1/2)

## -----------------------------------------------------------------------------
set.seed(2)
m <- 1e6
x1 <- runif(m, min=0, max=1)
theta_hat_1 <- exp(x1)

## -----------------------------------------------------------------------------
set.seed(3)
m <- 5e5
x2 <- runif(m, min=0, max=1);
theta_hat_2 <- (exp(x2) + exp(1-x2)) / 2

## -----------------------------------------------------------------------------
c(1-var(theta_hat_2) / var(theta_hat_1))

## -----------------------------------------------------------------------------
m<-1e4
g<-function(x){
  exp(-x^2/2)*x^2/sqrt(2*pi)*(x>1)
  }
  

#####using f1######
x<-rweibull(m,2,sqrt(2))
gf<-g(x)/dweibull(x,2,sqrt(2))
theta1<-mean(gf)
se1<-sd(gf)



####using f2######
x<-rgamma(m,3,2)
gf<-g(x)/dgamma(x,3,2)
theta2<-mean(gf)
se2<-sd(gf)

theta.hat<-c(theta1,theta2)

se<-c(se1,se2)
rbind(theta.hat,se)


## -----------------------------------------------------------------------------
t<-seq(1,10,0.1)
g<-exp(-t^2/2)*t^2/sqrt(2*pi)
f1<-dweibull(t,2,sqrt(2))
f2<-dgamma(t,3,2)
 
####get the curve of g,f1 and f2####
plot(t,g,type="l",col="black",main="compare g(t), f1(t) and f2(t) ")   
lines(t,f1,col="red")  
lines(t,f2,col="green")  
legend("topright",legend =c('g(t)','f1(t)',"f2(t)") ,lty=1,col=c("black","red","green")) 


## -----------------------------------------------------------------------------
t<-seq(1,10,0.1)
 g<-exp(-t^2/2)*t^2/sqrt(2*pi)
 f1<-dweibull(t,2,sqrt(2))
 f2<-dgamma(t,3,2)
r1<-g/f1
r2<-g/f2
plot(t,r1,col="red", type = "l")
lines(t,r2,col="green")
title(main="ratio function")

## -----------------------------------------------------------------------------
M <- 10000
N <- 50 
k <- 5 
r <- M/k 
T5 <- numeric(k)
est <- matrix(0, N, 2)
#use reverse transform method
g<-function(x,a,b) exp(-x)/(1+x^2)*(x>a)*(x<b)
h<-function(u,a,b) -log(exp(-a)-u*(exp(-a)-exp(-b)))
fg<-function(x,a,b) g(x,a,b)/(exp(-x)/(exp(-a)-exp(-b)))
for (i in 1:N) {
u<-runif(M)
u.s<-runif(M/k)
#importance sampling
est[i, 1] <- mean(fg(h(u,0,1),0,1))
#stratified importance sampling
for(j in 1:k) T5[j]<-mean(fg(h(u.s,(j-1)/k,j/k),(j-1)/k,j/k))
est[i, 2] <- sum(T5)
}
#apply(est,2,mean)
#apply(est,2,sd)
cat("theta.hat: ",apply(est,2,mean),"\n","se:",apply(est,2,sd))

## ---- cache=TRUE, fig.height=4, fig.width=5.5---------------------------------
mu<-1;
sigma<-1
n<-20
alpha<-0.05
UCL<-replicate(1000,expr = {
  x<-rlnorm(n,mu,sigma)
  y<-log(x)
  abs(sqrt(n/var(y))*(mean(y)-mu))
})
mean(UCL<qt(1-alpha/2,n-1))
cat("The empirical covering probability is", mean(UCL<qt(1-alpha/2,n-1)))

## -----------------------------------------------------------------------------
set.seed(12345)

n = 20
alpha = 0.05
m = 1000
UCLvar1 = UCLmean1 = UCLvar2 = UCLmean2 = numeric(m)
UCLvar = UCLmean = numeric(2)

UCLvar1 = replicate(1000, expr = {
x = rnorm(n, mean = 0, sd = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1) })
UCLvar[1] = mean(UCLvar1 > 4)

UCLvar2 = replicate(1000, expr = {
x = rchisq(n, df = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
UCLvar[2] = mean(UCLvar2 > 4)

UCLmean1 = replicate(1000,expr={
  y = rnorm(n, mean = 0, sd = 2)
  (mean(y)-sd(y)*qt(df=n-1,alpha)/sqrt(n))
  
})
UCLmean[1] = mean(UCLmean1 > 0)

UCLmean2 = replicate(1000,expr={
  y = rchisq(n,df=2)
  (mean(y)-sd(y)*qt(df=n-1,alpha)/sqrt(n))
  
})
UCLmean[2] = mean(UCLmean2 > 2)

f = data.frame(UCLmean,UCLvar,row.names = c("normal distribution","chi-square distribution"))
knitr::kable(f)

## ----beta---------------------------------------------------------------------

set.seed(12345)

sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

# beta(a,a)
pwr_beta = function(a){
 alpha = 0.1
 n = 20
 m = 1e4
 N = length(a)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rbeta(n, a[j], a[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr * (1-pwr) / m) 
 return(list(pwr = pwr,se = se))
}

 a = c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
 pwr = pwr_beta(a)$pwr
 # plot the power
 se = pwr_beta(a)$se
 plot(a, pwr, type = "b", xlab = "a", ylab = "pwr", pch=16)
 abline(h = 0.1, lty = 2)
 lines(a, pwr+se, lty = 4)
 lines(a, pwr-se, lty = 4)

## ----t------------------------------------------------------------------------

# t(v)
pwr_t = function(v){
 
 alpha = 0.1
 n = 20
 m = 1e3
 N = length(v)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rt(n,v[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr*(1-pwr) / m) 
  return(list(pwr = pwr,se = se))
}

v = seq(1,20)
pwr = pwr_t(v)$pwr
se = pwr_t(v)$se
# plot the power
plot(v, pwr, type = "b", xlab = "v", ylab = "pwr", ylim = c(0,1),pch=16)
abline(h = 0.1, lty = 2)
lines(v, pwr+se, lty = 4)
lines(v, pwr-se, lty = 4)


## -----------------------------------------------------------------------------
count5test <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        return(as.integer(max(c(outx, outy)) > 5))
}
set.seed(12345)
alpha.hat <- 0.055
n <- c(10, 20, 50, 100, 500, 1000)
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
m <- 1e4
result <- matrix(0, length(n), 2)
for (i in 1:length(n)){
  ni <- n[i]
  tests <- replicate(m, expr={
    x <- rnorm(ni, mu1, sigma1)
    y <- rnorm(ni, mu2, sigma2)
    Fp <- var.test(x, y)$p.value
    Ftest <- as.integer(Fp <= alpha.hat)
    c(count5test(x, y), Ftest)
    })
  result[i, ] <- rowMeans(tests)
}
data.frame(n=n, C5=result[, 1], Fp=result[, 2])


## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library("bootstrap")
n <- nrow(scor)
theta_hat <- cor(law$LSAT, law$GPA)
theta_jack <- numeric(n)
for(i in 1:n){
  x <- scor[-i,]
  theta_jack[i] <- cor(law$LSAT[-i],law$GPA[-i])
}
bias_jack <- (n-1)*(mean(theta_jack)-theta_hat)
se_jack <- sqrt((n-1)*mean((theta_jack-mean(theta_jack))^2))
print(round(c(bias_jack=bias_jack,se_jack=se_jack),4))

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
aircondit <- as.matrix(aircondit)
boot.mean <- function(x,i) mean(x[i])
boot.obj <- boot(aircondit, statistic=boot.mean, R=2000)
print(boot.ci(boot.obj, type = c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------

library(bootstrap)
set.seed(12345)

n = nrow(scor)
lambda_hat = eigen(cov(scor))$values
theta_hat = lambda_hat[1] / sum(lambda_hat)
theta_j = rep(0,n)

for (i in 1:n) {

x = scor [-i,]
lambda = eigen(cov(x))$values
theta_j[i] = lambda[1]/sum(lambda)

}
#estimated bias of theta_hat
bias_jack = (n-1)*(mean(theta_j)-theta_hat)
#estimated se of theta_hat
se_jack = (n-1)*sqrt(var(theta_j)/n)

print(round(c(bias_jack=bias_jack,se_jack=se_jack),4))


## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- numeric(n*(n-1)/2)
e2 <- numeric(n*(n-1)/2)
e3 <- numeric(n*(n-1)/2)
e4 <- numeric(n*(n-1)/2)
count <- 0
for (i in 1:(n-1))
  for (j in (i+1):n) {
    count <- count+1
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    P1 <- lm(y~x)
    y1_1 <- chemical[i]*P1$coef[2] + P1$coef[1]
    y1_2 <- chemical[j]*P1$coef[2] + P1$coef[1]
    e1[count] <- (magnetic[i]-y1_1)^2+(magnetic[j]-y1_2)^2
    
    P2 <- lm(y~x+I(x^2))
    y2_1 <- P2$coef[1] + P2$coef[2] * chemical[i] + P2$coef[3] * chemical[i]^2
    y2_2 <- P2$coef[1] + P2$coef[2] * chemical[j] + P2$coef[3] * chemical[j]^2
    e2[count] <- (magnetic[i]-y2_1)^2+(magnetic[j]-y2_2)^2
    
    P3 <- lm(log(y)~x)
    y3_1 <- exp(P3$coef[1] + P3$coef[2] * chemical[i])
    y3_2 <- exp(P3$coef[1] + P3$coef[2] * chemical[j])
    e3[count] <- (magnetic[i]-y3_1)^2+(magnetic[j]-y3_2)^2
    
    P4 <- lm(log(y)~log(x))
    y4_1 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[i]))
    y4_2 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[j]))
    e4[count] <- (magnetic[i]-y4_1)^2+(magnetic[j]-y4_2)^2
  }

e = c(mean(e1)/2,mean(e2)/2,mean(e3)/2,mean(e4)/2)
matrix(e, nrow=1,
       dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))
detach(ironslag)

## -----------------------------------------------------------------------------

set.seed(12345)

# Count Five test
count5test = function(x, y) {
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y))
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
# Count Five test permutation
count5test_permutation = function(z) {

n = length(z)
x = z[1:(n/2)]
y = z[-(1:(n/2))]
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y)) 
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > 5))
}
permutation = function(z,R) {
  n = length(z)
  out = numeric(R)
  for (r in 1: R){
      p = sample(1:n ,n ,replace = FALSE)
      out[r] = count5test_permutation(z[p])
  }
  sum(out)/R
}              


n1 = 20
n2 = 50
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1e3

alphahat1 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
count5test(x, y)
}))
alphahat2 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
permutation(z,1000) 
})<0.05)

round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)


## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0,0,0)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
  mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
rbimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1), mu1, sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2), mu2, sd2)
  return(x)
}
for(i in 1:m){
  mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
  mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=10
n2=100
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic)
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  param <- c(0, 1, 1)
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
      if (u[i] <= (dskewlap(y, param = param) / dskewlap(x[i-1], param = param)))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
    }
  return(list(x=x, k=k))
}
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
param <- c(0, 1, 1)
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N)

refline <- qskewlap(c(.025, .975),param = param)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
  plot((rw)[,j] ,type="l",
    xlab=bquote(sigma == .(round(sigma[j],3))),
    ylab="X" , ylim=range(rw[,j]))
  abline(h=refline)
}

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

k <- 4    # four chains
x0 <- c(-10,-5,5,10)    # overdispersed initial values
N <- 10000    # length of chains
b <- 200    # burn-in length



X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(0.5,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (1000+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(1000+1):N], type="l", xlab="sigma=0.5", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(1,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (500+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x2 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(500+1):N], type="l", xlab="sigma=1", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(4,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x3 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=4", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(16,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x4 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=16", ylab="R_hat")
abline(h=1.2, lty=2)

c(x2,x3,x4)

## -----------------------------------------------------------------------------
k = c(4:25,100,500,1000)
S = function(a,k){
 ck = sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

f = function(a,k){S(a,k)-S(a,k-1)}
#curve(f(x),xlim = c(0,sqrt(k)))
a <- seq(0, 4, by=0.01)
plot(a, f(a, k[23]), lty=1, col=1, type="l", xlim=c(0, 4), xlab="a", ylab="f(a|k)", main="f(a) with different k")
lines(a, f(a, k[24]), xlim = c(0, 4), lty=2, col=2)
lines(a, f(a, k[25]), xlim = c(0, 4), lty=3, col=3)
legend("topright", legend=c("k=100", "k=500", "k=1000"), col=1:3,lty=1:3)
# So the lower and upper bound in function uniroot should be 1 and 2 respectively

solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

root = matrix(0,2,length(k))

for (i in 1:length(k)){
  root[2,i]=round(solve(k[i]),4)
}

root[1,] = k
rownames(root) = c('k','A(k)')
root


## -----------------------------------------------------------------------------

library(nloptr)
# Mle 
eval_f0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  return(sum(x)-0.999999)
}

opts = list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle = NULL
r = matrix(0,1,2)
r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j = 2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res = nloptr( x0=c(0.2,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=444,n.B=132,nOO=361,nAB=63 )
j = j+1
r = rbind(r,res$solution)
mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
#the result of EM algorithm
r 
#the max likelihood values
plot(-mle,type = 'l')


## -----------------------------------------------------------------------------

attach(mtcars)

formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
f3 = vector("list", length(formulas))
for (i in seq_along(formulas)){
  f3[[i]] = lm(formulas[[i]], data = mtcars)
}
f3
#2 lapply
la3 = lapply(formulas, function(x) lm(formula = x, data = mtcars))
la3


## -----------------------------------------------------------------------------

set.seed(123)
trials = replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])


## -----------------------------------------------------------------------------
datalist <- list(mtcars, faithful)
lapply(datalist, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE) return(simplify2array(out))
  unlist(out)
}
mylapply(datalist, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

cppFunction('List rw_Metropolis(double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0;
  int k = 0;
  NumericVector u = runif(N);
  for(int i = 1;i<N;i++){
    NumericVector y(1);
    y=rnorm(1,x[i-1],sigma);
    if(u[i] <= (exp(abs(x[i-1])-abs(y[0])))){
      x[i] = y[0];
    }
    else{
      x[i] = x[i-1];
      k=k+1;
    }
  }
  return List::create(Named("x") = x, Named("k") = k);
}')


set.seed(3000)

lap_f = function(x) exp(-abs(x))

rw.Metropolis = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}
N = 2000
sigma = 2
x0 = 25
rw1 = rw.Metropolis(sigma,x0,N)
rw2 = rw_Metropolis(sigma,x0,N)
plot(rw1$x,type = "l",xlab =bquote(sigma==2),ylab = "rw1",ylim = range(rw1$x))
plot(rw2$x,type = "l",xlab =bquote(sigma==2),ylab = "rw2",ylim = range(rw2$x))
#qqplot
set.seed(0)
rw1 = rw.Metropolis(sigma,x0,N)$x[-(1:100)]
rw2 = rw_Metropolis(sigma,x0,N)$x[-(1:100)]
qqplot(rw1,rw2)
abline(a=0,b=1,col='black')
#compare the time
set.seed(0)
(time = microbenchmark(rw1=rw.Metropolis(sigma,x0,N),rw2=rw_Metropolis(sigma,x0,N)))

