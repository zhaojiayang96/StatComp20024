#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @title count 
//' @description count the number of rejected
//' @param sigma the variance of the normal distribution 
//' @param x0 the initial value
//' @param N number of trials
//' @return a vector
//' @export
//[[Rcpp::export]]
NumericVector Metropolis(double sigma, double x0, int N){
  NumericVector x(N);
  x(0)=x0;
  NumericVector u = as<NumericVector>(runif(N));
  for(int i=1; i<N; i++){
    double y = as<double>(rnorm(1, x(i), sigma));
    if (u(i) <= (exp(-abs(y))/exp(-abs(x(i-1))))){
      x(i) = y;
    }
    else {
      x(i) = x(i-1);
    }
  }
  return(x);
}