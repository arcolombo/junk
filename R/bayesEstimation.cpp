//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP bayesEstimation(SEXP Xs, SEXP Ys, SEXP Zs, SEXP Vs, SEXP Qs) {
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericMatrix Yr(Ys);
Rcpp::NumericVector Zr(Zs);
Rcpp::NumericMatrix Vr(Vs);
Rcpp::NumericVector Qr(Qs);
int n = Xr.nrow(), k = Xr.ncol();
arma::mat x(Xr.begin(),n,k,false);
arma::mat y(Yr.begin(),n,k,false);
arma::vec z(Zr.begin(),Zr.size(),false);
arma::mat v(Vr.begin(),n,k,false);
arma::vec q(Qr.begin(),Qr.size(),false);

arma::mat sda(n,1);
arma::vec dof(n);
arma::vec sd(n);

sd  = (x/y);
sda = sd/(v%z);
dof = q;
Rcpp::NumericVector Sd = Rcpp::wrap(sd);
Rcpp::NumericVector Sda = Rcpp::wrap(sda);
Rcpp::NumericVector DOF = Rcpp::wrap(dof);
   Sd.names() = Rcpp::List(Xr.attr("dimnames"))[0];
   Sda.names() = Rcpp::List(Yr.attr("dimnames"))[0];
 
return Rcpp::List::create( Rcpp::Named("SD") = Sd,
                           Rcpp::Named("DOF") =  DOF,
                           Rcpp::Named("sd.alpha") = Sda);
}
