//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP notbayesEstimation(SEXP Xs, SEXP Ys, SEXP Zs, SEXP Vs) {
Rcpp::NumericVector Xr(Xs);
Rcpp::NumericMatrix Yr(Ys);
Rcpp::NumericMatrix Zr(Zs);
Rcpp::NumericMatrix Vr(Vs);
int n = Yr.nrow(), k = Yr.ncol();
int l = Vr.nrow(), m = Vr.ncol();
arma::mat x(Xr.begin(),n,k,false); //fit2$sigma
arma::mat y(Yr.begin(),n,k,false); //fit2b$stdev.unscaled
arma::mat z(Zr.begin(),1,1,false); //min.variance.factor
arma::mat v(Vr.begin(),l,m,false); //fit2$df.residual


arma::mat sda(n,1);
arma::vec dof(n);
arma::vec sd(n);

sd  = sqrt(  (y%x)%(y%x) + as_scalar(z));
sda = sd/(x%y);
dof = v;
Rcpp::NumericVector Sd = Rcpp::wrap(sd);
Rcpp::NumericVector Sda = Rcpp::wrap(sda);
Rcpp::NumericVector DOF = Rcpp::wrap(dof);
   Sd.names() = Rcpp::List(Yr.attr("dimnames"))[0];
   Sda.names() = Rcpp::List(Yr.attr("dimnames"))[0];
 
return Rcpp::List::create( Rcpp::Named("SD") = Sd,
                           Rcpp::Named("DOF") =  DOF,
                           Rcpp::Named("sd.alpha") = Sda);
}
