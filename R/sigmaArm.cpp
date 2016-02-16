//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP sigmaArm(SEXP Xs, SEXP Ys, SEXP Zs, SEXP Vs) {
Rcpp::NumericMatrix Yr(Ys);
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector Zr(Zs);
Rcpp::NumericVector Vr(Vs);
int n = Xr.nrow(), k = Xr.ncol();
arma::mat x(Xr.begin(),n,k,false);
arma::mat y(Yr.begin(),n,k,false);
arma::mat z(Zr.begin(),1,1,false);
arma::mat r(Vr.begin(),1,1,false);
//arma::mat sbase = zeros<arma::mat>(n,k);
arma::mat sbase(n,k);
arma::vec mn(n);
arma::vec sd(n);
mn = (sum(y,1)-sum(x,1))/as_scalar(z);
for (int colm=0;colm<k;colm++){
sbase.col(colm) = x.col(colm) - y.col(colm) + mn;
}
sbase = sum(sbase%sbase,1)/(as_scalar(z)-1);
sd  = sqrt((sbase/as_scalar(z)) + as_scalar(r));

return Rcpp::List::create( Rcpp::Named("Mean") = Rcpp::wrap(mn),
                           Rcpp::Named("SD") = Rcpp::wrap(sd),
                           Rcpp::Named("DOF") = Rcpp::wrap(z));
}
