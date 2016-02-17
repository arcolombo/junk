//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP dof(SEXP SOlds, SEXP SYous, SEXP NOlds, SEXP NYoungs) {
Rcpp::NumericVector SOldr(SOlds);
Rcpp::NumericVector SYour(SYous);
Rcpp::NumericVector NOldr(NOlds);
Rcpp::NumericVector NYoungr(NYoungs);
arma::vec sold(SOldr.begin(),SOldr.size(),false);
arma::vec syoung(SYour.begin(),SYour.size(),false);
arma::vec Nol(NOldr.begin(),NOldr.size(),false);
arma::vec Nyo(NYoungr.begin(),NYoungr.size(),false);
arma::vec dof(SYour.size());
int co = as_scalar(Nol);
int cy = as_scalar(Nyo);
dof = pow((syoung/cy + sold/co),2)/( (syoung%syoung/((cy%cy)%(cy-1))) + (sold%sold/((co%co)%(co-1))) );
Rcpp::NumericVector DOF = Rcpp::wrap(dof);
return DOF;
}
