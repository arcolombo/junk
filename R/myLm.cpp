#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List mylm(NumericVector yr, NumericMatrix Xr) {

    int n = Xr.nrow(), k = Xr.ncol();

    arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
    arma::colvec y(yr.begin(), yr.size(), false);

    arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
    arma::colvec resid = y - X*coef;            // residuals

    double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-k) );
                                                // std.error of estimate
    arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = coef,
        Rcpp::Named("stderr")       = stderrest
    ) ;

}
