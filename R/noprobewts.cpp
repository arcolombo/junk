//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP noprobewts(SEXP Ms, SEXP coeffs, SEXP resids, SEXP effts, SEXP rnks, SEXP ftvls, SEXP assgs, SEXP qrs, SEXP dfrs) {
Rcpp::NumericMatrix M(Ms);
Rcpp::NumericMatrix Coeff(coeffs);
Rcpp::NumericMatrix Resid(resids);
Rcpp::NumericMatrix Efft(effts);
Rcpp::NumericVector Rnk(rnks);
Rcpp::NumericMatrix Ftvl(ftvls);
Rcpp::NumericVector Assg(assgs);
Rcpp::NumericMatrix Qr(qrs);
Rcpp::NumericVector Dfr(dfrs);

int ngenes = Xr.nrow(), k = Xr.ncol();
arma::vec dfr(Dfr.begin(),Dfr.size(),false);
if(dfr>0){
    //here we assume that fit$effects is a matrix. FIX ME: will need another script to design if fit$effects is not a matrix.

}
/*arma::mat y(Yr.begin(),n,k,false);
arma::mat z(Zr.begin(),1,1,false);
arma::mat r(Vr.begin(),1,1,false);
//arma::mat sbase = zeros<arma::mat>(n,k);
arma::mat sbase(n,k);
arma::vec mn(n);
arma::vec sd(n);
arma::vec sda(n);
mn = (sum(y,1)-sum(x,1))/as_scalar(z);
for (int colm=0;colm<k;colm++){
sbase.col(colm) = x.col(colm) - y.col(colm) + mn;
}
sbase = sum(sbase%sbase,1)/(as_scalar(z)-1);
sd  = sqrt((sbase/as_scalar(z)) + as_scalar(r));
sda = sd/sqrt(sbase/as_scalar(z));
Rcpp::NumericVector Mn = Rcpp::wrap(mn);
Rcpp::NumericVector Sd = Rcpp::wrap(sd);
Rcpp::NumericVector Sda = Rcpp::wrap(sda);
   Mn.names() =  Rcpp::List(Xr.attr("dimnames"))[0];
   Sd.names() = Rcpp::List(Xr.attr("dimnames"))[0];
   Sda.names() = Rcpp::List(Xr.attr("dimnames"))[0];
 
return Rcpp::List::create( Rcpp::Named("Mean") = Mn,
                           Rcpp::Named("SD") = Sd,
                           Rcpp::Named("DOF") =  Rcpp::wrap(z),
                           Rcpp::Named("SDAlpha") = Sda); */
}
