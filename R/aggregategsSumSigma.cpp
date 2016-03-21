//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP aggregategsSumSigma(SEXP Means, SEXP SDs, SEXP DOFs, SEXP COLSs, SEXP geneSets) {
Rcpp::NumericVector Mean(Means);
Rcpp::NumericVector SD(SDs);
Rcpp::NumericVector DOF(DOFs);
Rcpp::CharacterVector COLS(COLSs);
Rcpp::CharacterVector geneSet(geneSets);

int n = Mean.size();
int m = SD.size();
int l = COLS.size();
int o = geneSet.size();
Rcpp::NumericVector idx(o);



//set these objects as armadillo objects 
arma::mat mean(Mean.begin(),Mean.size(),false);
arma::mat sd(SD.begin(),SD.size(),false);
arma::colvec dof(DOF.begin(),DOF.size(),false);
//arma::colvec cols(COLS.begin(),COLS.size(),false);
//arma::uvec geneset(geneSet.begin(),geneSet.size(),false);
//arma::mat sbase = zeros<arma::mat>(n,k);
//arma::mat sbase(n,k);
//arma::vec mn(n);
//arma::vec sd(n);
//arma::vec sda(n);
/*arma::mat x(Xr.begin(),n,k,false);
arma::mat y(Yr.begin(),n,k,false);
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
 */
return Rcpp::List::create( Rcpp::Named("geneSets") = Rcpp::wrap(idx));

}
