//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include <math.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP aggregategsSumSigma( SEXP SDs, SEXP DOFs, SEXP geneSets) {
Rcpp::NumericVector SD(SDs);
Rcpp::NumericVector DOF(DOFs);
Rcpp::List geneSet(geneSets); 

//Rcpp::NumericVector geneSet(geneSets);
//note this function assumes that each input is not NA
//Fix Me: use sapply from sugar 
//for right now this holds for one gene set, need to test for multiple gene set case
int n = SD.size();
int m = DOF.size();
int o = geneSet.size(); //we assume non-empty (reduce complexity)
cout << " list size " << o << endl;
//Rcpp::NumericVector idx(n);

//for (int i =0; i<o; i++) {
 SEXP nn = geneSet[0];
 Rcpp::NumericVector index(nn);
 int p = index.size();
// arma::uvec idx = Rcpp::as<arma::uvec>(index);
 
//we subset before computing, and use the Rcpp index vector to avoid creating an unsigned vector as armadillo type
Rcpp::NumericVector sd = SD[index -1];
 Rcpp::NumericVector dof = DOF[index-1];
//Rcpp::NumericVector test= (sd*sd)*(dof/(dof-2));

//fast copy pointer address without data cache copy armadillo variables
arma::vec asd(sd.begin(),sd.size(),false);
arma::colvec adof(dof.begin(),dof.size(),false);


//arma::uvec idx = Rcpp::as<arma::uvec>(geneSet);
arma::vec test(p);
test =(asd%asd)%(adof/(adof-2));
arma::vec sumSigma(1); //there is a sumSigma for each geneSet
sumSigma = sqrt(sum(test)); //summing the subsets
arma::vec finalDof(1);
finalDof = floor(min(adof));

 
return Rcpp::List::create( Rcpp::Named("SumSigma") = Rcpp::wrap(sumSigma),
                           Rcpp::Named("MinDof") = Rcpp::wrap(finalDof));


//return Rcpp::List::create(Rcpp::Named("SumSigma") = sumSigma);
}
