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


//Rcpp::List geneSet(geneSets); I had trouble writing from a List perspective
Rcpp::NumericVector geneSet(geneSets);
//note this function assumes that each input is not NA
//Fix Me: use sapply from sugar 

//for right now this holds for one gene set, need to test for multiple gene set case
int n = SD.size();
int m = DOF.size();
int o = geneSet.size(); //we assume non-empty (reduce complexity)

//set these objects as armadillo objects 
arma::vec sd(SD.begin(),SD.size(),false);
arma::colvec dof(DOF.begin(),DOF.size(),false);
arma::uvec idx = Rcpp::as<arma::uvec>(geneSet);
arma::vec test(n);
test =(sd%sd)%(dof/(dof-2));
arma::vec sumSigma(1); //there is a sumSigma for each geneSet
sumSigma = sqrt(sum(test.elem(idx-1)));
arma::vec finalDof(1);
finalDof = floor(min(dof.elem(idx-1)));

 
return Rcpp::List::create( Rcpp::Named("SumSigma") = Rcpp::wrap(sumSigma),
                           Rcpp::Named("MinDof") = Rcpp::wrap(finalDof));

}
