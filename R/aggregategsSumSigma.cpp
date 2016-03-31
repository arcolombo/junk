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


//note this function assumes that each input is not NA
//calculates the sd and mindof in order , names are assigned in R

int n = SD.size();
int m = DOF.size();
int o = geneSet.size(); //we assume non-empty (reduce complexity)
arma::vec sumSigma(o); //there is a sumSigma for each geneSet
arma::vec finalDof(o);


//need to run the sapply function
for ( int i=0 ; i < o ; i++) { //running a for loop
 SEXP nn = geneSet[i];
 Rcpp::NumericVector index(nn);
 int p = index.size();
arma::vec test(p);
 //we subset before computing, and use the Rcpp index vector to avoid creating an unsigned vector as armadillo type
Rcpp::NumericVector sd = SD[index -1]; //converting to 0 based
Rcpp::NumericVector dof = DOF[index-1];
//fast copy pointer address without data cache copy armadillo variables
arma::vec asd(sd.begin(),sd.size(),false);
arma::colvec adof(dof.begin(),dof.size(),false);
test =(asd%asd)%(adof/(adof-2));
sumSigma(i) = sqrt(sum(test)); //summing the subsets
finalDof(i) = floor(min(adof)); 
}

  
return Rcpp::List::create( Rcpp::Named("SumSigma") = Rcpp::wrap(sumSigma),
                           Rcpp::Named("MinDof") = Rcpp::wrap(finalDof));
}
