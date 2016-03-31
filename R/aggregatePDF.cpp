//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include <math.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP aggregatePDF(SEXP MaxDiffs, SEXP geneSets, SEXP SDs, SEXP nPoints, SEXP DOFs) {

Rcpp::NumericVector MaxDiff(MaxDiffs);
Rcpp::List geneSet(geneSets);
Rcpp::NumericVector SD(SDs);
Rcpp::NumericVector nPoint(nPoints); //single constant vector lngth 1
Rcpp::NumericVector DOF(DOFs);


//note this function assumes that each input is not NA
//calculates the sd and mindof in order , names are assigned in R
//the geneSets non-empty are checked in R level

int n = SD.size();
int m = DOF.size();
int o = geneSet.size(); //we assume non-empty (reduce complexity)
int p = MaxDiff.size();   //o ==p  
Rcpp::NumericVector Norm(o);
Rcpp::NumericVector MaxDiffInt(o); //preallocating the size of list
//arma::vec Norm(o);
//arma::vec npoint(nPoint.begin(),nPoint.size(),false);
//arma::vec maxdiff(MaxDiff.begin(),MaxDiff.size(),false);
for ( int i =0; i<o ; i++){
 SEXP nn = geneSet[i];
 Rcpp::NumericVector index(nn);
int p = index.size(); //enteries in each geneSet list element

Rcpp::NumericVector sd = SD[index -1]; //converting to 0 based
Rcpp::NumericVector dof = DOF[index-1];

Norm[i]= (2*MaxDiff[i])/(nPoint[0]-1); //normalize

MaxDiffInt[i] = MaxDiff[i]/sd; 
Rcpp::NumericVector x1 = Rcpp::seq(-MaxDiffInt[i],MaxDiffInt[i]



} //for each list item

/*
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
*/
  
return Rcpp::List::create( Rcpp::Named("Norm") = Norm);
//                           Rcpp::Named("MinDof") = Rcpp::wrap(finalDof));
}
