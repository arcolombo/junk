#include <Rcpp.h>
using namespace Rcpp;

double sumsigma( NumericVector sd, NumericVector dof, NumericVector geneSets ) {
NumericVector Indexes = geneSets;
NumericVector SD = sd[Indexes-1];
NumericVector DOF = dof[Indexes-1];
return sqrt(sum((SD*SD)*(DOF/(DOF-2))));
} //sumsigma function

//[[Rcpp::export]]

List SumSigmaApply(NumericVector sd, NumericVector dof, List gList) {
return lapply(sd,dof,gList,sumsigma);

}
