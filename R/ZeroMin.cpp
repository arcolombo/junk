src2<- 'using namespace Rcpp;

double mindof( NumericVector dof, NumericVector Indexes ) {

NumericVector DOF = dof[Indexes-1];
return floor(min(DOF));

} //sumsigma function


RCPP_MODULE(minDof) {

    function("mindof", &mindof);
} //exposing module
'

fx2 <- cxxfunction(signature(), plugin="Rcpp", include=src2)
 ss2 <- Module("minDof", getDynLib(fx2))



