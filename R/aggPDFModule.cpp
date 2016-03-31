src<- 'using namespace Rcpp;

NumericVector sumsigma( NumericVector sd, NumericVector dof, List geneSets ) {
NumericVector output(geneSets.size());
for(int i = 0 ; i < geneSets.size() ; i++) {
SEXP nn = geneSets[i];
NumericVector Indexes(nn);
NumericVector SD = sd[Indexes-1];
NumericVector DOF = dof[Indexes-1];
output[i]= sqrt(sum((SD*SD)*(DOF/(DOF-2))));
   }
 return output;
} //sumsigma function

//List sumSigmaApply( List gList){

// return lapply(gList,sumsigma(sd, dof, );
//}

RCPP_MODULE(modsigma) {

    function("sumsigma", &sumsigma);
} //exposing module
'

fx <- cxxfunction(signature(), plugin="Rcpp", include=src)
 ss <- Module("modsigma", getDynLib(fx))


