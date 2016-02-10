#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector sigmasCpp(NumericMatrix x, NumericMatrix y, NumericVector v, NumericVector w, double z){
int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  
for (int i = 0; i < nrow ; i++){
    double total = 0;
    
for (int j =0; j < ncol ; j++){
     total+= pow((x(i,j)-y(i,j) - (v(i)-w(i))/z),2.0);
}
   out(i) = total/(z-1);
} 
return out);
} 
