#include<Rcpp.h>
using namespace Rcpp;
 class Zero {
public:
     Zero(NumericVector sd_, NumericVector dof_, List geneSets_) : sd(sd_) ,dof(dof_), geneSets(geneSets_){}

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


NumericVector mindof( NumericVector dof, List geneSets ) {
NumericVector outputDof(geneSets.size());
for(int i = 0 ; i < geneSets.size(); i++){
SEXP mm = geneSets[i];
NumericVector Indexes(mm);
NumericVector DOF = dof[Indexes-1];
outputDof[i] = floor(min(DOF));
  } 
return outputDof;

} //minDof



List sumSigma_minDof( NumericVector sd, NumericVector dof, List geneSets ) {
NumericVector output(geneSets.size());
NumericVector finalDof(geneSets.size());
for(int i = 0 ; i < geneSets.size() ; i++) {
SEXP nn = geneSets[i];
NumericVector Indexes(nn);
NumericVector SD = sd[Indexes-1];
NumericVector DOF = dof[Indexes-1];
output[i]= sqrt(sum((SD*SD)*(DOF/(DOF-2))));
finalDof[i] = floor(min(DOF));
   }

 return  Rcpp::List::create(Named("SumSigma") = output,
                                    Named("MinDof") = finalDof) 

} //sumsigma function



NumericVector sd, dof; 
List geneSets;
};

RCPP_MODULE(Sigma_Dof){
class_<Zero>("Zero")
.constructor<NumericVector,NumericVector,List>()
.field("sd", &Zero::sd)
.field("dof", &Zero::dof)
.field("geneSets", &Zero::geneSets)

.method("sumsigma", &Zero::sumsigma)
.method("mindof", &Zero::mindof)
.method("sumSigma_minDof", &Zero::sumSigma_minDof)
;
}



