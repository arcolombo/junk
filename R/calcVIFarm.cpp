//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP calcVIFarm(SEXP namesGrms, SEXP geneSets, SEXP rnEsets, SEXP esets, SEXP labels, SEXP sdAlphas) {
Rcpp::CharacterVector namesGrm(namesGrms);
Rcpp::NumericVector geneSet(geneSets);
Rcpp::CharacterVector rnEset(rnEsets);
Rcpp::NumericVector idx = geneSet-1; //0 based
Rcpp::CharacterVector GNames = namesGrm[idx];
Rcpp::NumericVector gs(0);
Rcpp::NumericVector grps(0); 
Rcpp::NumericMatrix eset(esets);
Rcpp::CharacterVector label(labels); //unique colnames(eset)
arma::mat covarMat;
Rcpp::NumericVector sdAlpha(sdAlphas); //must not be NULL check in R
arma::mat results;
arma::uvec agrps;
arma::vec vif;
int lengthGrps = label.size();
//cout << lengthGrps<<" " ;
int n = eset.nrow(), m = eset.ncol();
arma::mat est(eset.begin(),n,m,false);


//R is 1 based vectors , C++ is 0 based.  we subtract 1 from the R vectors to make 0 based, and extract rownames.  then we make gs 1 based when pushing back to R; this matches R compus


for(int i =0; i < rnEset.size();i++){
    for( int j=0; j<GNames.size();j++){
        if(rnEset(i)==GNames(j)){
         gs.push_back( i +1); // 1 based  
         
        } 
    }
 }
    if (gs.size() <2){
      cout << "GeneSet contains one of zero overlapping genes.  NAs produced";
      Rcpp::CharacterVector gs = "NA";
     }
  
//arma::vec ags(gs.begin(),gs.size(),false);
// ags = ags -1;//0 based
arma::uvec ags = Rcpp::as<arma::uvec>(gs) - 1; //0 based

//find groups to split the eset by labels using only unique labeling, where column names of eset is in terms of labels
Rcpp::CharacterVector colN = Rcpp::List(eset.attr("dimnames"))[1];
 covarMat = covarMat.zeros(gs.size(), gs.size()); // the covariance matrix is square by nrow X nrow by definition of the inner product of a matrix with itself
 
//index of columns matching label type
   for(int k =0; k<label.size();k++) { 
    Rcpp::NumericVector grps(0);
   for( int j = 0; j< eset.ncol();j++){
      if(colN(j) == label(k)){
        grps.push_back(j+1);  // 1 based
       
          }
      }
 //     arma::uvec atmp(tmp.begin(),tmp.size(),false);
   //         atmp = atmp - 1 ; //0 based
    // int tmpSize = tmp.size();
     arma::uvec agrps = Rcpp::as<arma::uvec>(grps) -1; 
    // cout << agrps;
      covarMat += cov(est.submat(ags,agrps).t()) * as_scalar(grps.size()-1) ;
   } //for each label 
 
    covarMat = covarMat/as_scalar(m - lengthGrps);

//FIX ME: check the sd.alpha computations
  
 //multiply matrix by the sd.alpha vectors
//Rcpp::CharacterVector a = rnEset[gs-1];
sdAlpha = sdAlpha[gs-1];
 //cout << covarMat;
arma::vec sdalpha(sdAlpha.begin(),sdAlpha.size(),false);

  

// FIX ME : fix the  covar.mat = t(covar.mat*a)*a computations
int row = size(covarMat)[0];
int col = size(covarMat)[1];
//cout << " row "<<row << " col "<<col; 
results = results.zeros(row,col); //storage object

 for(int i=0; i<row;i++){
    for(int j =0; j<col;j++){
     covarMat(i,j) = covarMat(i,j)*as_scalar(sdalpha(i));
   }
} 
 
 for(int i=0; i<row;i++){
    for(int j =0; j<col;j++){
     results(i,j) = covarMat(j,i)*as_scalar(sdalpha(i));
   }
}

 

//cout<<accu(results)<<" "<<accu(results.diag()) ;
vif = accu(results)/accu(results.diag());

//FIX ME: just return the vif
return Rcpp::List::create(Rcpp::Named("GNames") = GNames,
                         // Rcpp::Named("gs.i") = Rcpp::wrap(ags), 0 based
                          Rcpp::Named("gs.i") = gs,
                          Rcpp::Named("covar.mat") = Rcpp::wrap(covarMat),
                          Rcpp::Named("final.covar.mat") = Rcpp::wrap(results),
                          Rcpp::Named("a") = Rcpp::wrap(sdAlpha),
                          Rcpp::Named("vif") = Rcpp::wrap(vif));
                         
                           
                           
                           
}