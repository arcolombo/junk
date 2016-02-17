//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP sigmaSingle(SEXP Bs, SEXP Ps, SEXP Ms) {
//declar
Rcpp::NumericMatrix Pr(Ps);
Rcpp::NumericMatrix Br(Bs);
Rcpp::NumericVector Mr(Ms);
int nb = Br.nrow(), kb = Br.ncol();
int np = Pr.nrow(), kp = Pr.ncol();
arma::mat b(Br.begin(),nb,kb,false);
arma::mat p(Pr.begin(),nb,kb,false);
arma::mat mvf(Mr.begin(),1,1,false);

arma::mat sbase(nb,kb);
arma::mat spost(np,kp);
arma::vec mnb(nb);
arma::vec mnp(np);
arma::vec mn(nb);
arma::vec sd(nb);
arma::vec sda(nb);
arma::vec dof(nb);

//calc
mnp = sum(p,1)/kb;
mnb = sum(b,1)/kb;
mn = mnp - mnb;
for (int colm=0;colm<kb;colm++){
sbase.col(colm) = b(colm) - mnb;
spost.col(colm) = p(colm) - mnp;
}
sbase = sum(sbase%sbase,1)/(kb-1);
spost = sum(spost%spost,1)/(kp-1);

arma::vec syoung = sbase + as_scalar(mvf);
arma::vec sold = spost + as_scalar(mvf);


dof = pow((syoung/kb + sold/kp),2)/( (syoung%syoung/((kb%kb)%(kb-1))) + (sold%sold/((kp%kp)%(kp-1))) );

return Rcpp::List::create( Rcpp::Named("DOF") = Rcpp::wrap(dof),
                           Rcpp::Named("syoung") = Rcpp::wrap(syoung),
                           Rcpp::Named("sold") = Rcpp::wrap(sold));

/*sd  = sqrt((sbase/as_scalar(kb)) + as_scalar(mvf));
sda = sd/sqrt(sbase/as_scalar(kb));
Rcpp::NumericVector Mn = Rcpp::wrap(mn);
Rcpp::NumericVector Sd = Rcpp::wrap(sd);
Rcpp::NumericVector Sda = Rcpp::wrap(sda);
   Mn.names() =  Rcpp::List(Br.attr("dimnames"))[0];
   Sd.names() = Rcpp::List(Br.attr("dimnames"))[0];
   Sda.names() = Rcpp::List(Br.attr("dimnames"))[0];
 
return Rcpp::List::create( Rcpp::Named("Mean") = Mn,
                           Rcpp::Named("SD") = Sd,
                           Rcpp::Named("DOF") =  Rcpp::wrap(kb),
                           Rcpp::Named("SDAlpha") = Sda); */
}
