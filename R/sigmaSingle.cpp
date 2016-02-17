//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::export]]
extern "C" SEXP sigmaSingle(SEXP Bs, SEXP Ps, SEXP Ms) {
//declar
Rcpp::NumericMatrix Pr(Ps);
Rcpp::NumericMatrix Br(Bs);
Rcpp::NumericVector Mr(Ms);
double nb = Br.nrow(), kb = Br.ncol();
double np = Pr.nrow(), kp = Pr.ncol();
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
sbase.col(colm) = b.col(colm) - mnb;
spost.col(colm) = p.col(colm) - mnp;
}
sbase = sum(sbase%sbase,1)/(kb-1);
spost = sum(spost%spost,1)/(kp-1);

//arma::vec syoung = sbase + as_scalar(mvf);
//arma::vec sold = spost + as_scalar(mvf);

//dof = pow((syoung/kb + sold/kp),2)/((syoung%syoung/(as_scalar(kb*kb)*(as_scalar(kb-1))))+(sold%sold/(as_scalar(kp*kp)*(as_scalar(kp-1)))));

dof = pow(( (sbase + as_scalar(mvf))/kb + (spost + as_scalar(mvf))/kp),2)/(( ((sbase + as_scalar(mvf))%(sbase + as_scalar(mvf)) )/(kb*kb*(kb-1)))+( (spost + as_scalar(mvf))%(spost + as_scalar(mvf))/(kp*kp*(kp-1))));

sd  = sqrt((sbase/as_scalar(kb))+(spost/as_scalar(kp))+as_scalar(mvf));
sda = sd/sqrt(sbase/as_scalar(kb) + spost/as_scalar(kp));
Rcpp::NumericVector Mn = Rcpp::wrap(mn);
Rcpp::NumericVector Sd = Rcpp::wrap(sd);
Rcpp::NumericVector Sda = Rcpp::wrap(sda);
Rcpp::NumericVector DOF = Rcpp::wrap(dof);
   Mn.names() =  Rcpp::List(Br.attr("dimnames"))[0];
   Sd.names() = Rcpp::List(Br.attr("dimnames"))[0];
   Sda.names() = Rcpp::List(Br.attr("dimnames"))[0];
   DOF.names() = Rcpp::List(Br.attr("dimnames"))[0]; 
return Rcpp::List::create( Rcpp::Named("Mean") = Mn,
                           Rcpp::Named("SD") = Sd,
                           Rcpp::Named("DOF") = DOF,
                           Rcpp::Named("sd.alpha") = Sda); 
}
